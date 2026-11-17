use std::collections::{BTreeMap, HashMap, HashSet};

use anyhow::{ensure, Result};
use serde::Deserialize;

use crate::hmmer::DomainHit;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RulePack {
    pub version: String,
    /// HMM NAME to domain bit-score threshold, or null for the E-value fallback.
    pub domains: BTreeMap<String, Option<f64>>,
    /// Optional per-model full-sequence E-value thresholds. These take
    /// precedence over the bit-score/fallback decision for that model.
    #[serde(default)]
    pub domain_evalues: BTreeMap<String, f64>,
    pub fallback_evalue: f64,
    /// Ordered groups. By default, the first matching rule in each group wins.
    /// With `best_match`, the single matching rule with the smallest supporting
    /// full-sequence E-value wins across all groups (AnimalTFDB-style).
    #[serde(default)]
    pub best_match: bool,
    pub groups: Vec<Vec<FamilyRule>>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FamilyRule {
    pub family: String,
    pub subfamily: String,
    #[serde(default)]
    pub required: BTreeMap<String, usize>,
    #[serde(default)]
    pub any: Vec<String>,
    #[serde(default)]
    pub maximum: BTreeMap<String, usize>,
    #[serde(default)]
    pub only_listed_domains: bool,
}

pub struct ProteinClassification<'a> {
    pub protein_id: &'a str,
    pub families: Vec<&'a FamilyRule>,
    pub domains: Vec<&'a DomainHit>,
}

pub fn classify<'a>(
    pack: &'a RulePack,
    hits: &'a [DomainHit],
) -> Result<(Vec<ProteinClassification<'a>>, Vec<&'a str>)> {
    pack.validate()?;
    let mut proteins = BTreeMap::<&str, Vec<&DomainHit>>::new();
    let mut unknown = HashSet::new();
    for hit in hits {
        let Some(cutoff) = pack.domains.get(&hit.model_name) else {
            unknown.insert(hit.model_name.as_str());
            continue;
        };
        let passed = if let Some(evalue) = pack.domain_evalues.get(&hit.model_name) {
            hit.full_evalue < *evalue
        } else {
            match cutoff {
                Some(score) => hit.domain_score >= *score,
                None => hit.full_evalue < pack.fallback_evalue,
            }
        };
        if passed {
            proteins.entry(&hit.protein_id).or_default().push(hit);
        }
    }
    let mut classified = Vec::new();
    for (protein_id, domains) in proteins {
        let mut counts = HashMap::new();
        for hit in &domains {
            *counts.entry(hit.model_name.as_str()).or_insert(0) += 1;
        }
        let families: Vec<_> = if pack.best_match {
            let mut best: Option<(&FamilyRule, f64)> = None;
            for rule in pack.groups.iter().flatten() {
                if !rule.matches(&counts) {
                    continue;
                }
                let evalue = rule.support_evalue(&domains);
                if best.is_none_or(|(_, best_evalue)| evalue < best_evalue) {
                    best = Some((rule, evalue));
                }
            }
            best.map(|(rule, _)| rule).into_iter().collect()
        } else {
            pack.groups
                .iter()
                .filter_map(|group| group.iter().find(|rule| rule.matches(&counts)))
                .collect()
        };
        if !families.is_empty() {
            classified.push(ProteinClassification { protein_id, families, domains });
        }
    }
    let mut unknown: Vec<_> = unknown.into_iter().collect();
    unknown.sort_unstable();
    Ok((classified, unknown))
}

impl RulePack {
    fn validate(&self) -> Result<()> {
        ensure!(!self.version.trim().is_empty(), "missing rule pack version");
        ensure!(
            self.fallback_evalue.is_finite() && self.fallback_evalue > 0.0,
            "fallback E-value must be positive and finite"
        );
        ensure!(!self.domains.is_empty() && !self.groups.is_empty(), "empty rule pack");
        for (name, score) in &self.domains {
            ensure!(!name.trim().is_empty(), "empty model name");
            ensure!(score.is_none_or(f64::is_finite), "invalid threshold for {name}");
        }
        for (name, evalue) in &self.domain_evalues {
            ensure!(self.domains.contains_key(name), "unknown model in domain_evalues: {name}");
            ensure!(
                evalue.is_finite() && *evalue > 0.0,
                "invalid E-value threshold for {name}"
            );
        }
        for group in &self.groups {
            ensure!(!group.is_empty(), "empty rule group");
            for rule in group {
                ensure!(
                    !rule.family.trim().is_empty() && !rule.subfamily.trim().is_empty(),
                    "family and subfamily are required"
                );
                ensure!(
                    !rule.required.is_empty() || !rule.any.is_empty(),
                    "rule {} needs a positive domain condition", rule.subfamily
                );
                for (domain, count) in &rule.required {
                    ensure!(*count > 0, "invalid minimum count for {domain}");
                    ensure!(
                        rule.maximum.get(domain).is_none_or(|max| max >= count),
                        "maximum below minimum for {domain}"
                    );
                }
                for domain in rule.required.keys().chain(rule.any.iter()).chain(rule.maximum.keys()) {
                    ensure!(self.domains.contains_key(domain), "unknown model in rule: {domain}");
                }
            }
        }
        Ok(())
    }
}

impl FamilyRule {
    /// Best full-sequence E-value among the rule's required (defining) models.
    /// Shared context domains in `any` and veto domains in `maximum` do not
    /// contribute, so a generic backbone domain cannot win a tie for a
    /// subfamily rule.
    fn support_evalue(&self, domains: &[&DomainHit]) -> f64 {
        domains
            .iter()
            .filter(|hit| self.required.contains_key(hit.model_name.as_str()))
            .map(|hit| hit.full_evalue)
            .fold(f64::INFINITY, f64::min)
    }

    fn matches(&self, counts: &HashMap<&str, usize>) -> bool {
        let count = |domain: &String| counts.get(domain.as_str()).copied().unwrap_or(0);
        self.required.iter().all(|(domain, min)| count(domain) >= *min)
            && (self.any.is_empty() || self.any.iter().any(|domain| count(domain) > 0))
            && self.maximum.iter().all(|(domain, max)| count(domain) <= *max)
            && (!self.only_listed_domains || counts.keys().all(|domain| {
                self.required.contains_key(*domain) || self.any.iter().any(|listed| listed == *domain)
            }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(full_evalue: f64, domain_score: f64) -> DomainHit {
        DomainHit {
            protein_id: "P1".to_owned(),
            model_name: "Forkhead".to_owned(),
            model_accession: Some("PF00250.26".to_owned()),
            model_length: 88,
            full_evalue,
            domain_score,
            conditional_evalue: 1e-10,
            independent_evalue: 1e-10,
            model_start: 2,
            model_end: 88,
            protein_start: 10,
            protein_end: 96,
        }
    }

    fn pack() -> RulePack {
        RulePack {
            version: "test".to_owned(),
            domains: BTreeMap::from([("Forkhead".to_owned(), Some(999.0))]),
            domain_evalues: BTreeMap::from([("Forkhead".to_owned(), 1e-2)]),
            fallback_evalue: 1e-4,
            best_match: false,
            groups: vec![vec![FamilyRule {
                family: "Fork_head".to_owned(),
                subfamily: "Fork_head".to_owned(),
                required: BTreeMap::from([("Forkhead".to_owned(), 1)]),
                any: Vec::new(),
                maximum: BTreeMap::new(),
                only_listed_domains: false,
            }]],
        }
    }

    #[test]
    fn per_model_evalue_overrides_bit_score_threshold() {
        let pack = pack();
        let passing_hits = [hit(1e-3, 1.0)];
        let passed = classify(&pack, &passing_hits).unwrap().0;
        assert_eq!(passed.len(), 1);
        let failing_hits = [hit(1e-1, 5000.0)];
        let failed = classify(&pack, &failing_hits).unwrap().0;
        assert!(failed.is_empty());
    }

    #[test]
    fn best_match_chooses_smallest_supporting_evalue() {
        let mut pack = pack();
        pack.best_match = true;
        pack.domains.insert("Other".to_owned(), None);
        pack.domain_evalues.insert("Other".to_owned(), 1e-2);
        pack.groups.push(vec![FamilyRule {
            family: "Other".to_owned(),
            subfamily: "Other".to_owned(),
            required: BTreeMap::from([("Other".to_owned(), 1)]),
            any: Vec::new(),
            maximum: BTreeMap::new(),
            only_listed_domains: false,
        }]);
        let mut forkhead = hit(1e-3, 50.0);
        forkhead.model_name = "Forkhead".to_owned();
        let mut other = hit(1e-4, 1.0);
        other.model_name = "Other".to_owned();
        let hits = [forkhead, other];
        let classified = classify(&pack, &hits).unwrap().0;
        assert_eq!(classified.len(), 1);
        assert_eq!(classified[0].families[0].subfamily, "Other");
    }

    #[test]
    fn rejects_evalue_for_unknown_model() {
        let mut pack = pack();
        pack.domain_evalues.insert("NotInDomains".to_owned(), 1e-5);
        let err = match classify(&pack, &[]) {
            Ok(_) => panic!("expected an invalid rule pack"),
            Err(err) => err,
        };
        let err = format!("{err:#}");
        assert!(err.contains("unknown model in domain_evalues"), "{err}");
    }
}
