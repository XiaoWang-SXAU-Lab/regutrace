//! ReguTrace core: local TF–target evidence construction from genome, GFF3,
//! PWM motifs and user-provided evidence lists. Outputs evidence records with
//! provenance, not a definitive GRN.

pub mod fasta;
pub mod gff;
pub mod io;
pub mod promoter;
pub mod pwm;
pub mod scan;
pub mod store;
