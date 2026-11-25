#!/usr/bin/env python3
"""Dinucleotide-preserving promoter null on a representative subset."""
import argparse
import collections
import csv
import json
import random


def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split("|")[0]
                buf = []
            else:
                buf.append(line.strip().upper())
        if name:
            seqs[name] = "".join(buf)
    return seqs


def dinuc_shuffle(seq, rng):
    # Hierholzer Eulerian trail over the dinucleotide multigraph. This exactly
    # preserves every dinucleotide count while shuffling their order.
    edges = [seq[i:i+2] for i in range(len(seq)-1)]
    graph = collections.defaultdict(list)
    for a, b in edges:
        graph[a].append(b)
    indeg = collections.Counter(b for _, b in edges)
    outdeg = collections.Counter(a for a, _ in edges)
    starts = [n for n in graph if outdeg[n] > indeg[n]]
    start = starts[0] if starts else edges[0][0]
    stack = [start]
    circuit = []
    while stack:
        v = stack[-1]
        if graph.get(v):
            i = rng.randrange(len(graph[v]))
            u = graph[v].pop(i)
            stack.append(u)
        else:
            circuit.append(stack.pop())
    circuit.reverse()
    if len(circuit) != len(edges) + 1:
        raise RuntimeError("dinucleotide graph is disconnected")
    return "".join(n[0] for n in circuit[:-1]) + circuit[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--promoters", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--limit", type=int, default=200)
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args()
    seqs = read_fasta(args.promoters)
    rng = random.Random(args.seed)
    selected = [(g, s) for g, s in seqs.items() if len(s) == 1000][:args.limit]
    with open(args.out_fasta, "w") as f:
        for gene, seq in selected:
            f.write(f">{gene}\n")
            shuffled = dinuc_shuffle(seq, rng)
            for i in range(0, len(shuffled), 80):
                f.write(shuffled[i:i+80] + "\n")
    json.dump({"promoters": len(selected), "seed": args.seed,
               "rule": "dinucleotide counts preserved; first base sampled from dinucleotide starts"},
              open(args.out_json, "w"), indent=2)
    print(f"wrote {len(selected)} shuffled promoters")


if __name__ == "__main__":
    main()
