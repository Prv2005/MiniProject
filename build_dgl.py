#!/usr/bin/env python3
"""
De Bruijn Graph builder for TP53 mRNA analysis.

Pipeline:
  1. Parse reads.fastq  -> ignore quality scores, extract sequences only
  2. Build De Bruijn graph entirely from reads (k-mer counting + adjacency via 4-base extension)
  3. Label edges using reference.fasta (y=1 if (k+1)-mer appears in reference, else y=0)
  4. Compute rich node/edge features + Z-score normalisation inside DGL
  5. Save DGL graph as graph.bin
"""

import os
import re
import dgl
import torch
import numpy as np
import networkx as nx
from collections import Counter

os.makedirs("output", exist_ok=True)


# ============================================================
# 1.  PARSERS  (no BioPython dependency, handles raw FASTQ/FASTA)
# ============================================================

def parse_fasta(path: str) -> list[str]:
    """Return a list of sequence strings from a FASTA file."""
    sequences = []
    current_seq_parts = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_seq_parts:
                    sequences.append("".join(current_seq_parts).upper())
                    current_seq_parts = []
            else:
                current_seq_parts.append(line)
        if current_seq_parts:
            sequences.append("".join(current_seq_parts).upper())
    return sequences


def parse_fastq_sequences_only(path: str) -> list[str]:
    """
    Robust FASTQ parser that extracts sequence strings only.
    Quality scores are completely ignored.

    Strategy — state machine with 3 states:
      HEADER   : looking for a line starting with '@'
      SEQUENCE : collecting sequence lines until we hit a '+' separator
      QUALITY  : skipping quality lines until we see the next '@' header

    This correctly handles:
      - Standard 4-line FASTQ records
      - Multi-line sequences (sequence wrapped across 2+ lines)
      - Multi-line quality blocks
      - '.fq' and '.fastq' extensions
    """
    sequences = []
    state = "HEADER"
    current_seq_parts: list[str] = []

    with open(path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip()
            if not line:
                continue                          # skip blank lines everywhere

            if state == "HEADER":
                if line.startswith("@"):
                    current_seq_parts = []
                    state = "SEQUENCE"
                # anything else before the first header is ignored

            elif state == "SEQUENCE":
                if line.startswith("+"):          # separator line → quality block starts
                    if current_seq_parts:
                        sequences.append("".join(current_seq_parts).upper())
                    state = "QUALITY"
                else:
                    current_seq_parts.append(line) # accumulate (handles multi-line seq)

            elif state == "QUALITY":
                # Quality lines are completely ignored.
                # The only thing we watch for is the next record header.
                # A '@' at the start of a quality line is valid ASCII-33 quality,
                # BUT in practice we can distinguish it because after consuming
                # the expected number of quality characters we'd see the next '@'.
                # Simplest safe heuristic: once we hit '@' again, treat it as a header.
                if line.startswith("@"):
                    current_seq_parts = []
                    state = "SEQUENCE"
                # otherwise just skip — quality data discarded

    return sequences


# ============================================================
# 2.  K-MER COUNTING
# ============================================================

def count_kmers(sequences: list[str], k: int) -> Counter:
    """Count every k-mer across all sequences."""
    kmer_counts: Counter = Counter()
    valid = re.compile(r'^[ACGT]+$')          # skip k-mers with N or other chars
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer = seq[i: i + k]
            if valid.match(kmer):
                kmer_counts[kmer] += 1
    return kmer_counts


# ============================================================
# 3.  GRAPH BUILDER  (adjacency via 4-base extension, like Code-2)
# ============================================================

class DeBruijnGraphBuilder:
    BASES = ("A", "T", "C", "G")

    def __init__(self, kmer_counts: Counter, reference_sequences: list[str], k: int):
        self.k = k
        self.kmer_counts = kmer_counts
        self.reference_sequences = reference_sequences

        # Assign integer index to every k-mer (node)
        self.kmer_to_idx: dict[str, int] = {
            kmer: idx for idx, kmer in enumerate(kmer_counts.keys())
        }

        self.nodes: dict[int, dict] = {}   # idx -> {'kmer': str, 'abundance': int}
        self.edges: dict[tuple, dict] = {} # (src_idx, dst_idx) -> feature dict

    # ----------------------------------------------------------
    # 3a. Reference-based edge labelling
    # ----------------------------------------------------------
    def _label_edge(self, edge_seq: str) -> int:
        """Return 1 if the (k+1)-mer edge_seq appears in any reference sequence."""
        for ref in self.reference_sequences:
            if edge_seq in ref:
                return 1
        return 0

    # ----------------------------------------------------------
    # 3b. Build nodes dict
    # ----------------------------------------------------------
    def _build_nodes(self):
        print(f"  Building {len(self.kmer_to_idx)} nodes …")
        for kmer, idx in self.kmer_to_idx.items():
            self.nodes[idx] = {
                "kmer":      kmer,
                "abundance": self.kmer_counts[kmer],
            }

    # ----------------------------------------------------------
    # 3c. Build edges via 4-base extension (identical strategy to Code-2)
    # ----------------------------------------------------------
    def _build_edges(self):
        print("  Building edges via 4-base extension …")
        kmer_set = self.kmer_to_idx           # fast membership test

        for kmer, src_idx in kmer_set.items():
            suffix = kmer[1:]                 # (k-1)-mer suffix of current node

            for base in self.BASES:
                neighbor = suffix + base      # candidate successor k-mer
                if neighbor in kmer_set:
                    dst_idx  = kmer_set[neighbor]
                    edge_seq = kmer + base    # (k+1)-mer spanning the edge

                    # ---- Edge features (same philosophy as Code-2) ----
                    src_abund = self.kmer_counts[kmer]
                    dst_abund = self.kmer_counts[neighbor]

                    abundance_avg     = (src_abund + dst_abund) / 2.0
                    occurrence_sim    = abs(src_abund - dst_abund)       # lower = more similar
                    label             = self._label_edge(edge_seq)

                    self.edges[(src_idx, dst_idx)] = {
                        "edge_seq":        edge_seq,
                        "abundance_avg":   float(abundance_avg),
                        "occurrence_sim":  float(occurrence_sim),
                        "label":           float(label),
                    }

        print(f"  Found {len(self.edges)} edges.")

    # ----------------------------------------------------------
    # 3d. Build DGL graph with Z-score normalisation (PyTorch, like Code-2)
    # ----------------------------------------------------------
    def _build_dgl_graph(self) -> dgl.DGLGraph:
        print("  Assembling DGL graph …")

        # --- NetworkX directed graph as intermediate ---
        DiGraph = nx.DiGraph()

        for idx, attrs in self.nodes.items():
            DiGraph.add_node(idx, x=float(attrs["abundance"]))

        for (src, dst), eattrs in self.edges.items():
            DiGraph.add_edge(
                src, dst,
                e=[eattrs["abundance_avg"], eattrs["occurrence_sim"]],
                y=eattrs["label"],
            )

        g = dgl.from_networkx(DiGraph, node_attrs=["x"], edge_attrs=["e", "y"])

        # ---- Z-score normalisation: node features ----
        node_feats = g.ndata["x"].float().unsqueeze(1)  # shape [N, 1]
        n_mean = node_feats.mean(dim=0)
        n_std  = node_feats.std(dim=0).clamp(min=1e-9)
        g.ndata["x"] = torch.round((node_feats - n_mean) / n_std * 10000) / 10000

        # ---- Z-score normalisation: edge features (per-feature) ----
        edge_feats = g.edata["e"].float()               # shape [E, 2]

        e1 = edge_feats[:, 0]
        e2 = edge_feats[:, 1]

        e1_norm = (e1 - e1.mean()) / e1.std().clamp(min=1e-9)
        e2_norm = (e2 - e2.mean()) / e2.std().clamp(min=1e-9)

        e1_norm = torch.round(e1_norm * 10000) / 10000
        e2_norm = torch.round(e2_norm * 10000) / 10000

        g.edata["e"] = torch.stack([e1_norm, e2_norm], dim=1)
        # g.edata["y"] remains as raw binary labels (no normalisation on labels)

        return g

    # ----------------------------------------------------------
    # 3e. Public entry point
    # ----------------------------------------------------------
    def build(self) -> dgl.DGLGraph:
        self._build_nodes()
        self._build_edges()
        g = self._build_dgl_graph()
        return g


# ============================================================
# 4.  SAVE UTILITIES
# ============================================================

def save_outputs(
    g:            dgl.DGLGraph,
    nodes:        dict,
    edges:        dict,
    kmer_counts:  Counter,
    kmer_to_idx:  dict,
):
    # Node summary
    with open("output/nodes.txt", "w") as f:
        f.write("node_idx\tkmer\tabundance\n")
        for idx, attrs in nodes.items():
            f.write(f"{idx}\t{attrs['kmer']}\t{attrs['abundance']}\n")

    # Edge summary
    with open("output/edges.txt", "w") as f:
        f.write("src_idx\tdst_idx\tedge_seq\tabundance_avg\toccurrence_sim\tlabel\n")
        for (src, dst), eattrs in edges.items():
            f.write(
                f"{src}\t{dst}\t{eattrs['edge_seq']}\t"
                f"{eattrs['abundance_avg']:.4f}\t{eattrs['occurrence_sim']:.4f}\t"
                f"{int(eattrs['label'])}\n"
            )

    # Label distribution
    labels = [int(e["label"]) for e in edges.values()]
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    with open("output/label_stats.txt", "w") as f:
        f.write(f"Total edges  : {len(labels)}\n")
        f.write(f"Positive (y=1): {n_pos}  ({100*n_pos/max(len(labels),1):.1f}%)\n")
        f.write(f"Negative (y=0): {n_neg}  ({100*n_neg/max(len(labels),1):.1f}%)\n")

    dgl.save_graphs("output/graph2.bin", [g])
    print("Saved: output/nodes.txt, output/edges.txt, output/label_stats.txt, output/graph.bin")


# ============================================================
# 5.  PIPELINE
# ============================================================

def build_all(reference_fasta: str, reads_fastq: str, k: int = 26):
    print("=" * 60)
    print(f"De Bruijn Graph  |  k={k}")
    print("=" * 60)

    # --- Parse ---
    print("\n[1] Parsing reference FASTA …")
    reference_seqs = parse_fasta(reference_fasta)
    print(f"    {len(reference_seqs)} reference sequence(s) loaded.")

    print("\n[2] Parsing reads FASTQ (quality scores ignored) …")
    read_seqs = parse_fastq_sequences_only(reads_fastq)
    print(f"    {len(read_seqs)} reads loaded.")

    # --- Count k-mers from READS only ---
    print(f"\n[3] Counting {k}-mers from reads …")
    kmer_counts = count_kmers(read_seqs, k)
    print(f"    Unique {k}-mers (nodes): {len(kmer_counts)}")

    # --- Build graph ---
    print("\n[4] Building De Bruijn graph …")
    builder = DeBruijnGraphBuilder(kmer_counts, reference_seqs, k)
    g = builder.build()

    # --- Save ---
    print("\n[5] Saving outputs …")
    save_outputs(g, builder.nodes, builder.edges, kmer_counts, builder.kmer_to_idx)

    return g, builder


# ============================================================
# 6.  MAIN
# ============================================================

if __name__ == "__main__":
    g, builder = build_all(
        reference_fasta="./curly-octo-train/reference.fasta",
        reads_fastq="./curly-octo-train/hifi_sim_small_0001.fq",
        k=26,
    )

    print("\n" + "=" * 60)
    print("GRAPH SUMMARY")
    print("=" * 60)
    print(f"  Nodes            : {g.num_nodes()}")
    print(f"  Edges            : {g.num_edges()}")
    print(f"  Node feature 'x' : {g.ndata['x'].shape}  (Z-score abundance)")
    print(f"  Edge feature 'e' : {g.edata['e'].shape}  (Z-score [avg_abund, occur_sim])")
    print(f"  Edge label   'y' : {g.edata['y'].shape}  (binary: 1=in ref, 0=not)")

    # Quick label stats from tensor
    y = g.edata["y"]
    print(f"\n  Positive edges (y=1) : {int(y.sum())} / {len(y)}")
    print(f"  Negative edges (y=0) : {int((y == 0).sum())} / {len(y)}")
    print("\nAll outputs saved to ./output/")
