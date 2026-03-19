#!/usr/bin/env python3
"""
Visualize the De Bruijn graph stored in output/graph.bin

Four visualizations produced:
  1. graph_structure.png      — NetworkX spring layout, nodes coloured by abundance
  2. node_abundance_hist.png  — histogram of z-score node features
  3. edge_feature_scatter.png — scatter of avg_abundance vs occurrence_similarity
  4. label_distribution.png   — bar chart of positive vs negative edges
  5. subgraph_top50.png       — zoomed view of 50 highest-abundance nodes + their edges

Install dependencies:
    pip install dgl torch networkx matplotlib seaborn
"""

import dgl
import torch
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors


# ── Load ──────────────────────────────────────────────────────────────────────
graphs, _ = dgl.load_graphs("output/graph2.bin")
g = graphs[0]

print("Graph loaded:")
print(f"  Nodes : {g.num_nodes()}")
print(f"  Edges : {g.num_edges()}")
print(f"  ndata : {dict(g.ndata)}")
print(f"  edata : {dict(g.edata)}")

node_z      = g.ndata["x"].squeeze().numpy()          # [N]      z-score abundance
edge_e      = g.edata["e"].numpy()                    # [E, 2]   z-score [avg_abund, occur_sim]
edge_labels = g.edata["y"].numpy().astype(int)        # [E]      0 or 1

src, dst = g.edges()
src = src.numpy()
dst = dst.numpy()

# ── Shared style ──────────────────────────────────────────────────────────────
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "axes.grid": True, "grid.color": "#e0e0e0", "font.size": 11})
FIGSIZE = (10, 7)

# =============================================================================
# 1. Full graph structure — spring layout, node colour = z-score abundance
# =============================================================================
print("\n[1] Drawing full graph structure …")

G_nx = dgl.to_networkx(g).to_undirected()   # convert to networkx for layout

# Spring layout can be slow on large graphs; cap at 1000 nodes for speed
MAX_NODES = 1000
if g.num_nodes() > MAX_NODES:
    print(f"    Graph has >{MAX_NODES} nodes — sampling {MAX_NODES} highest-abundance nodes.")
    top_idx = np.argsort(node_z)[-MAX_NODES:]
    G_nx    = G_nx.subgraph(top_idx)
    pass  # subgraph already set
else:
    pass  # G_nx already has all nodes

pos = nx.spring_layout(G_nx, seed=42, k=0.5)

# Always index node_z by original node ID to avoid re-index mismatch
nodes_in_sub = list(G_nx.nodes())
sub_z        = node_z[nodes_in_sub]
norm         = mcolors.Normalize(vmin=sub_z.min(), vmax=sub_z.max())
cmap         = cm.RdYlGn
node_colors  = [cmap(norm(node_z[n])) for n in nodes_in_sub]

fig, ax = plt.subplots(figsize=FIGSIZE)
nx.draw_networkx(
    G_nx, pos,
    node_color=node_colors,
    node_size=60,
    with_labels=False,
    edge_color="#aaaaaa",
    width=0.5,
    arrows=False,
    ax=ax,
)
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
plt.colorbar(sm, ax=ax, label="Node feature (z-score abundance)")
ax.set_title("De Bruijn Graph — spring layout\n(node colour = z-score k-mer abundance)")
ax.axis("off")
plt.tight_layout()
plt.savefig("output/graph_structure.png", dpi=150)
plt.close()
print("    Saved: output/graph_structure.png")


# =============================================================================
# 2. Node feature distribution — histogram of z-score abundances
# =============================================================================
print("[2] Node abundance histogram …")

fig, ax = plt.subplots(figsize=FIGSIZE)
ax.hist(node_z, bins=40, color="#4C72B0", edgecolor="white", linewidth=0.5)
ax.axvline(0, color="red", linestyle="--", linewidth=1.2, label="Mean (z=0)")
ax.set_xlabel("Node feature (z-score abundance)")
ax.set_ylabel("Number of nodes")
ax.set_title("Distribution of Node Abundances (Z-score)")
ax.legend()
plt.tight_layout()
plt.savefig("output/node_abundance_hist.png", dpi=150)
plt.close()
print("    Saved: output/node_abundance_hist.png")


# =============================================================================
# 3. Edge feature scatter — avg_abundance vs occurrence_similarity, coloured by label
# =============================================================================
print("[3] Edge feature scatter plot …")

e_avg   = edge_e[:, 0]   # z-score average abundance
e_sim   = edge_e[:, 1]   # z-score occurrence similarity
colors  = np.where(edge_labels == 1, "#2ecc71", "#e74c3c")   # green=positive, red=negative

fig, ax = plt.subplots(figsize=FIGSIZE)
scatter_neg = ax.scatter(e_avg[edge_labels == 0], e_sim[edge_labels == 0],
                          c="#e74c3c", alpha=0.5, s=15, label="y=0 (not in reference)")
scatter_pos = ax.scatter(e_avg[edge_labels == 1], e_sim[edge_labels == 1],
                          c="#2ecc71", alpha=0.7, s=15, label="y=1 (in reference)")
ax.set_xlabel("Edge feature 1 — z-score avg abundance")
ax.set_ylabel("Edge feature 2 — z-score occurrence similarity")
ax.set_title("Edge Features coloured by Label")
ax.legend()
plt.tight_layout()
plt.savefig("output/edge_feature_scatter.png", dpi=150)
plt.close()
print("    Saved: output/edge_feature_scatter.png")


# =============================================================================
# 4. Label distribution — bar chart
# =============================================================================
print("[4] Label distribution bar chart …")

n_pos = int((edge_labels == 1).sum())
n_neg = int((edge_labels == 0).sum())
total = len(edge_labels)

fig, ax = plt.subplots(figsize=(6, 5))
bars = ax.bar(
    ["y=0\n(not in reference)", "y=1\n(in reference)"],
    [n_neg, n_pos],
    color=["#e74c3c", "#2ecc71"],
    edgecolor="white",
    width=0.5,
)
for bar, count in zip(bars, [n_neg, n_pos]):
    pct = 100 * count / total
    ax.text(bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.5,
            f"{count}\n({pct:.1f}%)",
            ha="center", va="bottom", fontsize=11)
ax.set_ylabel("Number of edges")
ax.set_title("Edge Label Distribution")
ax.set_ylim(0, max(n_neg, n_pos) * 1.2)
plt.tight_layout()
plt.savefig("output/label_distribution.png", dpi=150)
plt.close()
print("    Saved: output/label_distribution.png")


# =============================================================================
# 5. Subgraph — top 50 highest-abundance nodes, directed edges coloured by label
# =============================================================================
print("[5] Subgraph of top-50 nodes …")

TOP_N = min(50, g.num_nodes())
top50_idx = set(np.argsort(node_z)[-TOP_N:].tolist())

# Keep only edges where BOTH endpoints are in top50
mask      = np.array([(s in top50_idx and d in top50_idx) for s, d in zip(src, dst)])
sub_src   = src[mask]
sub_dst   = dst[mask]
sub_label = edge_labels[mask]

G_sub = nx.DiGraph()
G_sub.add_nodes_from(top50_idx)
for s, d, lbl in zip(sub_src, sub_dst, sub_label):
    G_sub.add_edge(int(s), int(d), label=int(lbl))

pos_sub    = nx.spring_layout(G_sub, seed=7, k=1.2)
sub_node_z = np.array([node_z[n] for n in G_sub.nodes()])
norm_sub   = mcolors.Normalize(vmin=sub_node_z.min(), vmax=sub_node_z.max())
nc         = [cmap(norm_sub(node_z[n])) for n in G_sub.nodes()]

edge_colors_sub = ["#2ecc71" if d["label"] == 1 else "#e74c3c"
                   for _, _, d in G_sub.edges(data=True)]

fig, ax = plt.subplots(figsize=(12, 9))
nx.draw_networkx_nodes(G_sub, pos_sub, node_color=nc, node_size=300, ax=ax)
nx.draw_networkx_labels(G_sub, pos_sub, font_size=7, ax=ax)
nx.draw_networkx_edges(G_sub, pos_sub, edge_color=edge_colors_sub,
                       arrows=True, arrowsize=15, width=1.5, ax=ax,
                       connectionstyle="arc3,rad=0.1")

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#2ecc71", label="y=1 (in reference)"),
    Patch(facecolor="#e74c3c", label="y=0 (not in reference)"),
]
ax.legend(handles=legend_elements, loc="upper left")
sm2 = cm.ScalarMappable(cmap=cmap, norm=norm_sub)
sm2.set_array([])
plt.colorbar(sm2, ax=ax, label="Node z-score abundance")
ax.set_title(f"Subgraph — Top {TOP_N} highest-abundance nodes\n(directed, edge colour = label)")
ax.axis("off")
plt.tight_layout()
plt.savefig("output/subgraph_top50.png", dpi=150)
plt.close()
print("    Saved: output/subgraph_top50.png")

# =============================================================================
# Done
# =============================================================================
print("\n✓ All visualizations saved to ./output/")
print("""
output/
├── graph_structure.png       — full graph, spring layout, colour = z-score abundance
├── node_abundance_hist.png   — histogram of node z-score features
├── edge_feature_scatter.png  — edge features scatter, colour = label
├── label_distribution.png    — bar chart: positive vs negative edges
└── subgraph_top50.png        — top-50 nodes subgraph with directed edges
""")
