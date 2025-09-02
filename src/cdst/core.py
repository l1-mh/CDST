"""
Core functions for CDST (CoDing Sequence Typer).
Implements MD5-based hashing of CDS sequences, distance matrix generation,
minimum spanning tree (MST) construction, hierarchical clustering (HC) trees,
and comparison utilities.
"""

import hashlib
import json
import os
from typing import Dict, List, Tuple

import pandas as pd
import networkx as nx
from Bio import SeqIO
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform

# src/cdst/core.py  (append at the end)

def run_full_pipeline(fasta_files, output_dir, min_cds_len=201, tree_mode="both", verbose=False):
    """
    Run the full CDST pipeline:
    1. Generate MD5 hashes from FASTA files
    2. Create comparison and difference matrices
    3. Optionally generate MST and/or HC trees

    Parameters
    ----------
    fasta_files : list of str
        List of input CDS FASTA files (.ffn)
    output_dir : str
        Directory to write results
    min_cds_len : int
        Minimum CDS length (default=201)
    tree_mode : str
        One of {"mst", "hc", "both"}
    verbose : bool
        Verbose logging
    """
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: generate md5
    md5_dict = {}
    for fasta_file in fasta_files:
        md5_dict[fasta_file] = generate_md5_for_fasta(
            fasta_file, min_cds_len=min_cds_len, verbose=verbose
        )
    json_output_path = os.path.join(output_dir, "md5_hashes.json")
    with open(json_output_path, "w") as f:
        json.dump(md5_dict, f, indent=4)
    print(f"[cdst] MD5 hashes written to {json_output_path}")

    # Step 2: comparison + difference matrices
    comparison_matrix = generate_comparison_matrix(md5_dict, verbose=verbose)
    comp_path = os.path.join(output_dir, "comparison_matrix.csv")
    comparison_matrix.to_csv(comp_path)
    print(f"[cdst] Comparison matrix written to {comp_path}")

    diff_matrix = calculate_difference_matrix(comparison_matrix)
    diff_path = os.path.join(output_dir, "difference_matrix.csv")
    diff_matrix.to_csv(diff_path)
    print(f"[cdst] Difference matrix written to {diff_path}")

    # Step 3: MST
    if tree_mode in ("mst", "both"):
        edge_list = generate_edge_list(diff_matrix)
        mst_edges = generate_mst(edge_list)
        mst_csv = os.path.join(output_dir, "mst.csv")
        with open(mst_csv, "w") as f:
            f.write("Node1,Node2,Distance\n")
            for u, v, data in mst_edges:
                f.write(f"{u},{v},{data['weight']}\n")
        print(f"[cdst] MST edges written to {mst_csv}")
        newick_str = mst_to_newick(mst_edges, list(diff_matrix.index))
        mst_newick = os.path.join(output_dir, "mst.newick")
        with open(mst_newick, "w") as f:
            f.write(newick_str)
        print(f"[cdst] MST Newick tree written to {mst_newick}")

    # Step 4: HC
    if tree_mode in ("hc", "both"):
        hc_tree = generate_hc_tree(diff_matrix)
        leaf_names = list(diff_matrix.index)
        newick_str = tree_to_newick(hc_tree, "", hc_tree.dist, leaf_names)
        hc_newick = os.path.join(output_dir, "hc.newick")
        with open(hc_newick, "w") as f:
            f.write(newick_str)
        print(f"[cdst] HC Newick tree written to {hc_newick}")


def generate_md5_for_fasta(fasta_file: str, min_cds_len: int = 201, verbose: bool = False) -> List[str]:
    """
    Read CDS FASTA (e.g., Prodigal .ffn) and return a list of MD5 hashes
    after filtering for ambiguous bases and minimum length.
    """
    md5_list = []
    kept, skipped_len, skipped_ambig = 0, 0, 0
    if verbose:
        print(f"[cdst] Processing file: {fasta_file} (min_cds_len={min_cds_len})")
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        if any(ch not in "ATCGatcg" for ch in seq):
            skipped_ambig += 1
            continue
        if len(seq) < max(0, int(min_cds_len)):
            skipped_len += 1
            continue
        seq = seq.upper()
        md5_hash = hashlib.md5(seq.encode()).hexdigest()
        md5_list.append(md5_hash)
        kept += 1
    if verbose:
        print(f"[cdst] Kept: {kept}, Skipped (len): {skipped_len}, Skipped (ambiguous): {skipped_ambig}")
    return md5_list


def generate_comparison_matrix(md5_dict: Dict[str, List[str]], verbose: bool = False) -> pd.DataFrame:
    """
    Generate a pairwise comparison matrix of shared CDS counts.
    """
    files = list(md5_dict.keys())
    matrix = []
    total_comparisons = len(files) * len(files)
    comparison_count = 0
    for file1 in files:
        row = []
        set1 = set(md5_dict[file1])
        for file2 in files:
            common_md5s = set1 & set(md5_dict[file2])
            row.append(len(common_md5s))
            if verbose:
                comparison_count += 1
                if comparison_count % 100 == 0:
                    print(f"[cdst] Comparing {comparison_count}/{total_comparisons} ...", end="\r")
        matrix.append(row)
    if verbose:
        print("\n[cdst] Comparison completed.")
    return pd.DataFrame(matrix, index=files, columns=files)


def calculate_difference_matrix(comparison_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Compute relative differences and symmetrize the distance matrix.
    """
    diff_matrix = comparison_matrix.copy().astype(float)
    for idx, row in comparison_matrix.iterrows():
        self_comp = row[idx]
        diff_matrix.loc[idx] = (self_comp - row) / self_comp
    for i in range(len(diff_matrix)):
        for j in range(i + 1, len(diff_matrix)):
            m = min(diff_matrix.iloc[i, j], diff_matrix.iloc[j, i])
            diff_matrix.iloc[i, j] = diff_matrix.iloc[j, i] = m
    return diff_matrix


def generate_edge_list(diff_matrix: pd.DataFrame) -> List[Tuple[str, str, float]]:
    """
    Convert a distance matrix into an edge list.
    """
    edge_list = []
    samples = diff_matrix.index
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            s1 = samples[i]
            s2 = samples[j]
            d = diff_matrix.loc[s1, s2]
            edge_list.append((s1, s2, d))
    return edge_list


def generate_mst(edge_list: List[Tuple[str, str, float]]) -> List[Tuple[str, str, dict]]:
    """
    Build a minimum spanning tree (MST) from an edge list.
    """
    G = nx.Graph()
    G.add_weighted_edges_from(edge_list)
    mst = nx.minimum_spanning_tree(G)
    return list(mst.edges(data=True))


def mst_to_newick(mst_edges: List[Tuple[str, str, dict]], leaf_names: List[str]) -> str:
    """
    Convert MST edges into a Newick-formatted string.
    """
    connections = {name: [] for name in leaf_names}
    for u, v, data in mst_edges:
        connections[u].append((v, data['weight']))
        connections[v].append((u, data['weight']))

    def build_newick(node, parent=None):
        children = [n for n, _ in connections[node] if n != parent]
        if not children:
            return node
        subtrees = []
        for child in children:
            w = [w for n, w in connections[node] if n == child][0]
            subtrees.append(build_newick(child, node) + ":%f" % w)
        return "(" + ",".join(subtrees) + ")" + node

    root = leaf_names[0]
    return build_newick(root) + ";"


def generate_hc_tree(diff_matrix: pd.DataFrame):
    """
    Build a hierarchical clustering tree using average linkage.
    """
    sym = diff_matrix.copy()
    for i in range(len(sym)):
        for j in range(i + 1, len(sym)):
            m = min(sym.iloc[i, j], sym.iloc[j, i])
            sym.iloc[i, j] = sym.iloc[j, i] = m
    condensed = squareform(sym)
    Z = linkage(condensed, method='average')
    tree, _ = to_tree(Z, rd=True)
    return tree


def tree_to_newick(node, newick: str, parentdist: float, leaf_names: List[str]) -> str:
    """
    Recursively convert a hierarchical clustering tree into Newick format.
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = tree_to_newick(node.get_left(), newick, node.dist, leaf_names)
        newick = tree_to_newick(node.get_right(), ",%s" % newick, node.dist, leaf_names)
        newick = "(%s" % newick
        return newick
