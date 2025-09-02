#!/usr/bin/env python

# ============================================
# Version: v0.2.0
# Changelog:
# - Add CDS length filtering (default: keep CDS length >= 201 nt)
# - CLI: new option --min-cds-len/-L for subcommands that read FASTA
# ============================================

__version__ = "0.2.0"

import hashlib
import argparse
import json
import os
import pandas as pd
from Bio import SeqIO
import networkx as nx
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
from io import StringIO

def generate_md5_for_fasta(fasta_file, min_cds_len=201, verbose=False):
    """
    Read CDS FASTA (e.g., Prodigal .ffn) and return MD5 list after filtering.
    - Filters out sequences containing non-ATCG
    - Filters out sequences with length < min_cds_len (nt)
    """
    md5_list = []
    kept, skipped_len, skipped_ambig = 0, 0, 0
    if verbose:
        print(f"[cdst] Processing file: {fasta_file} (min_cds_len={min_cds_len})")
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        # Ambiguous bases filter
        if any(ch not in "ATCGatcg" for ch in seq):
            skipped_ambig += 1
            continue
        # Length filter
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

def generate_comparison_matrix(md5_dict, verbose=False):
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

def calculate_difference_matrix(comparison_matrix):
    # Relative distance (directional), then symmetrize by min()
    diff_matrix = comparison_matrix.copy().astype(float)
    for idx, row in comparison_matrix.iterrows():
        self_comp = row[idx]
        diff_matrix.loc[idx] = (self_comp - row) / self_comp
    # Symmetrize by minimum (v0.1.1 behavior)
    for i in range(len(diff_matrix)):
        for j in range(i + 1, len(diff_matrix)):
            m = min(diff_matrix.iloc[i, j], diff_matrix.iloc[j, i])
            diff_matrix.iloc[i, j] = diff_matrix.iloc[j, i] = m
    return diff_matrix

def generate_edge_list(diff_matrix):
    edge_list = []
    samples = diff_matrix.index
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            s1 = samples[i]
            s2 = samples[j]
            d = diff_matrix.loc[s1, s2]
            edge_list.append((s1, s2, d))
    return edge_list

def generate_mst(edge_list):
    G = nx.Graph()
    G.add_weighted_edges_from(edge_list)
    mst = nx.minimum_spanning_tree(G)
    return list(mst.edges(data=True))

def mst_to_newick(mst_edges, leaf_names):
    connections = {name: [] for name in leaf_names}
    for u, v, data in mst_edges:
        connections[u].append((v, data['weight']))
        connections[v].append((u, data['weight']))
    def build_newick(node, parent=None):
        children = [n for n, _ in connections[node] if n != parent]
        if not children:
            return node
        subtrees = []
        for i, child in enumerate(children):
            w = [w for n, w in connections[node] if n == child][0]
            subtrees.append(build_newick(child, node) + ":%f" % w)
        return "(" + ",".join(subtrees) + ")" + node
    root = leaf_names[0]
    return build_newick(root) + ";"

def generate_hc_tree(diff_matrix):
    # Ensure symmetry (min)
    sym = diff_matrix.copy()
    for i in range(len(sym)):
        for j in range(i + 1, len(sym)):
            m = min(sym.iloc[i, j], sym.iloc[j, i])
            sym.iloc[i, j] = sym.iloc[j, i] = m
    condensed = squareform(sym)
    Z = linkage(condensed, method='average')
    tree, _ = to_tree(Z, rd=True)
    return tree

def tree_to_newick(node, newick, parentdist, leaf_names):
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

# --------- Subcommand implementations ---------

def generate_md5(args):
    md5_dict = {}
    for fasta_file in args.input:
        md5_hashes = generate_md5_for_fasta(fasta_file, min_cds_len=args.min_cds_len, verbose=args.verbose)
        md5_dict[fasta_file] = md5_hashes
    json_output_path = os.path.join(args.output, "md5_hashes.json")
    with open(json_output_path, 'w') as f:
        json.dump(md5_dict, f, indent=4)
    print(f"[cdst] MD5 hashes written: {json_output_path}")

def generate_matrices(args):
    with open(args.json, 'r') as f:
        md5_dict = json.load(f)
    comparison_matrix = generate_comparison_matrix(md5_dict, verbose=args.verbose)
    csv_output_path = os.path.join(args.output, "comparison_matrix.csv")
    comparison_matrix.to_csv(csv_output_path)
    print(f"[cdst] Comparison matrix: {csv_output_path}")
    diff_matrix = calculate_difference_matrix(comparison_matrix)
    diff_matrix_output_path = os.path.join(args.output, "difference_matrix.csv")
    diff_matrix.to_csv(diff_matrix_output_path)
    print(f"[cdst] Difference matrix: {diff_matrix_output_path}")

def generate_mst_files(args):
    diff_matrix = pd.read_csv(args.matrix, index_col=0)
    edge_list = generate_edge_list(diff_matrix)
    edge_list_output_path = os.path.join(args.output, "edge_list.csv")
    with open(edge_list_output_path, 'w') as f:
        f.write("Sample1,Sample2,Distance\n")
        for s1, s2, d in edge_list:
            f.write(f"{s1},{s2},{d}\n")
    print(f"[cdst] Edge list: {edge_list_output_path}")
    mst_edges = generate_mst(edge_list)
    mst_csv_output_path = os.path.join(args.output, "mst.csv")
    with open(mst_csv_output_path, 'w') as f:
        f.write("Node1,Node2,Distance\n")
        for u, v, data in mst_edges:
            f.write(f"{u},{v},{data['weight']}\n")
    print(f"[cdst] MST edges: {mst_csv_output_path}")
    leaf_names = list(diff_matrix.index)
    newick_str = mst_to_newick(mst_edges, leaf_names)
    mst_newick_output_path = os.path.join(args.output, "mst.newick")
    with open(mst_newick_output_path, 'w') as f:
        f.write(newick_str)
    print(f"[cdst] MST (Newick): {mst_newick_output_path}")

def generate_hc_tree_file(args):
    diff_matrix = pd.read_csv(args.matrix, index_col=0)
    hc_tree = generate_hc_tree(diff_matrix)
    leaf_names = list(diff_matrix.index)
    newick_str = tree_to_newick(hc_tree, "", hc_tree.dist, leaf_names)
    hc_output_path = os.path.join(args.output, "hc.newick")
    with open(hc_output_path, 'w') as f:
        f.write(newick_str)
    print(f"[cdst] HC tree (Newick): {hc_output_path}")

def run_full_pipeline(args):
    md5_dict = {}
    for fasta_file in args.input:
        md5_hashes = generate_md5_for_fasta(fasta_file, min_cds_len=args.min_cds_len, verbose=args.verbose)
        md5_dict[fasta_file] = md5_hashes
    json_output_path = os.path.join(args.output, "md5_hashes.json")
    with open(json_output_path, 'w') as f:
        json.dump(md5_dict, f, indent=4)
    print(f"[cdst] MD5 hashes written: {json_output_path}")

    comparison_matrix = generate_comparison_matrix(md5_dict, verbose=args.verbose)
    csv_output_path = os.path.join(args.output, "comparison_matrix.csv")
    comparison_matrix.to_csv(csv_output_path)
    print(f"[cdst] Comparison matrix: {csv_output_path}")

    diff_matrix = calculate_difference_matrix(comparison_matrix)
    diff_matrix_output_path = os.path.join(args.output, "difference_matrix.csv")
    diff_matrix.to_csv(diff_matrix_output_path)
    print(f"[cdst] Difference matrix: {diff_matrix_output_path}")

    if args.tree in ['mst', 'both']:
        edge_list = generate_edge_list(diff_matrix)
        edge_list_output_path = os.path.join(args.output, "edge_list.csv")
        with open(edge_list_output_path, 'w') as f:
            f.write("Sample1,Sample2,Distance\n")
            for s1, s2, d in edge_list:
                f.write(f"{s1},{s2},{d}\n")
        print(f"[cdst] Edge list: {edge_list_output_path}")

        mst_edges = generate_mst(edge_list)
        mst_csv_output_path = os.path.join(args.output, "mst.csv")
        with open(mst_csv_output_path, 'w') as f:
            f.write("Node1,Node2,Distance\n")
            for u, v, data in mst_edges:
                f.write(f"{u},{v},{data['weight']}\n")
        print(f"[cdst] MST edges: {mst_csv_output_path}")

        leaf_names = list(diff_matrix.index)
        newick_str = mst_to_newick(mst_edges, leaf_names)
        mst_newick_output_path = os.path.join(args.output, "mst.newick")
        with open(mst_newick_output_path, 'w') as f:
            f.write(newick_str)
        print(f"[cdst] MST (Newick): {mst_newick_output_path}")

    if args.tree in ['hc', 'both']:
        hc_tree = generate_hc_tree(diff_matrix)
        leaf_names = list(diff_matrix.index)
        newick_str = tree_to_newick(hc_tree, "", hc_tree.dist, leaf_names)
        hc_output_path = os.path.join(args.output, "hc.newick")
        with open(hc_output_path, 'w') as f:
            f.write(newick_str)
        print(f"[cdst] HC tree (Newick): {hc_output_path}")

def merge_matrices(existing_matrices, new_md5_dict, verbose=False):
    all_samples = set()
    for matrix in existing_matrices:
        all_samples.update(matrix.index)
    all_samples.update(new_md5_dict.keys())
    all_samples = sorted(all_samples)

    combined_comparison_matrix = pd.DataFrame(0, index=all_samples, columns=all_samples)

    for matrix in existing_matrices:
        for i in matrix.index:
            for j in matrix.columns:
                combined_comparison_matrix.loc[i, j] = matrix.loc[i, j]

    for new_sample in new_md5_dict.keys():
        set_new = set(new_md5_dict[new_sample])
        for existing_sample in combined_comparison_matrix.index:
            if new_sample == existing_sample:
                continue
            set_exist = set(new_md5_dict.get(existing_sample, []))
            common = set_new & set_exist
            combined_comparison_matrix.loc[new_sample, existing_sample] = len(common)
            combined_comparison_matrix.loc[existing_sample, new_sample] = len(common)
            if verbose:
                print(f"[cdst] Comparing {new_sample} vs {existing_sample}")

    return combined_comparison_matrix

def join_json_files(input_dirs, output_dir, generate_matrix_bool=False, generate_mst_bool=False, verbose=False):
    combined_md5_dict = {}
    existing_comparison_matrices = []

    for input_dir in input_dirs:
        json_file_path = os.path.join(input_dir, "md5_hashes.json")
        if os.path.exists(json_file_path):
            with open(json_file_path, 'r') as f:
                md5_dict = json.load(f)
                combined_md5_dict.update(md5_dict)

        if generate_matrix_bool:
            comparison_matrix_path = os.path.join(input_dir, "comparison_matrix.csv")
            if os.path.exists(comparison_matrix_path):
                existing_comparison_matrices.append(pd.read_csv(comparison_matrix_path, index_col=0))

    combined_json_output_path = os.path.join(output_dir, "combined_md5_hashes.json")
    with open(combined_json_output_path, 'w') as f:
        json.dump(combined_md5_dict, f, indent=4)
    print(f"[cdst] Combined MD5 hashes: {combined_json_output_path}")

    if generate_matrix_bool:
        combined_comparison_matrix = merge_matrices(existing_comparison_matrices, combined_md5_dict, verbose=verbose)
        comparison_matrix_output_path = os.path.join(output_dir, "combined_comparison_matrix.csv")
        combined_comparison_matrix.to_csv(comparison_matrix_output_path)
        print(f"[cdst] Combined comparison matrix: {comparison_matrix_output_path}")

        combined_diff_matrix = calculate_difference_matrix(combined_comparison_matrix)
        diff_matrix_output_path = os.path.join(output_dir, "combined_difference_matrix.csv")
        combined_diff_matrix.to_csv(diff_matrix_output_path)
        print(f"[cdst] Combined difference matrix: {diff_matrix_output_path}")

        if generate_mst_bool:
            edge_list = generate_edge_list(combined_diff_matrix)
            edge_list_output_path = os.path.join(output_dir, "combined_edge_list.csv")
            with open(edge_list_output_path, 'w') as f:
                f.write("Sample1,Sample2,Distance\n")
                for s1, s2, d in edge_list:
                    f.write(f"{s1},{s2},{d}\n")
            print(f"[cdst] Combined edge list: {edge_list_output_path}")

            mst_edges = generate_mst(edge_list)
            mst_csv_output_path = os.path.join(output_dir, "combined_mst.csv")
            with open(mst_csv_output_path, 'w') as f:
                f.write("Node1,Node2,Distance\n")
                for u, v, data in mst_edges:
                    f.write(f"{u},{v},{data['weight']}\n")
            print(f"[cdst] Combined MST edges: {mst_csv_output_path}")

            leaf_names = list(combined_diff_matrix.index)
            newick_str = mst_to_newick(mst_edges, leaf_names)
            mst_newick_output_path = os.path.join(output_dir, "combined_mst.newick")
            with open(mst_newick_output_path, 'w') as f:
                f.write(newick_str)
            print(f"[cdst] Combined MST (Newick): {mst_newick_output_path}")

def compare_new_samples_with_existing(new_md5_dict, existing_md5_dict, verbose=False):
    comparison_results_normalized = []
    comparison_results_unnormalized = []
    all_distances_normalized = []
    all_distances_unnormalized = []

    for new_sample, new_md5s in new_md5_dict.items():
        closest_norm = None
        closest_unnorm = None
        min_norm = float('inf')
        min_unnorm = float('inf')
        set_new = set(new_md5s)

        for existing_sample, existing_md5s in existing_md5_dict.items():
            set_exist = set(existing_md5s)
            common = set_new & set_exist

            # normalized
            d_a_b_norm = (len(set_new) - len(common)) / len(set_new) if set_new else 1.0
            d_b_a_norm = (len(set_exist) - len(common)) / len(set_exist) if set_exist else 1.0
            d_norm = min(d_a_b_norm, d_b_a_norm)

            # unnormalized
            d_a_b_unn = len(set_new) - len(common)
            d_b_a_unn = len(set_exist) - len(common)
            d_unn = min(d_a_b_unn, d_b_a_unn)

            all_distances_normalized.append((new_sample, existing_sample, d_norm))
            all_distances_unnormalized.append((new_sample, existing_sample, d_unn))

            if d_norm < min_norm:
                min_norm = d_norm
                closest_norm = existing_sample
            if d_unn < min_unnorm:
                min_unnorm = d_unn
                closest_unnorm = existing_sample

            if verbose:
                print(f"[cdst] {new_sample} vs {existing_sample}: norm={d_norm}, unnorm={d_unn}")

        comparison_results_normalized.append((new_sample, closest_norm, min_norm))
        comparison_results_unnormalized.append((new_sample, closest_unnorm, min_unnorm))

    # write simple CSVs in CWD (legacy behavior)
    with open("result.csv", 'w') as f:
        f.write("NewSample,ClosestSample,NormalizedDistance\n")
        for a, b, d in comparison_results_normalized:
            f.write(f"{a},{b},{d}\n")
    with open("result_unnormalized.csv", 'w') as f:
        f.write("NewSample,ClosestSample,UnnormalizedDistance\n")
        for a, b, d in comparison_results_unnormalized:
            f.write(f"{a},{b},{d}\n")
    with open("distances.csv", 'w') as f:
        f.write("Sample1,Sample2,NormalizedDistance\n")
        for a, b, d in all_distances_normalized:
            f.write(f"{a},{b},{d}\n")
    with open("distances_unnormalized.csv", 'w') as f:
        f.write("Sample1,Sample2,UnnormalizedDistance\n")
        for a, b, d in all_distances_unnormalized:
            f.write(f"{a},{b},{d}\n")

    print("[cdst] Wrote: result.csv, result_unnormalized.csv, distances.csv, distances_unnormalized.csv")
    return comparison_results_normalized, comparison_results_unnormalized, all_distances_normalized, all_distances_unnormalized

def test_new_samples(args):
    new_md5_dict = {}
    for fasta_file in args.input:
        md5_hashes = generate_md5_for_fasta(fasta_file, min_cds_len=args.min_cds_len, verbose=args.verbose)
        new_md5_dict[fasta_file] = md5_hashes

    existing_md5_dict = {}
    if os.path.exists(args.json):
        with open(args.json, 'r') as f:
            existing_md5_dict = json.load(f)

    comparison_results_normalized, comparison_results_unnormalized, all_distances_normalized, all_distances_unnormalized = compare_new_samples_with_existing(
        new_md5_dict, existing_md5_dict, verbose=args.verbose)

    # write into output dir
    os.makedirs(args.output, exist_ok=True)
    out_norm = os.path.join(args.output, "comparison_results_normalized.csv")
    out_unn = os.path.join(args.output, "comparison_results_unnormalized.csv")
    with open(out_norm, 'w') as f:
        f.write("NewSample,ClosestSample,NormalizedDistance\n")
        for a, b, d in comparison_results_normalized:
            f.write(f"{a},{b},{d}\n")
    with open(out_unn, 'w') as f:
        f.write("NewSample,ClosestSample,UnnormalizedDistance\n")
        for a, b, d in comparison_results_unnormalized:
            f.write(f"{a},{b},{d}\n")
    print(f"[cdst] Wrote: {out_norm} , {out_unn}")

    dist_out_norm = os.path.join(args.output, "distances_normalized.csv")
    dist_out_unn = os.path.join(args.output, "distances_unnormalized.csv")
    with open(dist_out_norm, 'w') as f:
        f.write("Sample1,Sample2,NormalizedDistance\n")
        for a, b, d in all_distances_normalized:
            f.write(f"{a},{b},{d}\n")
    with open(dist_out_unn, 'w') as f:
        f.write("Sample1,Sample2,UnnormalizedDistance\n")
        for a, b, d in all_distances_unnormalized:
            f.write(f"{a},{b},{d}\n")
    print(f"[cdst] Wrote: {dist_out_norm} , {dist_out_unn}")

    # optional MST update (kept as-is)

# --------- CLI ---------

def parse_arguments():
    parser = argparse.ArgumentParser(description="CDST: MD5 hash-based CDS comparison and clustering.")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    subparsers = parser.add_subparsers(dest='command', required=True)

    # generate
    generate_parser = subparsers.add_parser('generate', help='Generate MD5 hashes from CDS FASTA files')
    generate_parser.add_argument('-i', '--input', type=str, nargs='+', required=True, help="Input CDS FASTA files (e.g., .ffn)")
    generate_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder")
    generate_parser.add_argument('-L', '--min-cds-len', type=int, default=201, help="Minimum CDS length to keep (nt; default: 201)")
    generate_parser.add_argument('-v', '--verbose', action='store_true', help="Verbose output")

    # matrix
    matrix_parser = subparsers.add_parser('matrix', help='Generate comparison/difference matrices from JSON')
    matrix_parser.add_argument('-j', '--json', type=str, required=True, help="Input JSON with MD5 hashes")
    matrix_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder")
    matrix_parser.add_argument('-v', '--verbose', action='store_true', help="Verbose output")

    # mst
    mst_parser = subparsers.add_parser('mst', help='Generate MST and edge list from difference matrix')
    mst_parser.add_argument('-m', '--matrix', type=str, required=True, help="Difference matrix CSV")
    mst_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder")

    # hc
    hc_parser = subparsers.add_parser('hc', help='Generate HC tree from difference matrix')
    hc_parser.add_argument('-m', '--matrix', type=str, required=True, help="Difference matrix CSV")
    hc_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder")

    # run
    run_parser = subparsers.add_parser('run', help='Run full pipeline from CDS FASTA to trees')
    run_parser.add_argument('-i', '--input', type=str, nargs='+', required=True, help="Input CDS FASTA files (e.g., .ffn)")
    run_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder")
    run_parser.add_argument('-L', '--min-cds-len', type=int, default=201, help="Minimum CDS length to keep (nt; default: 201)")
    run_parser.add_argument('-v', '--verbose', action='store_true', help="Verbose output")
    run_parser.add_argument('-T', '--tree', choices=['mst', 'hc', 'both'], help="Tree to generate")

    # join
    join_parser = subparsers.add_parser('join', help='Join multiple JSONs; optionally build matrices and MST')
    join_parser.add_argument('-d', '--inputdirs', type=str, nargs='+', required=True, help="Input directories containing md5_hashes.json")
    join_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder")
    join_parser.add_argument('--matrix', action='store_true', help="Generate combined matrices")
    join_parser.add_argument('--mst', action='store_true', help="Generate combined MST")

    # test
    test_parser = subparsers.add_parser('test', help='Compare new samples against existing JSON')
    test_parser.add_argument('-i', '--input', type=str, nargs='+', required=True, help="Input CDS FASTA files (e.g., .ffn)")
    test_parser.add_argument('-j', '--json', type=str, required=True, help="Existing JSON with MD5 hashes")
    test_parser.add_argument('--mst', type=str, help="Optional existing MST CSV file")
    test_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder")
    test_parser.add_argument('-L', '--min-cds-len', type=int, default=201, help="Minimum CDS length to keep (nt; default: 201)")
    test_parser.add_argument('-v', '--verbose', action='store_true', help="Verbose output")

    return parser.parse_args()

def main():
    args = parse_arguments()
    os.makedirs(args.output, exist_ok=True)

    if args.command == 'generate':
        generate_md5(args)
    elif args.command == 'matrix':
        generate_matrices(args)
    elif args.command == 'mst':
        generate_mst_files(args)
    elif args.command == 'hc':
        generate_hc_tree_file(args)
    elif args.command == 'run':
        run_full_pipeline(args)
    elif args.command == 'join':
        join_json_files(args.inputdirs, args.output, generate_matrix_bool=args.matrix, generate_mst_bool=args.mst)
    elif args.command == 'test':
        test_new_samples(args)

if __name__ == "__main__":
    main()
