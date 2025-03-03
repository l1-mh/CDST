#!/usr/bin/env python

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

def generate_md5_for_fasta(fasta_file, verbose=False):
    md5_list = []
    if verbose:
        print(f"Processing file: {fasta_file}")
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        if any(char not in "ATCGatcg" for char in sequence):
            continue 
        sequence = sequence.upper()
        md5_hash = hashlib.md5(sequence.encode()).hexdigest()
        md5_list.append(md5_hash)
    return md5_list

def generate_comparison_matrix(md5_dict, verbose=False):
    files = list(md5_dict.keys())
    matrix = []
    total_comparisons = len(files) * len(files)
    comparison_count = 0
    for file1 in files:
        row = []
        for file2 in files:
            common_md5s = set(md5_dict[file1]) & set(md5_dict[file2])
            row.append(len(common_md5s))
            if verbose:
                comparison_count += 1
                print(f"Computing {comparison_count}/{total_comparisons}...", end='\r')
        matrix.append(row)
    if verbose:
        print("Comparison completed.")
    return pd.DataFrame(matrix, index=files, columns=files)



def calculate_difference_matrix(comparison_matrix):
    diff_matrix = comparison_matrix.copy().astype(float)  
    for index, row in comparison_matrix.iterrows():
        self_comparison_value = row[index]
        diff_matrix.loc[index] = (self_comparison_value - row) / self_comparison_value
    return diff_matrix


def generate_edge_list(diff_matrix):
    edge_list = []
    samples = diff_matrix.index
    for i, sample1 in enumerate(samples):
        for j, sample2 in enumerate(samples):
            if i < j:
                distance = min(diff_matrix.loc[sample1, sample2], diff_matrix.loc[sample2, sample1])
                edge_list.append((sample1, sample2, distance))
    return edge_list

def generate_mst(edge_list):
    G = nx.Graph()
    G.add_weighted_edges_from(edge_list)
    mst = nx.minimum_spanning_tree(G)
    return list(mst.edges(data=True))

#def mst_to_newick(mst_edges, leaf_names):
#    connections = {name: [] for name in leaf_names}
#    for u, v, data in mst_edges:
#        connections[u].append((v, data['weight']))
#        connections[v].append((u, data['weight']))
#    def build_newick(node, parent=None):
#        children = [n for n, _ in connections[node] if n != parent]
#        if not children:
#            return node
#        subtrees = [build_newick(child, node) + ":%f" % connections[node][i][1] for i, child in enumerate(children)]
#        return "(" + ",".join(subtrees) + ")" + node
#    root = leaf_names[0]
#    return build_newick(root) + ";"

def generate_hc_tree(diff_matrix):
    sym_diff_matrix = diff_matrix.copy()
    for i in range(len(sym_diff_matrix)):
        for j in range(i + 1, len(sym_diff_matrix)):
            sym_diff_matrix.iloc[i, j] = sym_diff_matrix.iloc[j, i] = min(sym_diff_matrix.iloc[i, j], sym_diff_matrix.iloc[j, i])
    condensed_dist = squareform(sym_diff_matrix)
    Z = linkage(condensed_dist, method='average')
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

def generate_md5(args):
    md5_dict = {}
    for fasta_file in args.input:
        md5_hashes = generate_md5_for_fasta(fasta_file, verbose=args.verbose)
        md5_dict[fasta_file] = md5_hashes
    json_output_path = os.path.join(args.output, "md5_hashes.json")
    with open(json_output_path, 'w') as f:
        json.dump(md5_dict, f, indent=4)
    print(f"MD5 hashes have been written to {json_output_path}")

def generate_matrices(args):
    with open(args.json, 'r') as f:
        md5_dict = json.load(f)
    comparison_matrix = generate_comparison_matrix(md5_dict, verbose=args.verbose)
    csv_output_path = os.path.join(args.output, "comparison_matrix.csv")
    comparison_matrix.to_csv(csv_output_path)
    print(f"Comparison matrix has been written to {csv_output_path}")
    diff_matrix = calculate_difference_matrix(comparison_matrix)
    diff_matrix_output_path = os.path.join(args.output, "difference_matrix.csv")
    diff_matrix.to_csv(diff_matrix_output_path)
    print(f"Difference matrix has been written to {diff_matrix_output_path}")

def generate_mst_files(args):
    diff_matrix = pd.read_csv(args.matrix, index_col=0)
    edge_list = generate_edge_list(diff_matrix)
    edge_list_output_path = os.path.join(args.output, "edge_list.csv")
    with open(edge_list_output_path, 'w') as f:
        f.write("Sample1,Sample2,Distance\n")
        for sample1, sample2, distance in edge_list:
            f.write(f"{sample1},{sample2},{distance}\n")
    print(f"Edge list has been written to {edge_list_output_path}")
    mst_edges = generate_mst(edge_list)
    mst_csv_output_path = os.path.join(args.output, "mst.csv")
    with open(mst_csv_output_path, 'w') as f:
        f.write("Node1,Node2,Distance\n")
        for u, v, data in mst_edges:
            f.write(f"{u},{v},{data['weight']}\n")
    print(f"MST edge list has been written to {mst_csv_output_path}")
#    leaf_names = list(diff_matrix.index)
#    newick_str = mst_to_newick(mst_edges, leaf_names)
#    mst_newick_output_path = os.path.join(args.output, "mst.newick")
#    with open(mst_newick_output_path, 'w') as f:
#        f.write(newick_str)
#    print(f"Minimum Spanning Tree in Newick format has been written to {mst_newick_output_path}")

def generate_hc_tree_file(args):
    diff_matrix = pd.read_csv(args.matrix, index_col=0)
    hc_tree = generate_hc_tree(diff_matrix)
    leaf_names = list(diff_matrix.index)
    newick_str = tree_to_newick(hc_tree, "", hc_tree.dist, leaf_names)
    hc_output_path = os.path.join(args.output, "hc.newick")
    with open(hc_output_path, 'w') as f:
        f.write(newick_str)
    print(f"Hierarchical Clustering Tree has been written to {hc_output_path}")


def run_full_pipeline(args):

    md5_dict = {}
    for fasta_file in args.input:
        md5_hashes = generate_md5_for_fasta(fasta_file, verbose=args.verbose)
        md5_dict[fasta_file] = md5_hashes
    json_output_path = os.path.join(args.output, "md5_hashes.json")
    with open(json_output_path, 'w') as f:
        json.dump(md5_dict, f, indent=4)
    print(f"MD5 hashes have been written to {json_output_path}")

    comparison_matrix = generate_comparison_matrix(md5_dict, verbose=args.verbose)
    csv_output_path = os.path.join(args.output, "comparison_matrix.csv")
    comparison_matrix.to_csv(csv_output_path)
    print(f"Comparison matrix has been written to {csv_output_path}")

    diff_matrix = calculate_difference_matrix(comparison_matrix)
    diff_matrix_output_path = os.path.join(args.output, "difference_matrix.csv")
    diff_matrix.to_csv(diff_matrix_output_path)
    print(f"Difference matrix has been written to {diff_matrix_output_path}")

    if args.tree in ['mst', 'both']:
        edge_list = generate_edge_list(diff_matrix)
        edge_list_output_path = os.path.join(args.output, "edge_list.csv")
        with open(edge_list_output_path, 'w') as f:
            f.write("Sample1,Sample2,Distance\n")
            for sample1, sample2, distance in edge_list:
                f.write(f"{sample1},{sample2},{distance}\n")
        print(f"Edge list has been written to {edge_list_output_path}")

        mst_edges = generate_mst(edge_list)
        mst_csv_output_path = os.path.join(args.output, "mst.csv")
        with open(mst_csv_output_path, 'w') as f:
            f.write("Node1,Node2,Distance\n")
            for u, v, data in mst_edges:
                f.write(f"{u},{v},{data['weight']}\n")
        print(f"MST edge list has been written to {mst_csv_output_path}")

#        leaf_names = list(diff_matrix.index)
#        newick_str = mst_to_newick(mst_edges, leaf_names)
#        mst_newick_output_path = os.path.join(args.output, "mst.newick")
#        with open(mst_newick_output_path, 'w') as f:
#            f.write(newick_str)
#        print(f"Minimum Spanning Tree in Newick format has been written to {mst_newick_output_path}")

    if args.tree in ['hc', 'both']:
        hc_tree = generate_hc_tree(diff_matrix)
        leaf_names = list(diff_matrix.index)
        newick_str = tree_to_newick(hc_tree, "", hc_tree.dist, leaf_names)
        hc_output_path = os.path.join(args.output, "hc.newick")
        with open(hc_output_path, 'w') as f:
            f.write(newick_str)
        print(f"Hierarchical Clustering Tree has been written to {hc_output_path}")


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
        for existing_sample in combined_comparison_matrix.index:
            if new_sample != existing_sample:
                common_md5s = set(new_md5_dict[new_sample]) & set(new_md5_dict.get(existing_sample, []))
                combined_comparison_matrix.loc[new_sample, existing_sample] = len(common_md5s)
                combined_comparison_matrix.loc[existing_sample, new_sample] = len(common_md5s)
                if verbose:
                    print(f"Comparing {new_sample} with {existing_sample}")

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
    print(f"Combined MD5 hashes have been written to {combined_json_output_path}")

    if generate_matrix_bool:
        combined_comparison_matrix = merge_matrices(existing_comparison_matrices, combined_md5_dict, verbose=verbose)
        comparison_matrix_output_path = os.path.join(output_dir, "combined_comparison_matrix.csv")
        combined_comparison_matrix.to_csv(comparison_matrix_output_path)
        print(f"Combined comparison matrix has been written to {comparison_matrix_output_path}")

        combined_diff_matrix = calculate_difference_matrix(combined_comparison_matrix)
        diff_matrix_output_path = os.path.join(output_dir, "combined_difference_matrix.csv")
        combined_diff_matrix.to_csv(diff_matrix_output_path)
        print(f"Combined difference matrix has been written to {diff_matrix_output_path}")

        if generate_mst_bool:
            edge_list = generate_edge_list(combined_diff_matrix)
            edge_list_output_path = os.path.join(output_dir, "combined_edge_list.csv")
            with open(edge_list_output_path, 'w') as f:
                f.write("Sample1,Sample2,Distance\n")
                for sample1, sample2, distance in edge_list:
                    f.write(f"{sample1},{sample2},{distance}\n")
            print(f"Combined edge list has been written to {edge_list_output_path}")

            mst_edges = generate_mst(edge_list)
            mst_csv_output_path = os.path.join(output_dir, "combined_mst.csv")
            with open(mst_csv_output_path, 'w') as f:
                f.write("Node1,Node2,Distance\n")
                for u, v, data in mst_edges:
                    f.write(f"{u},{v},{data['weight']}\n")
            print(f"Combined MST edge list has been written to {mst_csv_output_path}")

#            leaf_names = list(combined_diff_matrix.index)
#            newick_str = mst_to_newick(mst_edges, leaf_names)
#            mst_newick_output_path = os.path.join(output_dir, "combined_mst.newick")
#            with open(mst_newick_output_path, 'w') as f:
#                f.write(newick_str)
#            print(f"Combined Minimum Spanning Tree in Newick format has been written to {mst_newick_output_path}")



def compare_new_samples_with_existing(new_md5_dict, existing_md5_dict, verbose=False):
    comparison_results_normalized = []
    comparison_results_unnormalized = []
    
    all_distances_normalized = []
    all_distances_unnormalized = []
    
    for new_sample, new_md5s in new_md5_dict.items():
        closest_sample_normalized = None
        closest_sample_unnormalized = None
        min_distance_normalized = float('inf')
        min_distance_unnormalized = float('inf')

        for existing_sample, existing_md5s in existing_md5_dict.items():
            common_md5s = set(new_md5s) & set(existing_md5s)
            
            # Calculate normalized distance
            distance_a_b_normalized = (len(new_md5s) - len(common_md5s)) / len(new_md5s)
            distance_b_a_normalized = (len(existing_md5s) - len(common_md5s)) / len(existing_md5s)
            distance_normalized = min(distance_a_b_normalized, distance_b_a_normalized)
            
            # Calculate unnormalized distance
            distance_a_b_unnormalized = len(new_md5s) - len(common_md5s)
            distance_b_a_unnormalized = len(existing_md5s) - len(common_md5s)
            distance_unnormalized = min(distance_a_b_unnormalized, distance_b_a_unnormalized)

            # Store distances for full comparison table
            all_distances_normalized.append((new_sample, existing_sample, distance_normalized))
            all_distances_unnormalized.append((new_sample, existing_sample, distance_unnormalized))

            if distance_normalized < min_distance_normalized:
                min_distance_normalized = distance_normalized
                closest_sample_normalized = existing_sample

            if distance_unnormalized < min_distance_unnormalized:
                min_distance_unnormalized = distance_unnormalized
                closest_sample_unnormalized = existing_sample

            if verbose:
                print(f"Comparing {new_sample} with {existing_sample}: Normalized Distance = {distance_normalized}, Unnormalized Distance = {distance_unnormalized}")

        comparison_results_normalized.append((new_sample, closest_sample_normalized, min_distance_normalized))
        comparison_results_unnormalized.append((new_sample, closest_sample_unnormalized, min_distance_unnormalized))

    result_output_path_normalized = "result.csv"
    result_output_path_unnormalized = "result_unnormalized.csv"
    
    with open(result_output_path_normalized, 'w') as f:
        f.write("NewSample,ClosestSample,NormalizedDistance\n")
        for new_sample, closest_sample, distance in comparison_results_normalized:
            f.write(f"{new_sample},{closest_sample},{distance}\n")

    with open(result_output_path_unnormalized, 'w') as f:
        f.write("NewSample,ClosestSample,UnnormalizedDistance\n")
        for new_sample, closest_sample, distance in comparison_results_unnormalized:
            f.write(f"{new_sample},{closest_sample},{distance}\n")

    distances_output_path_normalized = "distances.csv"
    distances_output_path_unnormalized = "distances_unnormalized.csv"
    
    with open(distances_output_path_normalized, 'w') as f:
        f.write("Sample1,Sample2,NormalizedDistance\n")
        for sample1, sample2, distance in all_distances_normalized:
            f.write(f"{sample1},{sample2},{distance}\n")

    with open(distances_output_path_unnormalized, 'w') as f:
        f.write("Sample1,Sample2,UnnormalizedDistance\n")
        for sample1, sample2, distance in all_distances_unnormalized:
            f.write(f"{sample1},{sample2},{distance}\n")

    print(f"Comparison results have been written to {result_output_path_normalized}, {result_output_path_unnormalized}, {distances_output_path_normalized}, and {distances_output_path_unnormalized}")
    
    return comparison_results_normalized, comparison_results_unnormalized, all_distances_normalized, all_distances_unnormalized


def test_new_samples(args):
    new_md5_dict = {}
    for fasta_file in args.input:
        md5_hashes = generate_md5_for_fasta(fasta_file, verbose=args.verbose)
        new_md5_dict[fasta_file] = md5_hashes

    existing_md5_dict = {}
    if os.path.exists(args.json):
        with open(args.json, 'r') as f:
            existing_md5_dict = json.load(f)

    comparison_results_normalized, comparison_results_unnormalized, all_distances_normalized, all_distances_unnormalized = compare_new_samples_with_existing(new_md5_dict, existing_md5_dict, verbose=args.verbose)

    comparison_output_path_normalized = os.path.join(args.output, "comparison_results_normalized.csv")
    comparison_output_path_unnormalized = os.path.join(args.output, "comparison_results_unnormalized.csv")

    with open(comparison_output_path_normalized, 'w') as f:
        f.write("NewSample,ClosestSample,NormalizedDistance\n")
        for new_sample, closest_sample, distance in comparison_results_normalized:
            f.write(f"{new_sample},{closest_sample},{distance}\n")
    
    with open(comparison_output_path_unnormalized, 'w') as f:
        f.write("NewSample,ClosestSample,UnnormalizedDistance\n")
        for new_sample, closest_sample, distance in comparison_results_unnormalized:
            f.write(f"{new_sample},{closest_sample},{distance}\n")

    print(f"Comparison results have been written to {comparison_output_path_normalized} and {comparison_output_path_unnormalized}")
    
    distances_output_path_normalized = os.path.join(args.output, "distances_normalized.csv")
    distances_output_path_unnormalized = os.path.join(args.output, "distances_unnormalized.csv")

    with open(distances_output_path_normalized, 'w') as f:
        f.write("Sample1,Sample2,NormalizedDistance\n")
        for sample1, sample2, distance in all_distances_normalized:
            f.write(f"{sample1},{sample2},{distance}\n")
    
    with open(distances_output_path_unnormalized, 'w') as f:
        f.write("Sample1,Sample2,UnnormalizedDistance\n")
        for sample1, sample2, distance in all_distances_unnormalized:
            f.write(f"{sample1},{sample2},{distance}\n")

    print(f"Full distance comparison tables have been written to {distances_output_path_normalized} and {distances_output_path_unnormalized}")

    if args.mst and os.path.exists(args.mst):
        diff_matrix = pd.read_csv(args.mst.replace('mst.csv', 'difference_matrix.csv'), index_col=0)
        for new_sample, closest_sample, distance in comparison_results_normalized:
            diff_matrix.loc[new_sample] = float('inf')
            diff_matrix.loc[:, new_sample] = float('inf')
            diff_matrix.loc[new_sample, closest_sample] = distance
            diff_matrix.loc[closest_sample, new_sample] = distance

        edge_list = generate_edge_list(diff_matrix)
        mst_edges = generate_mst(edge_list)

        mst_csv_output_path = os.path.join(args.output, "updated_mst.csv")
        with open(mst_csv_output_path, 'w') as f:
            f.write("Node1,Node2,Distance\n")
            for u, v, data in mst_edges:
                f.write(f"{u},{v},{data['weight']}\n")
        print(f"Updated MST edge list has been written to {mst_csv_output_path}")



def parse_arguments():
    parser = argparse.ArgumentParser(description="Process FASTA files and generate various outputs.")
    subparsers = parser.add_subparsers(dest='command', required=True)

    generate_parser = subparsers.add_parser('generate', help='Generate MD5 hashes from FASTA files')
    generate_parser.add_argument('-i', '--input', type=str, nargs='+', required=True, help="Input FASTA files")
    generate_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder to save results")
    generate_parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output")

    matrix_parser = subparsers.add_parser('matrix', help='Generate comparison and difference matrices from JSON')
    matrix_parser.add_argument('-j', '--json', type=str, required=True, help="Input JSON file with MD5 hashes")
    matrix_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder to save results")
    matrix_parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output")

    mst_parser = subparsers.add_parser('mst', help='Generate MST and edge list from difference matrix')
    mst_parser.add_argument('-m', '--matrix', type=str, required=True, help="Input CSV file with difference matrix")
    mst_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder to save results")

    hc_parser = subparsers.add_parser('hc', help='Generate HC tree from difference matrix')
    hc_parser.add_argument('-m', '--matrix', type=str, required=True, help="Input CSV file with difference matrix")
    hc_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder to save results")

    run_parser = subparsers.add_parser('run', help='Run full pipeline from FASTA to tree generation')
    run_parser.add_argument('-i', '--input', type=str, nargs='+', required=True, help="Input FASTA files")
    run_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder to save results")
    run_parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output")
    run_parser.add_argument('-T', '--tree', choices=['mst', 'hc', 'both'], help="Generate tree: mst for Minimum Spanning Tree, hc for Hierarchical Clustering, both for both trees")

    join_parser = subparsers.add_parser('join', help='Join multiple JSON files and optionally generate matrices and MST')
    join_parser.add_argument('-d', '--inputdirs', type=str, nargs='+', required=True, help="Input directories containing JSON files")
    join_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder to save results")
    join_parser.add_argument('--matrix', action='store_true', help="Generate combined comparison and difference matrices")
    join_parser.add_argument('--mst', action='store_true', help="Generate combined MST")

    test_parser = subparsers.add_parser('test', help='Test new samples against existing JSON data')
    test_parser.add_argument('-i', '--input', type=str, nargs='+', required=True, help="Input FASTA files")
    test_parser.add_argument('-j', '--json', type=str, required=True, help="Existing JSON file with MD5 hashes")
    test_parser.add_argument('--mst', type=str, help="Optional existing MST CSV file")
    test_parser.add_argument('-o', '--output', type=str, required=True, help="Output folder to save results")
    test_parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose output")

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
