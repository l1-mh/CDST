
import argparse
import json
import os
import pandas as pd

from . import core


def main():
    parser = argparse.ArgumentParser(
        prog="cdst",
        description="CoDing Sequence Typer (CDST): MD5-based bacterial typing and clustering."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    
    # Subcommand: run (full pipeline)
    run = subparsers.add_parser("run", help="Run full CDST pipeline.")
    run.add_argument("-i", "--input", nargs="+", required=True, help="Input CDS FASTA files")
    run.add_argument("-o", "--output", required=True, help="Output directory")
    run.add_argument("-L", "--min-cds-len", type=int, default=201, help="Minimum CDS length")
    run.add_argument("-T", "--tree", choices=["mst", "hc", "both"], default="both",
                     help="Tree type to generate (default: both)")
    run.add_argument("-v", "--verbose", action="store_true", help="Verbose output")


    # Subcommand: generate
    gen = subparsers.add_parser("generate", help="Generate MD5 hash lists from FASTA files.")
    gen.add_argument("-i", "--input", nargs="+", required=True, help="Input CDS FASTA files")
    gen.add_argument("-o", "--output", required=True, help="Output directory")
    gen.add_argument("-L", "--min-cds-len", type=int, default=201, help="Minimum CDS length")
    gen.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    # Subcommand: matrix
    mat = subparsers.add_parser("matrix", help="Generate comparison and difference matrices.")
    mat.add_argument("-j", "--json", required=True, help="Input JSON file of MD5 hash lists")
    mat.add_argument("-o", "--output", required=True, help="Output directory")
    mat.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    # Subcommand: mst
    mst = subparsers.add_parser("mst", help="Generate MST and Newick tree from a difference matrix.")
    mst.add_argument("-m", "--matrix", required=True, help="Input difference matrix CSV")
    mst.add_argument("-o", "--output", required=True, help="Output directory")

    # Subcommand: hc
    hc = subparsers.add_parser("hc", help="Generate hierarchical clustering tree from a difference matrix.")
    hc.add_argument("-m", "--matrix", required=True, help="Input difference matrix CSV")
    hc.add_argument("-o", "--output", required=True, help="Output directory")

    args = parser.parse_args()

    if args.command == "generate":
        os.makedirs(args.output, exist_ok=True)
        md5_dict = {}
        for fasta_file in args.input:
            md5_dict[fasta_file] = core.generate_md5_for_fasta(
                fasta_file, min_cds_len=args.min_cds_len, verbose=args.verbose
            )
        json_output_path = os.path.join(args.output, "md5_hashes.json")
        with open(json_output_path, "w") as f:
            json.dump(md5_dict, f, indent=4)
        print(f"[cdst] MD5 hashes written to {json_output_path}")

    elif args.command == "matrix":
        with open(args.json, "r") as f:
            md5_dict = json.load(f)
        comparison_matrix = core.generate_comparison_matrix(md5_dict, verbose=args.verbose)
        os.makedirs(args.output, exist_ok=True)
        comp_path = os.path.join(args.output, "comparison_matrix.csv")
        comparison_matrix.to_csv(comp_path)
        print(f"[cdst] Comparison matrix written to {comp_path}")
        diff_matrix = core.calculate_difference_matrix(comparison_matrix)
        diff_path = os.path.join(args.output, "difference_matrix.csv")
        diff_matrix.to_csv(diff_path)
        print(f"[cdst] Difference matrix written to {diff_path}")

    elif args.command == "mst":
        diff_matrix = pd.read_csv(args.matrix, index_col=0)
        edge_list = core.generate_edge_list(diff_matrix)
        mst_edges = core.generate_mst(edge_list)
        os.makedirs(args.output, exist_ok=True)
        mst_csv = os.path.join(args.output, "mst.csv")
        with open(mst_csv, "w") as f:
            f.write("Node1,Node2,Distance\n")
            for u, v, data in mst_edges:
                f.write(f"{u},{v},{data['weight']}\n")
        print(f"[cdst] MST edges written to {mst_csv}")
        newick_str = core.mst_to_newick(mst_edges, list(diff_matrix.index))
        mst_newick = os.path.join(args.output, "mst.newick")
        with open(mst_newick, "w") as f:
            f.write(newick_str)
        print(f"[cdst] MST Newick tree written to {mst_newick}")

    elif args.command == "hc":
        diff_matrix = pd.read_csv(args.matrix, index_col=0)
        hc_tree = core.generate_hc_tree(diff_matrix)
        leaf_names = list(diff_matrix.index)
        newick_str = core.tree_to_newick(hc_tree, "", hc_tree.dist, leaf_names)
        os.makedirs(args.output, exist_ok=True)
        hc_newick = os.path.join(args.output, "hc.newick")
        with open(hc_newick, "w") as f:
            f.write(newick_str)
        print(f"[cdst] HC Newick tree written to {hc_newick}")
        
    elif args.command == "run":
        core.run_full_pipeline(
            args.input,
            args.output,
            min_cds_len=args.min_cds_len,
            tree_mode=args.tree,
            verbose=args.verbose,
        )

