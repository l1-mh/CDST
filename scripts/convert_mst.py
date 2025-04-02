import pandas as pd
import numpy as np
import networkx as nx
from collections import defaultdict
import argparse
import os


def load_mst(file):
    mst = pd.read_csv(file)
    if len(mst.columns) < 3:
        raise ValueError("The input MST file must have at least 3 columns: Source, Target, Weight")
    mst.columns = ['Source', 'Target', 'Weight']
    return mst


def load_classification(file):
    classification = pd.read_csv(file)
    if len(classification.columns) < 2:
        raise ValueError("The classification file must have at least 2 columns: Sample, Category")
    classification.columns = ['Sample', 'Category']
    return classification.set_index('Sample').to_dict()['Category']


def validate_mst(mst):
    G = nx.Graph()
    G.add_weighted_edges_from(mst.values)
    if not nx.is_tree(G):
        raise ValueError("The input MST is not valid. The input graph is not a tree.")


def merge_by_category(mst, classification, epsilon):
    category_edges = defaultdict(list)
    category_counts = defaultdict(lambda: defaultdict(int))

    for _, row in mst.iterrows():
        source, target, weight = row['Source'], row['Target'], row['Weight']
        cat1 = classification.get(source, 'Unknown')
        cat2 = classification.get(target, 'Unknown')

        if cat1 == 'Unknown' or cat2 == 'Unknown':
            continue  # Skip if the category is not found

        if cat1 != cat2:
            sorted_pair = tuple(sorted([cat1, cat2]))
            category_edges[sorted_pair].append(weight)
            category_counts[cat1][cat2] += 1
            category_counts[cat2][cat1] += 1

    result = []

    for (cat1, cat2), weights in category_edges.items():
        avg_distance = np.mean(weights)
        num_connections = len(weights)
        weight = 1 / (avg_distance + epsilon)
        result.append([cat1, cat2, avg_distance, num_connections, weight])

    result_df = pd.DataFrame(result, columns=['Source_Category', 'Target_Category', 'Average_Distance', 'Num_Connections', 'Weight'])
    return result_df


def process_mst(mst, classification, epsilon):
    if classification:
        result_df = merge_by_category(mst, classification, epsilon)

        # Check if the merged result is a valid MST
        G = nx.Graph()
        for _, row in result_df.iterrows():
            G.add_edge(row['Source_Category'], row['Target_Category'], weight=row['Average_Distance'])

        if not nx.is_tree(G):
            mst_tree = nx.minimum_spanning_tree(G, weight='weight')
            edges = list(mst_tree.edges(data=True))
            result = []

            for edge in edges:
                source, target, data = edge
                avg_distance = data['weight']
                weight = 1 / (avg_distance + epsilon)
                result.append([source, target, avg_distance, 1, weight])

            result_df = pd.DataFrame(result, columns=['Source_Category', 'Target_Category', 'Average_Distance', 'Num_Connections', 'Weight'])
            print("Merged result was not a valid MST. Converted to MST using NetworkX.")

    else:
        # Calculate weight and include in output file
        mst['Weight'] = 1 / (mst['Weight'] + epsilon)
        result_df = mst[['Source', 'Target', 'Weight']].copy()
        result_df['Distance'] = mst['Weight']
        result_df = result_df[['Source', 'Target', 'Distance', 'Weight']]

    return result_df


def main(args):
    mst = load_mst(args.input)
    validate_mst(mst)

    classification = load_classification(args.classification) if args.classification else None

    result_df = process_mst(mst, classification, args.epsilon)
    result_df.to_csv(args.output, index=False)

    if args.classification:
        print(f"Category-merged output saved to {args.output}")
    else:
        print(f"Original MST output saved to {args.output} with calculated weights.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and merge MST by category, if provided.")
    parser.add_argument('-i', '--input', required=True, help="Input MST file (CSV format)")
    parser.add_argument('-o', '--output', required=True, help="Output file (CSV format)")
    parser.add_argument('-c', '--classification', help="Classification file (CSV format, with columns: Sample, Category)")
    parser.add_argument('-e', '--epsilon', type=float, default=0.001, help="Small epsilon value to avoid zero division (default: 0.001)")

    args = parser.parse_args()
    main(args)
