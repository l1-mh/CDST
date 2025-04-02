import pandas as pd
import numpy as np
from collections import defaultdict
import argparse
import os
import networkx as nx


def validate_mst(mst_df):
    """
    Validate if the input MST file is a valid MST.
    The MST must be a connected graph without any cycles.
    """
    G = nx.Graph()

    for _, row in mst_df.iterrows():
        G.add_edge(row['Node1'], row['Node2'], weight=row['Distance'])

    if not nx.is_connected(G):
        raise ValueError("The input MST is not a connected graph.")
    
    if G.number_of_edges() != G.number_of_nodes() - 1:
        raise ValueError("The input MST contains cycles or disconnected components.")
    
    print("The input MST is valid.")


def merge_by_category(mst, classification, epsilon):
    category_edges = defaultdict(list)
    category_counts = defaultdict(int)
    category_nodes = defaultdict(set)

    for _, row in mst.iterrows():
        source, target, distance = row['Node1'], row['Node2'], row['Distance']
        cat1 = classification.get(source, 'Unknown')
        cat2 = classification.get(target, 'Unknown')

        category_nodes[cat1].add(source)
        category_nodes[cat2].add(target)

        if cat1 == 'Unknown' or cat2 == 'Unknown':
            continue

        if cat1 != cat2:
            sorted_pair = tuple(sorted([cat1, cat2]))
            category_edges[sorted_pair].append(distance)
            category_counts[sorted_pair] += 1

    # Create edge list
    edge_result = []
    for (cat1, cat2), distances in category_edges.items():
        avg_distance = np.mean(distances)
        num_connections = category_counts[(cat1, cat2)]
        weight = 1 / (avg_distance + epsilon)
        edge_result.append([cat1, cat2, avg_distance, num_connections, weight])

    edge_df = pd.DataFrame(edge_result, columns=['Source', 'Target', 'Distance', 'Num_Connections', 'Weight'])

    # Create category nodes list
    node_result = []
    for category, nodes in category_nodes.items():
        node_result.append([category, len(nodes)])

    node_df = pd.DataFrame(node_result, columns=['ID', 'Num_Node'])

    return edge_df, node_df


def read_classification_file(classification_file):
    classification = {}
    with open(classification_file, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 2:
                classification[parts[0]] = parts[1]
    return classification


def main(args):
    # Read MST file
    try:
        mst = pd.read_csv(args.input)
    except Exception as e:
        raise ValueError(f"Failed to read the input file: {e}")
    
    # If no columns provided, assume the first three columns are Node1, Node2, Distance
    if not set(['Node1', 'Node2', 'Distance']).issubset(mst.columns):
        if len(mst.columns) >= 3:
            mst.columns = ['Node1', 'Node2', 'Distance']
        else:
            raise ValueError("Input file must have at least three columns for Node1, Node2, and Distance.")

    # Validate the input MST
    validate_mst(mst)

    # Read classification file if provided
    if args.classification:
        classification = read_classification_file(args.classification)
        edge_df, node_df = merge_by_category(mst, classification, args.epsilon)
        
        os.makedirs(args.output, exist_ok=True)
        
        edge_output = os.path.join(args.output, 'merged_edges.csv')
        node_output = os.path.join(args.output, 'category_nodes.csv')

        edge_df.to_csv(edge_output, index=False)
        node_df.to_csv(node_output, index=False)

        print(f"Merged edges saved to: {edge_output}")
        print(f"Category nodes saved to: {node_output}")
    else:
        # If no classification file, just save the original MST with weights
        mst['Weight'] = 1 / (mst['Distance'] + args.epsilon)
        mst_output = os.path.join(args.output, 'edges_with_weights.csv')
        mst.to_csv(mst_output, index=False)
        print(f"Edges with weights saved to: {mst_output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate and merge MST by category")
    parser.add_argument('-i', '--input', required=True, help="Input MST CSV file")
    parser.add_argument('-c', '--classification', help="Optional classification file")
    parser.add_argument('-o', '--output', required=True, help="Output folder for results")
    parser.add_argument('-e', '--epsilon', type=float, default=0.001, help="Small value to avoid division by zero")

    args = parser.parse_args()
    main(args)
