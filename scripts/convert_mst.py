import pandas as pd
import numpy as np
from collections import defaultdict
import argparse
import os


def merge_by_category(mst, classification, epsilon):
    category_edges = defaultdict(list)
    category_counts = defaultdict(int)
    category_nodes = defaultdict(set)

    for _, row in mst.iterrows():
        source, target, distance = row['Source'], row['Target'], row['Distance']
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

    edge_df = pd.DataFrame(edge_result, columns=['Source_Category', 'Target_Category', 'Average_Distance', 'Num_Connections', 'Weight'])

    # Create category nodes list
    node_result = []
    for category, nodes in category_nodes.items():
        node_result.append([category, len(nodes)])

    node_df = pd.DataFrame(node_result, columns=['Category', 'Num_Nodes'])

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
    mst = pd.read_csv(args.input)
    classification = read_classification_file(args.classification)

    epsilon = args.epsilon

    edge_df, node_df = merge_by_category(mst, classification, epsilon)

    os.makedirs(args.output, exist_ok=True)
    
    edge_output = os.path.join(args.output, 'merged_edges.csv')
    node_output = os.path.join(args.output, 'category_nodes.csv')

    edge_df.to_csv(edge_output, index=False)
    node_df.to_csv(node_output, index=False)

    print(f"Merged edges saved to: {edge_output}")
    print(f"Category nodes saved to: {node_output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge MST edges by category and output node counts")
    parser.add_argument('-i', '--input', required=True, help="Input MST CSV file")
    parser.add_argument('-c', '--classification', required=True, help="Classification file")
    parser.add_argument('-o', '--output', required=True, help="Output folder for results")
    parser.add_argument('-e', '--epsilon', type=float, default=0.001, help="Small value to avoid division by zero")

    args = parser.parse_args()
    main(args)
