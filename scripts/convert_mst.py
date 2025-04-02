import pandas as pd
import networkx as nx
import argparse

def validate_mst(edges):
    """
    Validate if the given edge list forms a valid Minimum Spanning Tree (MST).
    """
    # Create an undirected graph
    G = nx.Graph()
    G.add_weighted_edges_from(edges)

    # Check if the graph is connected
    if not nx.is_connected(G):
        print("The graph is NOT connected. Therefore, it is NOT a valid MST.")
        return False

    # Check if there are any cycles
    try:
        cycle = nx.find_cycle(G)
        print(f"The graph contains a cycle: {cycle}. Therefore, it is NOT a valid MST.")
        return False
    except nx.exception.NetworkXNoCycle:
        pass  # No cycle found

    # Check the number of edges
    num_nodes = len(G.nodes)
    num_edges = len(G.edges)

    if num_edges != num_nodes - 1:
        print(f"The graph has {num_nodes} nodes and {num_edges} edges, which is incorrect for a valid MST.")
        return False

    print("The graph is a valid MST.")
    return True


def convert_to_cytoscape_format(input_file, output_file, epsilon=0.0001, source_col=None, target_col=None, weight_col=None):
    """
    Convert the MST file to Cytoscape-compatible format by adjusting edge weights.
    """
    # Load the original MST CSV file
    df = pd.read_csv(input_file)

    # If user hasn't specified columns, assume the first three columns are used
    if source_col is None or target_col is None or weight_col is None:
        if df.shape[1] < 3:
            raise ValueError("Input file must have at least three columns.")
        
        # Automatically detect the first three columns as Source, Target, Weight
        source_col, target_col, weight_col = df.columns[:3]
        print(f"Auto-detected columns: Source = {source_col}, Target = {target_col}, Weight = {weight_col}")

    # Check if the specified columns exist
    if not all(col in df.columns for col in [source_col, target_col, weight_col]):
        raise ValueError(f"Columns {source_col}, {target_col}, {weight_col} not found in the input file.")

    # Extract necessary columns
    df = df[[source_col, target_col, weight_col]]
    df.columns = ['Source', 'Target', 'Weight']

    # Convert weight to float
    df['Weight'] = df['Weight'].astype(float)
    
    # Generate the edge list for MST validation
    edges = list(zip(df['Source'], df['Target'], df['Weight']))
    
    # Validate the MST
    if not validate_mst(edges):
        print("MST validation failed. Exiting the program.")
        return

    # Adjust edge weights for Cytoscape compatibility
    df['Adjusted_Weight'] = 1 / (df['Weight'] + epsilon)
    df = df[['Source', 'Target', 'Adjusted_Weight']]
    df.columns = ['Source', 'Target', 'Weight']

    # Save to Cytoscape-compatible file
    df.to_csv(output_file, index=False)
    print(f"Successfully saved as {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate MST and Convert to Cytoscape-compatible format.")
    parser.add_argument("-i", "--input", required=True, help="Input MST CSV file (e.g., mst.csv)")
    parser.add_argument("-o", "--output", required=True, help="Output file name (e.g., cytoscape_network.csv)")
    parser.add_argument("-e", "--epsilon", type=float, default=0.0001, help="Small value to prevent division by zero (default: 0.0001)")
    parser.add_argument("--source_col", type=str, help="Name of the source column (default: first column)")
    parser.add_argument("--target_col", type=str, help="Name of the target column (default: second column)")
    parser.add_argument("--weight_col", type=str, help="Name of the weight column (default: third column)")

    args = parser.parse_args()

    convert_to_cytoscape_format(args.input, args.output, args.epsilon, args.source_col, args.target_col, args.weight_col)
