# CDST
CoDing Sequence Typer (CDST) is a simple, efficient, decentralized, and easily shareable genome typing and clustering method similar to cg/wgMLST, based on MD5 hash mapping of coding sequences (CDS) from genome assemblies.

----------------------------------------------------
DEPENDENCIES
----------------------------------------------------
Before running CDST, ensure that the following dependencies are installed:

Python Packages:
- argparse
- hashlib
- json
- pandas
- biopython
- networkx
- scipy

Install them using:
$ pip install biopython pandas networkx scipy

----------------------------------------------------
INSTALLATION
----------------------------------------------------
Clone this repository and navigate into the project folder:
$ git clone https://github.com/l1-mh/cdst.git
$ cd cdst

Make the script executable:
$ chmod +x cdst.py

Alternatively, you can run it directly using Python:
$ python cdst.py --help

----------------------------------------------------
USAGE
----------------------------------------------------
CDST provides multiple subcommands for different analysis steps.

1. Generate MD5 Hashes from FASTA Files:
$ python cdst.py generate -i sample1.ffn sample2.ffn -o output/

2. Compute Genetic Distance Matrices:
$ python cdst.py matrix -j output/md5_hashes.json -o output/

3. Generate Minimum Spanning Tree (MST):
$ python cdst.py mst -m output/difference_matrix.csv -o output/

4. Generate Hierarchical Clustering Tree:
$ python cdst.py hc -m output/difference_matrix.csv -o output/

Run the Full Pipeline Above:
$ python cdst.py run -i sample1.ffn sample2.ffn -o output/ -T both

Merge Databases:
$ python cdst.py join -d dir1/ dir2/ -o merged_output/ --matrix --mst

Compare New Samples Against an Existing Dataset:
$ python cdst.py test -i new_sample.ffn -j output/md5_hashes.json -o output/

----------------------------------------------------
OUTPUT FILES
----------------------------------------------------
Depending on the commands used, the following files will be generated:
- md5_hashes.json → Stores MD5 hashes for CDS sequences.
- comparison_matrix.csv → Number of shared hashes between samples.
- difference_matrix.csv → Distance matrix based on hash differences.
- edge_list.csv → Edge list representation of pairwise distances.
- mst.csv → Minimum Spanning Tree (MST) edge list.
- hc.newick → Hierarchical Clustering tree in Newick format.
- comparison_results.csv → Closest matches for new samples.

