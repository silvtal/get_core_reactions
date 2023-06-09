## `core_functions_full.py`

### Create file with X% consensus annotations ("PCG pangenome"). Save to files parsed list in BiGG/KEGG reaction format.

This script generates a list of core functions for a given phylogenetic
core group (PCG). This list is considered that PCG's pangenome, defined
as annotations present in at least 90% of available genomes for each PCG.
The input for this tool is a path containing .tsv annotation files
generated with eggNOG-mapper, for which a consensus annotation file will
be generated. The core reactions can be retrieved in KEGG, EC and in BiGG
format. A file containing descriptions corresponding to EC and KEGG entries
is also generated.

_More reactions can be added manually using tools like MetaNetX or by creating custom ko/EC files_


## `Print_unique_functions.ipynb`

### For each PCG, print all pangenome reactions exclusive to itself.

This script reads the files (either BiGG or KEGG) from the first step and goes through each PCG, compareing its reactions with the reactions of other groups, and keeping track of the unique reactions for each group in a dictionary.

## `Print_all_descriptions.ipynb`

### Print all pangenome reactions for each PCG, exclusive or not, alongside a description.
