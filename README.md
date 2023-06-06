
The `core_functions_full.py` script generates a list of core functions for a given phylogenetic
core group (PCG). This list is considered that PCG's pangenome, defined
as annotations present in at least 90% of available genomes for each PCG.
The input for this tool is a path containing .tsv annotation files
generated with eggNOG-mapper, for which a consensus annotation file will
be generated. The core reactions can be retrieved in KEGG, EC and in BiGG
format. A file containing descriptions corresponding to EC and KEGG entries
is also generated.

_More reactions can be added manually using tools like MetaNetX or by creating custom ko/EC files_
