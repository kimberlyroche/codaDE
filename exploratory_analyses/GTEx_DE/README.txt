`mart_export.txt` contains the Ensembl IDs of protein coding genes output from: http://useast.ensembl.org/biomart/martview/52de929db589b1609871b7be8e94feb6
Using settings:
  Human genes (GRCh38.p13)
  Filters: protein_coding
  Attributes: Gene stable ID
Note: Had to remove the spaces in the header labels

Use GTEx to determine the approximate proportion of genes differentially expressed between a random pair
of tissues. Print this to STDOUT.

Usage:
(1) Run sbatch file
(2) Move *.out to slurm_output
(3) Run parse.pl
