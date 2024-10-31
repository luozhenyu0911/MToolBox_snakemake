the way to run MToolBox_snakemake:

1. cp -rf ./* /your/path/

2. you just need to change the following three files:
data/analysis.tab
data/datasets.tab
data/reference_genomes.tab

3. to sub. to run
yhbatch -N 1 -n 24 -p rhenv run1_MToolBox-variant-calling.sh

