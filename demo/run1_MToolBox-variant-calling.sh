#!/bin/bash
# yhbatch -N 1 -n 24 -p rhenv
source /BIGDATA2/gzfezx_shhli_2/miniconda3/etc/profile.d/conda.sh
mtoolbox-activate
# ln -s /BIGDATA2/gzfezx_shhli_2/software/MToolBox_snakemake/modules .
MToolBox-variant-calling  -rkp --core 24 --scheduler greedy
