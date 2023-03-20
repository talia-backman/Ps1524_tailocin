# install VIBRANT to the conda environment 'annotate' and load that environment
# installation of VIBRANT can be found here: https://github.com/AnantharamanLab/VIBRANT
conda activate annotate
cd /uufs/chpc.utah.edu/common/home/karasov-group1/phages/genome_annotation/genome_quality/genomes
for i in *.fasta; do VIBRANT_run.py -i $i; done
