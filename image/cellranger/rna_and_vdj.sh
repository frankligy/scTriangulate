#BSUB -W 24:00
#BSUB -M 132000
#BSUB -n 20
#BSUB -R "span[hosts=1]"
#BSUB -J snaf_sc
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err

# module load sratoolkit/2.10.4
# cat tcr_srr_list.txt | xargs -L 1 -P 5 fasterq-dump -e 20
# for file in *.fastq; do echo $file; done | xargs -L 1 -P 10 gzip

# cd reference
# wget https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
# tar -zxvf refdata-cellranger-GRCh38-3.0.0.tar.gz

# module load cellranger/6.0.1

# function run_cellranger_count() {
# cellranger count --id=$1 \
#                  --fastqs=/data/salomonis2/FASTQs/NCI-R01/SNAF_sc/rna_fastq \
#                  --transcriptome=/data/salomonis2/FASTQs/NCI-R01/SNAF_sc/reference/refdata-cellranger-GRCh38-3.0.0 \
#                  --sample=$1 \
#                  --jobmode=lsf \
#                  --localmem=132 \
#                  --maxjobs=20
# }

# declare -a index_array
# index_array+=(BC09_TUMOR2 BC10_TUMOR1 BC11_TUMOR1 BC11_TUMOR2)
# for item in ${index_array[@]}
# do
#     run_cellranger_count $item
# done

# # genotyping
# export PATH=$(pwd)/hla_genotying:$PATH
# function run_schlacount(){
# gunzip -c /data/salomonis2/FASTQs/NCI-R01/SNAF_sc/cellranger_rna/$1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > /data/salomonis2/FASTQs/NCI-R01/SNAF_sc/cellranger_rna/$1/outs/barcodes.tsv
# sc_hla_count_linux -b /data/salomonis2/FASTQs/NCI-R01/SNAF_sc/cellranger_rna/$1/outs/possorted_genome_bam.bam \
#                    -d /data/salomonis2/FASTQs/NCI-R01/SNAF_sc/hla_genotying/reference \
#                    -c /data/salomonis2/FASTQs/NCI-R01/SNAF_sc/cellranger_rna/$1/outs/barcodes.tsv \
#                    --out-dir /data/salomonis2/FASTQs/NCI-R01/SNAF_sc/cellranger_rna/$1/outs/hla-typer-results \
#                    --pl-tmp /data/salomonis2/FASTQs/NCI-R01/SNAF_sc/cellranger_rna/$1/outs/pseudoaligner_tmp
# }
# declare -a index_array
# index_array+=(BC09_TUMOR1 BC09_TUMOR2 BC10_TUMOR1 BC11_TUMOR1 BC11_TUMOR2)
# for item in ${index_array[@]}
# do
#     run_schlacount $item
# done


# cd reference
# wget https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
# tar -zxvf refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz


module load cellranger/6.0.1
function run_cellranger_vdj() {
cellranger vdj --id=$1 \
               --fastqs=/data/salomonis2/FASTQs/NCI-R01/SNAF_sc/tcr_fastq \
               --reference=/data/salomonis2/FASTQs/NCI-R01/SNAF_sc/reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 \
               --sample=$1 \
               --localmem=132 \
               --localcores=20
}

declare -a index_array
index_array+=(BC09_TUMOR2 BC10_TUMOR1 BC11_TUMOR1 BC11_TUMOR2)
for item in ${index_array[@]}
do
    run_cellranger_vdj $item
done










































