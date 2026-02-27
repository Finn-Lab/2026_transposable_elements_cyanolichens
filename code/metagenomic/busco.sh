# example busco run - fungi
singularity exec containers/busco_v5.7.0_cv1.sif busco -i genomes -m genome -l fungi_odb10 -c 10  

# example busco run - chlorophyte
singularity exec containers/busco_v5.7.0_cv1.sif busco -i genomes -m genome -l chlorophyta_odb10 -c 10  