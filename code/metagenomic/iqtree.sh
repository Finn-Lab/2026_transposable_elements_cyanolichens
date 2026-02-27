# example iqtree
# eukaryotes used BUSCO phylogenomics output alignment; bacteria used GTDBtk alignments

iqtree -s ${alignment}.aln -T 8 -m TEST -mset LG,WAG,JTT -bb 1000