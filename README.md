# Transposable Element Diversification and the evolution of Peltigerales Lichen Symbionts
Scripts for analysis presented in Cameron et al. (2026) [Transposable Element Diversification and the Evolution of Peltigerales Lichen Symbionts](https://www.biorxiv.org/content/10.64898/2026.02.14.705750v1). 

## Abstract
Lichens are composite organisms formed through the symbiotic association between fungi, algae and/or bacteria. Multiple independent origins of the lichenized lifestyle have been reported in both fungal and algal lineages, but the molecular mechanisms and evolution underpinning these symbiotic relationships remain largely unknown. In this study, we performed long-read metagenomic sequencing on 11 Peltigerales lichen species to characterize the genomic content of lichen symbionts via metagenome assembled genomes (MAGs). Peltigerales genomes generated in this work represent the largest Lecanoromycetes genome sequenced to date, driven by high transposable element content. Transposable elements (TEs) are known to drive genome evolution in other symbioses but have been underexplored in lichen symbionts due technological limitations. Transcriptomics revealed that many genes associated with adaptations to the lichenized lifestyle are associated with TEs suggesting that they may play a key role in the evolution of lichenization.

## Code
Scripts for analysis of raw reads and generation of MAGs from long-read lichen metagenomes. 

## R 
Downstream analysis of results for bacterial, algal and fungal MAGs was completed in R for taxonomic profiling, functional enrichment analysis and deeper exploration of genomic characteristics including TE dynamics. Main areas featured in the Results section of Cameron et al., 2026 are highlighted below. 

Select intermediate and input files are provided in this repository. Please contact authors for access to other inputs and intermediate files found in the R analysis. 

[Microbial Community Analysis](https://github.com/Finn-Lab/2026_transposable_elements_cyanolichens/tree/main/R/01-microbial_community): Analysis of lichen metagenome composition and exploration of cyanobacterial photobionts.  
[Chlorophyte Photobiont](https://github.com/Finn-Lab/2026_transposable_elements_cyanolichens/tree/main/R/03-chlorophyta_photobiont): Analysis of chlorophyte photobiont MAG detected in _Lobaria pulmonaria_ and _Ricasolia virens_.  
[Mycobiont](https://github.com/Finn-Lab/2026_transposable_elements_cyanolichens/tree/main/R/04-fungi): Analysis of fungal MAGs and exploration of TE expression dynamics in lichenized fungi.   



