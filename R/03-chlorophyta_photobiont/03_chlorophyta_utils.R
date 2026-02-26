library(circlize)
library(Biostrings)
library(universalmotif)
library(tidyverse)
library(ape)
library(ggtree)
library(aplot)
library(tximportData)
library(tximport)
library(rtracklayer)
library(edgeR)
library(topGO)

te_colours <- c("DNA" = "#0f6385",
  "LINE" = "#358292",
  "LTR" = "#88aca2",
  "Other (Simple Repeat, Microsatellite, RNA)" = "#ecd9af",
  "Penelope" = "#d6a379",
  "Rolling Circle" = "#a77c65",
  "SINE" = "#745644",
  "Unclassified" = "#573d31")

alt12palette <- c("#372D2C",
                "#562c29",
                 "#ab5852",
                 "#cb9979",
                 "#eadaa0",
                 "#d69e49",
                "#898152",
                "#849949",
                 "#838469",
                 "#657268",
                 "#476066",
                "#272F30")

kimura_df <- function(divsumpath, species, metadata){
  dist <- as.data.frame(t(fread(divsumpath, skip = "Div", header = TRUE)))
  dist$species <- species
  dist$type <- rownames(dist)
  dist <- dist[! dist$type == "X",]
  colnames(dist) <- c(dist[1,1:71], "species", "type")
  dist <- dist[! dist$type == "Div",]
  rownames(dist) <- 1:nrow(dist)
  
  dist$classif <- dist$type
  dist$classif <- gsub("^DNA.*", "DNA", dist$classif)
  dist$classif <- gsub("^LTR.*", "LTR", dist$classif)
  dist$classif <- gsub(".*Penelope|^PLE.*", "Penelope", dist$classif)
  dist$classif <- gsub("^LINE.*", "LINE", dist$classif)
  dist$classif <- gsub("^SINE.*", "SINE", dist$classif)
  dist$classif <- gsub("^RC.*", "Rolling Circle", dist$classif)
  dist$classif <- gsub("^Unknown.*|Retroposon.*|Unspecified.*", "Unclassified", dist$classif)
  dist$classif <- gsub(".*RNA.*|^Satellite.*|^Simple_repeat.*|^Low_complexity.*|^ARTEFACT.*|^Other.*|^Segmenta.*", "Other (Simple Repeat, Microsatellite, RNA)", dist$classif)
  
  piv <- dist[,c(1:72, 74)] %>% pivot_longer(!c("species", "classif"), names_to = "Divergence", values_to = "count") %>% group_by(species, classif, Divergence) %>% summarise(Count = sum(count))
  piv$Divergence <- as.numeric(piv$Divergence)
  
  gen <- metadata[metadata$accession == species,]$assembly_length
  
  genomeSizes <- as.data.frame(matrix(c(species, gen),
                                      byrow = TRUE, ncol = 2))
  
  colnames(genomeSizes) <- c("species", "genomeSize")
  piv <- merge(piv, genomeSizes)
  
  piv$genomeSize <- as.numeric(piv$genomeSize)
  
  piv$proportion <- piv$Count / piv$genomeSize

  for (i in 1:nrow(piv)) {
    if (piv$Count[i] < 0) {
      piv$Count[i] <- 0
    }
  }
  
  for (i in 1:length(piv$proportion)){
    if (piv$proportion[i] < 0) {
      piv$proportion[i] <- 0
    }
    if (i %% 100 == 0) { #if x is divisible by 100
      cat(sprintf("Completed: %s\r", i)) #print Completed: %s\r
    }
  }
  
  piv <- piv %>% distinct()
  piv2 <- piv %>% group_by(species, Divergence, classif) %>% summarise(proportion = sum(proportion))
  piv2
}

make_circos_plot <- function(fasta_file, gff_file, gff_te_file){
  fasta <- readDNAStringSet(fasta_file)
  fasta_bed <- data.frame("seq" = names(fasta), "start" = 1, end = width(fasta)) %>% arrange(desc(end))
  contigcount <- length(fasta)
  fasta_newnames <- structure(paste0("ctg_", seq(from = 1, to = contigcount)), names = fasta_bed$seq)
  fasta_bed$seq <- fasta_newnames[fasta_bed$seq]
  
  g <- as.data.frame(get_bkg(fasta, k = 1, RC = F, window = TRUE, window.size = 10000, merge.res = FALSE))
  gC <- g[g$klet == "C",]
  gG <- g[g$klet == "G", ]
  gCG <- data.frame(seq = gC$sequence, start = gC$start, end = gC$stop, score = (gC$probability + gG$probability)/2)
  gCG$seq <- fasta_newnames[gCG$seq]
  
  gff_colnames <- c("seq","source","feature","start","end","score","strand","frame","attribute")
  
  gff <- read.table(gff_file)
  colnames(gff) <- gff_colnames
  gff_genes <- gff[gff$feature %in% c("gene"),]
  gff_genes$seq <- fasta_newnames[gff_genes$seq]
  gene_bed <- gff_genes[,c(1,4,5,6)]
  
  te_gff <- read.table(gff_te_file)
  colnames(te_gff) <- gff_colnames
  te_gff$seq <- fasta_newnames[te_gff$seq]
  te_bed <- te_gff[,c(1,4,5,6)]
  
  circos.clear()
  circos.par(start.degree = 100)
  circos.genomicInitialize(fasta_bed, plotType = "axis")
  circos.track(ylim = c(0, 1), 
               bg.col = c(rep(c("#afa59f", "white"),contigcount)), 
               bg.border = "grey20", track.height = 0.05)
  circos.genomicDensity(te_bed, col = c("#0f6385"), track.height = 0.2)
  circos.genomicDensity(gene_bed, col = c("#cd8157"), track.height = 0.2)
  circos.genomicTrack(gCG, panel.fun = function(region, value, ...) {
    i = getI(...)
    circos.genomicLines(region, value, col = "grey65", border = "grey65", area = TRUE, track.height = 0.1,...)
  })
  circos.clear()
}

# salmon RNAseq stuff
get_tx2gene <- function(x){
  gene2tx <- read.csv(x, sep = "\t", header = FALSE, col.names = c("gene", "tx"))
  tx2gene <- gene2tx[,c(2,1)]
  tx2gene
}

get_tximport <- function(salmonfile, genemapfile, txOut = FALSE){
  tx2gene <- get_tx2gene(genemapfile)
  tximport <- tximport(salmonfile, type = "salmon", tx2gene = tx2gene, txOut = txOut)
  tximport
}
bulktximport <- function(salmondir, genemapdir, txOut = FALSE){
  speciesnames <- list.files(salmondir)
  
  # samplepath <- file.path(salmondir, speciesnames,"transcripts_quant",samplenames)
  
  genemapfiles <- list.files(genemapdir)
  genemapid <- strsplit(genemapfiles, ".", fixed = TRUE)
  genemapfilepaths <- file.path(genemapdir, genemapfiles)
  
  tximport_list <- vector("list", length = length(genemapfiles))
  names(tximport_list) <- speciesnames
  
  for (i in 1:length(tximport_list)){
    speciesid <- genemapid[[i]][1]
    
    samplenames <- list.files(file.path(salmondir,speciesnames[i],"transcripts_quant"))
    samplepath <- file.path(salmondir, speciesnames[i],"transcripts_quant",samplenames)
    
    salmonfilepaths <- file.path(samplepath,"quant.sf")
    names(salmonfilepaths) <- samplenames
    
    # print(salmonfilepaths)
    # print(genemapfilepaths)
    tximport_list[[i]] <- get_tximport(salmonfilepaths, genemapfilepaths[i], txOut)
  }
  tximport_list
}

calc_topgo_fisher <- function(markergenes, ontology){
  markerList <- names(markergenes)
  testList <- factor(as.integer(testGenes %in% markerList))
  names(testList) <- testGenes
  
  testData <- new("topGOdata", ontology = ontology,
                  allGenes = testList, annot = annFUN.gene2GO,
                  gene2GO = testGene2GO)
  
  testFisher <- new("classicCount", testStatistic = GOFisherTest,
                    name = "fisher")
  testRes <- getSigGroups(testData, testFisher)
  testTab <- GenTable(testData, FisherPval = testRes,
                      orderBy = "FisherPval", topNodes = 1000)
  testTab
}

#bincat_colours_trebo_old <- c("20000" = "#0B0405FF",
#                    "10001d20000" = "#382A54FF",
#                    "5001d10000" = "#395D9CFF",
#                    "1001d5000" = "#3497A9FF",
#                    "1d1000" = "#60CEACFF",
#                    "0d" = "#B2EDBF")

bincat_colours_trebo <- c("20000" = "#440154FF",
                          "10001d20000" = "#414487FF",
                          "5001d10000" = "#2A788EFF",
                          "1001d5000" = "#22A884FF",
                          "1d1000" = "#7AD151FF",
                          "0d" = "#FDE725FF")

bincat_expression_colours_trebo <- c("20000_0" = "#9A79A3","20000_>0" = "#440154FF",
                               "10001d20000_0" = "#7E80AB","10001d20000_>0" = "#414487FF",
                               "5001d10000_0" ="#6F99A3","5001d10000_>0" = "#2A788EFF",
                               "1001d5000_0" = "#BCE3D9","1001d5000_>0" = "#22A884FF",
                               "1d1000_0" = "#C1E8B0","1d1000_>0" = "#7AD151FF",
                               "0d_0" = "#FCF29A","0d_>0" = "#FDE725FF")





