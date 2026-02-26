txt_from_folder <- function(path, ext = ".fa"){
  files <- list.files(path)
  filepaths <- paste0(path,files)
  filelist <- vector("list",length(files))
  for (i in 1:length(files)){
    table <- read.table(filepaths[i], sep = "\t", header = TRUE)
    filename <- strsplit(files[i], ".high")[[1]][1]
    table$sample_id <- paste0(filename, ext)
    filelist[[i]] <- table
  }
  table_final <- do.call(rbind, filelist)
  table_final
}

fastastats_from_folder <- function(path){
  files <- list.files(path)
  filepaths <- paste0(path,files)
  filelist <- vector("list",length(files))
  for (i in 1:length(files)){
    seq <- readDNAStringSet(filepaths[i])
    table <- data.frame(sample_id = files[i], length = sum(width(seq)), contig_count = length(seq))
    filelist[[i]] <- table
  }
  table_final <- do.call(rbind, filelist)
  table_final
}

gff_from_folder <- function(path){
  files <- list.files(path)
  filepaths <- paste0(path,files)
  filelist <- vector("list",length(files))
  for (i in 1:length(files)){
    #print(filepaths[i])
    gff <- read.table(filepaths[i], header = FALSE)
    colnames(gff) <- c("seq","source","feature","start","end","score","strand","frame","attribute")
    gff$sample <- gsub("[.]gff","",files[i])
    filelist[[i]] <- gff
  }
  table_final <- do.call(rbind, filelist)
  table_final
}

get_gene2te_dist_df <- function(gene_gff, te_gff, type = "transcript"){
  genes <- rtracklayer::import(gene_gff)
  genes <- genes[as.character(genes$type) == type]
  tes <- rtracklayer::import(te_gff)
  tes <- tes[width(tes) > 100]
  tes$ID2 <- paste0("TE_", 1:length(tes))
  dd <- distanceToNearest(genes, tes, ignore.strand = TRUE)
  genes$NearestTE <- NA
  genes$NearestTEDist <- NA
  genes$NearestTEDist_class <- NA
  genes$NearestTE[queryHits(dd)] <- tes$ID2[subjectHits(dd)]
  genes$NearestTEDist[queryHits(dd)] <- mcols(dd)$distance
  
  genes$NearestGene <- NA
  gd <- distanceToNearest(genes, ignore.strand = TRUE)
  genes$NearestGene[queryHits(gd)] <- mcols(gd)$distance
  #print(head(tes$type))
  #print(tes$type[subjectHits(dd)])
  genes$NearestTEDist_class[queryHits(dd)] <- as.character(tes$type)[subjectHits(dd)]
  genes_df <- as.data.frame(genes)
  genes_df
}

gene2te_bulk <- function(genedirpath, tedirpath){
  genefiles <- list.files(genedirpath)
  genefilepaths <- paste0(genedirpath,genefiles)
  genefilelist <- vector("list",length(genefiles))
  
  tefiles <- list.files(tedirpath)
  tefilepaths <- paste0(tedirpath,tefiles)
  tefilelist <- vector("list",length(tefiles))
  
  gene2te_list <- vector("list", length(genefiles))
  
  
  for (i in 1:length(genefiles)){
    for (j in 1:length(tefiles)){
      genesample <- gsub("[.]gtf", "",genefiles[i])
      #genesample <- gsub("[.].+","",genesample)
      
      tesample <- gsub("[.]filteredRepeats[.]gff","",tefiles[j])
      #tesample <- gsub("_.+","",tesample)
      #tesample <- gsub("^gl","", tesample)
      
      print(genesample)
      print(tesample)
      
      if (genesample == tesample){
        gene2te_df <- get_gene2te_dist_df(genefilepaths[i],tefilepaths[j])
        gene2te_df$sample_id <- gsub("[.]gtf","",genefiles[i])
        gene2te_list[[i]] <- gene2te_df
      }
    }
  }
  gene2te_bulk_df <- do.call(rbind, gene2te_list)
  gene2te_bulk_df
}

gene2te_bulk_gff <- function(genedirpath, tedirpath){
  genefiles <- list.files(genedirpath)
  genefilepaths <- paste0(genedirpath,genefiles)
  genefilelist <- vector("list",length(genefiles))
  
  tefiles <- list.files(tedirpath)
  tefilepaths <- paste0(tedirpath,tefiles)
  tefilelist <- vector("list",length(tefiles))
  
  gene2te_list <- vector("list", length(genefiles))
  
  
  for (i in 1:length(genefiles)){
    for (j in 1:length(tefiles)){
      genesample <- gsub("[.]gff", "",genefiles[i])
      #genesample <- gsub("[.].+","",genesample)
      
      tesample <- gsub("[.]filteredRepeats[.]gff","",tefiles[j])
      #tesample <- gsub("_.+","",tesample)
      #tesample <- gsub("^gl","", tesample)
      
      print(genesample)
      print(tesample)
      
      if (genesample == tesample){
        gene2te_df <- get_gene2te_dist_df(genefilepaths[i],tefilepaths[j], "gene")
        gene2te_df$sample_id <- gsub("[.]gff","",genefiles[i])
        gene2te_list[[i]] <- gene2te_df
      }
    }
  }
  gene2te_bulk_df <- do.call(rbind, gene2te_list)
  gene2te_bulk_df
}

get_kofamscan_df <- function(kofampath, kofamseqnamespath){
  kofamscan_df <- read.csv(kofampath, sep = "\t", header = FALSE)
  colnames(kofamscan_df) <- c("kofam_contig_id", "KEGG")
  kofamscan_df <- kofamscan_df[kofamscan_df$KEGG != "",]
  
  kofam_seqnames <- read.csv(kofamseqnamespath, sep = "\t", header = FALSE)
  colnames(kofam_seqnames) <- c("metaeuk_id", "kofam_contig_id")
  
  targets <- c()
  contigs <- c()
  targetlist <- strsplit(kofam_seqnames$metaeuk_id, "|", fixed = TRUE)
  contiglist <- strsplit(kofam_seqnames$metaeuk_id, "|", fixed = TRUE)
  
  for (i in 1:nrow(kofam_seqnames)){
    targets <- c(targets, targetlist[[i]][1])
    contigs <- c(contigs, contiglist[[i]][2])
  }
  
  kofam_seqnames$contig_id <- contigs
  kofam_seqnames$Target_ID <- targets
  
  kofam_full_df <- merge(kofamscan_df, kofam_seqnames)
  kofam_full_df
}

kofamscan_bulk <- function(kofamdirpath, kofamseqnamesdirpath){
  kofamfiles <- list.files(kofamdirpath)
  kofamfilepaths <- paste0(kofamdirpath,kofamfiles)
  kofamfilelist <- vector("list",length(kofamfiles))
  
  knamesfiles <- list.files(kofamseqnamesdirpath)
  knamesfilepaths <- paste0(kofamseqnamesdirpath,knamesfiles)
  knamesfilelist <- vector("list",length(knamesfiles))
  
  kofam_full_list <- vector("list", length(kofamfiles))
  kofam_files_broken <- c()
  for (i in 1:length(kofamfilelist)){
    kofamsample <- gsub("_result[.]txt", "",kofamfiles[i])
    knamesample <- gsub("_kofam_contigrefs[.]tsv","",knamesfiles[i])
    
    if (kofamsample == knamesample){
      if (ncol(read.csv(kofamfilepaths[i], header = FALSE, sep = "\t")) == 2){
        kofam_full_df <- get_kofamscan_df(kofamfilepaths[i],knamesfilepaths[i])
        kofam_full_df$sample_id <- kofamsample
        
        kofam_full_list[[i]] <- kofam_full_df
      }
    }
  }
  #kofam_full_list
  kofam_full_df <- do.call(rbind, kofam_full_list)
  kofam_full_df
}

get_iproscan <- function(ipropath){
  ipro_df <- read.csv(ipropath, sep = "\t", header = FALSE)
  colnames(ipro_df) <- c("metaeuk_id", "md5", "length", "analysis",
                         "signature_accession","signature_description", "start", "end",
                         "score", "status", "run_date", "ipro_accesion","ipro_description",
                         "go_term","pathways")
  
  targets <- c()
  contigs <- c()
  seqlist <- strsplit(ipro_df$metaeuk_id, ".", fixed = TRUE)
  
  for (i in 1:nrow(ipro_df)){
    targets <- c(targets, seqlist[[i]][1])
    contigs <- c(contigs, seqlist[[i]][2])
  }
  
  ipro_df$contig_id <- contigs
  ipro_df$Target_ID <- targets
  
  ipro_df
}

iproscan_bulk <- function(iprodirpath){
  iprofiles <- list.files(iprodirpath)
  iprofilepaths <- paste0(iprodirpath,iprofiles)
  iprofilelist <- vector("list",length(iprofiles))
  
  ipro_full_list <- vector("list", length(iprofiles))
  
  for (i in 1:length(iprofilelist)){
    iprosample <- gsub("[.]faa[.]tsv", "",iprofiles[i])
    ipro_full_df <- get_iproscan(iprofilepaths[i])
    ipro_full_df$sample_id <- iprosample
    ipro_full_list[[i]] <- ipro_full_df
  }
  ipro_full_df <- do.call(rbind, ipro_full_list)
  ipro_full_df
}

bincat_expression_colours_old <- c("#26445d","#7491a9",
                               "#4d4376","#97677c",
                               "#007a7a","#83acac",
                               "#6d8a70","#94bc98",
                               "#6b5151","#9a8686",
                               "#ef894a","#edb758")

bincat_expression_colours <- c("20000_0" = "#A2A2CA","20000_>0" = "#121254",
"10001d20000_0" = "#A192B0","10001d20000_>0" = "#320a5e",
"5001d10000_0" ="#C896C2","5001d10000_>0" = "#781c6d",
"1001d5000_0" = "#BB7B89","1001d5000_>0" = "#bc3754",
"1d1000_0" = "#ECB193","1d1000_>0" = "#ed6925",
"0d_0" = "#F9D88F","0d_>0" = "#fbb61a")

bincat_colours <- c("20000" = "#0C0C34",
                    "10001d20000" = "#320a5e",
                    "5001d10000" = "#781c6d",
                    "1001d5000" = "#bc3754",
                    "1d1000" = "#ed6925",
                    "0d" = "#fbb61a")

calc_topgo_fisher <- function(markergenes, genelist, gene2go,ontology, topNodes = 100, pval = 0.05){
  #markergeneList <- rownames(markergenes[markergenes$p_val_adj <= 0.05,])
  testList <- factor(as.integer(genelist %in% markergenes))
  names(testList) <- genelist
  
  testData <- new("topGOdata", ontology = ontology,
                  allGenes = testList, annot = annFUN.gene2GO,
                  gene2GO = gene2go)
  
  testFisher <- new("classicCount", testStatistic = GOFisherTest,
                    name = "fisher")
  testRes <- getSigGroups(testData, testFisher)
  testTab <- GenTable(testData, FisherPval = testRes,
                      orderBy = "FisherPval", topNodes = topNodes)
  
  testGenes <- lapply(testTab$GO.ID,
    function(x) as.character(unlist(genesInTerm(object = testData, whichGO = x))))
  testTab$Genes <- unlist(lapply(testGenes,
    function(x) paste0(x[x %in% markergenes],
                       collapse = "|")))
  testTab <- testTab[testTab$FisherPval <= pval,]
  testTab
}

species_mag_codes <- c("glLepBurg3_bin.203" = "Leptogium burgessi", 
                       "glLobPulm2_bin.43" = "Lobaria pulmonaria", 
                       "glNepLaev11_merged.0" = "Nephroma leavigatum", 
                       "glPelHori1_bin.14" = "Peltigera horizontalis", 
                       "glPelHyme1_bin.206" = "Peltigera hymenina",
                       "glPelMemb1_bin.35" = "Peltigera membranacea",
                       "glPelPrae3_bin.50" = "Peltigera praetextata", 
                       "glPseNorv1_bin.13" = "Pseudocyphellaria norvegica",
                       "glRicVire11_merged.0" = "Ricasolia virens")

species_sort_codes <- c("LepBurg3" = "Leptogium burgessi", 
                       "LobPulm2" = "Lobaria pulmonaria", 
                       "NepLaev11" = "Nephroma leavigatum", 
                       "PelHori1" = "Peltigera horizontalis", 
                       "PelHyme" = "Peltigera hymenina",
                       "PelMemb" = "Peltigera membranacea",
                       "PelPrae" = "Peltigera praetextata", 
                       "PseNorv" = "Pseudocyphellaria norvegica",
                       "RicVire" = "Ricasolia virens")
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
  #print(piv)
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

te_colours <- c("DNA" = "#0f6385",
                "LINE" = "#358292",
                "LTR" = "#88aca2",
                "Other (Simple Repeat, Microsatellite, RNA)" = "#ecd9af",
                "Penelope" = "#d6a379",
                "Rolling Circle" = "#a77c65",
                "SINE" = "#745644",
                "Unclassified" = "#573d31")

tesub_colours <- c("DNA/TcMar-Fot1" = "#03a1e1",
                   "DNA/TcMar-Sagan" = "#035f84",
                "LINE/L1" = "#93d3e1",
                "LINE/Tad1" = "#358292",
                "LTR/Gypsy" = "#cdd6d3",
                "LTR/Copia" = "#88aca2",
                "Satellite" = "#ecd9af",
                "Unknown" = "#573d31")


cazy12_palette <- c("AA" = "#372D2C",
                  "CB" = "#562c29",
                  "CB,AA" = "#ab5852",
                  "CB,CE" = "#cb9979",
                  "CB,GH" = "#eadaa0",
                  "CB,GT" = "#d69e49",
                  "CE" = "#898152",
                  "GH" = "#849949",
                  "GH,CB" = "#838469",
                  "GT" = "#657268",
                  "GT,GH" = "#476066",
                  "PL" = "#272F30")

