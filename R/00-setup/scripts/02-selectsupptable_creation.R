# the following is for processing and creating some of the supplementary tables
# that have outputs from various bioinformatic tools. This is just for demonstration
# and the raw files are not included in the GitHub but are available upon request
# from authors. 

# EukCC Supplementary table
eukcc_files <- list.files("eukcc/", full.names = T)
eukcc_fileid <- list.files("eukcc") %>% word(sep = fixed("_"))

eukcc_list <- vector("list", length = length(eukcc_files))
names(eukcc_list) <- eukcc_fileid

for (i in 1:length(eukcc_files)){
  eukcc_df <- read.csv(eukcc_files[i], sep = "\t")
  eukcc_df$sample_id <- names(eukcc_list)[i]
  eukcc_df$sample_bin_id <- paste0(eukcc_df$sample_id, "_", eukcc_df$bin)
  eukcc_list[[i]] <- eukcc_df
}

eukcc_df_all <- do.call(rbind, eukcc_list)

write.csv(eukcc_df_all,"eukcc.csv", quote = F, row.names = F)

# Distance to TE probability testing permutations
te_probability_df <- read.csv("permutations/permtest_results.txt", sep = "\t") %>% na.omit()
te_probability_df$sig_alt <- te_probability_df$rand_alt
te_probability_df[te_probability_df$rand_pval > 0.05,]$sig_alt <- "not"

te_probability_df_test <- te_probability_df[te_probability_df$sample_id == "glLepBurg3_bin.203",]
te_probability_df_test2 <- te_probability_df[te_probability_df$sample_id == "glPelHori1_bin.14",]


ggplot(te_probability_df_test, aes(x = rand_zscore, y = -log(rand_pval), color = rand_alt))+
  geom_point()
ggplot(te_probability_df_test, aes(x = rand_zscore, y = -log(rand_pval), color = sig_alt))+
  geom_point()


ggplot(te_probability_df_test[te_probability_df_test$sig_alt != "not",], aes(x = rand_zscore, y = description, size = -log(rand_pval), color = sig_alt))+
  geom_point(alpha = 0.6)+theme_bw()
ggplot(te_probability_df_test2, aes(x = rand_zscore, y = description, size = -log(rand_pval), color = rand_alt))+
  geom_point(alpha = 0.6)+theme_bw()

# bbmap genome stats
bbmap_genomestats_files <- list.files("genome_stats/", full.names = T)
bbmap_genomestats_genomeid <- list.files("genome_stats/")

bbmap_genomestats_list <- vector("list", length = length(bbmap_genomestats_genomeid))
names(bbmap_genomestats_list) <- bbmap_genomestats_genomeid

for (i in 1:length(bbmap_genomestats_list)){
  bbmap_genomestats <- read.csv(bbmap_genomestats_files[i], sep = "\t")
  sampleid <- bbmap_genomestats_genomeid[i]
  contig_num <- bbmap_genomestats[3,2]
  genomesize <- bbmap_genomestats[5,2]
  gc_content <- bbmap_genomestats[1,8]
  contig_n50 <- word(bbmap_genomestats[7,2], 2, sep =fixed("/")) 
  contig_l50 <- word(bbmap_genomestats[7,2], sep = fixed("/"))
  contig_n90 <- word(bbmap_genomestats[9,2],2, sep = fixed("/"))
  contig_l90 <- word(bbmap_genomestats[9,2], sep = fixed("/"))
  bbmap_genomestats_list[[i]] <- data.frame("bin_id" = sampleid, "size" = genomesize,
                                            "contig_number" = contig_num,
                                            "gc_content" = gc_content, "contig_l50" = contig_l50,
                                            "contig_n50" = contig_n50, "contig_l90" =contig_l90,
                                            "contig_n90" = contig_n90)
}

bbmap_genomestats_df <- do.call(rbind, bbmap_genomestats_list)

write.csv(bbmap_genomestats_df, "bbmap_genomestats.csv",
          row.names = F, quote = F)
# bbmap assembly stats
bbmap_assemblystats_files <- list.files("assembly_stats/", full.names = T)
bbmap_assemblystats_genomeid <- list.files("assembly_stats/")

bbmap_assemblystats_list <- vector("list", length = length(bbmap_assemblystats_genomeid))
names(bbmap_assemblystats_list) <- bbmap_assemblystats_genomeid

for (i in 1:length(bbmap_assemblystats_list)){
  bbmap_genomestats <- read.csv(bbmap_assemblystats_files[i], sep = "\t")
  sampleid <- bbmap_assemblystats_files[i]
  contig_num <- bbmap_genomestats[3,2]
  genomesize <- bbmap_genomestats[5,2]
  gc_content <- bbmap_genomestats[1,8]
  contig_n50 <- word(bbmap_genomestats[7,2], 2, sep =fixed("/")) 
  contig_l50 <- word(bbmap_genomestats[7,2], sep = fixed("/"))
  contig_n90 <- word(bbmap_genomestats[9,2],2, sep = fixed("/"))
  contig_l90 <- word(bbmap_genomestats[9,2], sep = fixed("/"))
  bbmap_assemblystats_list[[i]] <- data.frame("bin_id" = sampleid, "size" = genomesize,
                                            "contig_number" = contig_num,
                                            "gc_content" = gc_content, "contig_l50" = contig_l50,
                                            "contig_n50" = contig_n50, "contig_l90" =contig_l90,
                                            "contig_n90" = contig_n90)
}

bbmap_assemblystats_df <- do.call(rbind, bbmap_assemblystats_list)

write.csv(bbmap_assemblystats_df, "bbmap_assemblystats.csv",
          row.names = F, quote = F)
