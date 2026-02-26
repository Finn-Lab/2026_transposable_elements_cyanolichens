library(ggtree)
library(aplot)
library(ggplot2)
library(tidyverse)

# input files not included in GitHub are available from authors upon request
fungi_busco_tree <- readRDS("00-setup/RDS/fungi_busco_ggtree.RDS")

fungi_earlgrey_high_summ_table <- txt_from_folder("earlgrey/fungi/summary_highlevelcounts/")
fungi_earlgrey_high_summ_table$accession <- fungi_earlgrey_high_summ_table$sample_id

fungiref_earlgrey_high_summ_table <- txt_from_folder("earlgrey/fungi_refs/summary_highlevelcounts/")
fungiref_earlgrey_high_summ_table$accession <- substr(fungiref_earlgrey_high_summ_table$sample_id, 1, 15)

fungi_earlgrey_high_summ_table <- rbind(fungi_earlgrey_high_summ_table, fungiref_earlgrey_high_summ_table)
fungi_metadata <- read.csv("fungi_meta_habitat.csv")
fungi_earlgrey_highsumm_metadata <- merge(fungi_earlgrey_high_summ_table, fungi_metadata)

fungi_earlgrey_highsumm_metadata_busco <- fungi_earlgrey_highsumm_metadata[fungi_earlgrey_highsumm_metadata$accession %in% fungi_busco_tree_tips,]
fungi_earlgrey_busco_table_summary <- fungi_earlgrey_highsumm_metadata_busco %>% group_by(sampleid) %>% summarize(prop = sum(proportion*100))

fungi_te_comp <- ggplot(fungi_earlgrey_highsumm_metadata_busco, aes(x = proportion*100, y = sampleid, fill = tclassif))+
  geom_bar(stat = "identity", color = "black")+
  geom_vline(xintercept = median(fungi_earlgrey_busco_table_summary$prop), linetype = 2, colour = "firebrick3", linewidth = 3.5)+
  theme_bw()+scale_x_continuous(expand = c(0,0), limits = c(0,75))+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), legend.position = "none")+
  scale_fill_manual(values = te_colours)+
  labs(fill = "TE Classification", x = "Relative Genome Content (%)")+
  theme(axis.text.x = element_text(size = 32), axis.title.x = element_text(size = 48))

fungi_metadata_assembly <- read.csv("fungi_metadata.csv", header = TRUE)
fungi_metadata_assembly_busco <- merge(fungi_metadata_assembly, fungi_metadata, by = "accession")
fungi_metadata_assembly_busco <- merge(fungi_metadata_assembly_busco, fungi_earlgrey_busco_table_summary)
fungi_metadata_assembly_busco$core_genome <- fungi_metadata_assembly_busco$assembly_length *((100-fungi_metadata_assembly_busco$prop)/100)
fungi_metadata_assembly_busco$te_genome <- fungi_metadata_assembly_busco$assembly_length*(fungi_metadata_assembly_busco$prop / 100)
fungi_metadata_assembly_busco_v2 <- pivot_longer(fungi_metadata_assembly_busco, cols = c("core_genome", "te_genome"), names_to = "genometype", values_to = "genome_prop")
fungi_metadata_assembly_busco_v2[fungi_metadata_assembly_busco_v2$genometype == "core_genome",]$genometype <- fungi_metadata_assembly_busco_v2[fungi_metadata_assembly_busco_v2$genometype == "core_genome",]$life.y
fungi_metadata_assembly_busco_v2$genometype <- factor(fungi_metadata_assembly_busco_v2$genometype,
                                                      levels = c("te_genome", "lichen (me)", "lichen", "not lichenized (me)", "not lichenized")) 
fungi_tree_colours_tetype <- c(fungi_tree_colours,"te_genome" = "#d2c9ba")
fungi_genomesize_bar <- ggplot(fungi_metadata_assembly_busco_v2, aes(x = as.numeric(genome_prop), y = sampleid, fill = genometype))+
  geom_bar(stat = "identity", colour = "black")+
  geom_vline(xintercept = median(fungi_metadata_assembly_busco$assembly_length), linetype = 2, colour = "firebrick3", linewidth = 3.5)+
  theme_bw()+scale_x_continuous(expand = c(0,0), limits = c(0, 150000000))+
  scale_fill_manual(values =fungi_tree_colours_tetype)+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), legend.position = "none")+
  labs(x = "Genome Size (bp)")+
  theme(axis.text.x = element_text(size = 32), axis.title.x = element_text(size = 48))

#fungi_te_tree_genomes <- fungi_genomesize_bar %>% insert_left(fungi_busco_tree) %>% insert_right(fungi_te_comp)
fungi_te_tree_genomes <- fungi_genomesize_bar %>% insert_left(fungi_tree) %>% insert_right(fungi_te_comp)

ggsave("figures/fungi/fungitree_tes_v2.pdf", fungi_te_tree_genomes, 
       dpi = "retina", height = 45, width = 45)
kimuradir <- "earlgrey/fungi/kimura_distances/"
divsumfiles <- list.files(kimuradir)
divsum_accessions <- paste0(substr(divsumfiles, 1,(nchar(divsumfiles)-7)),".fa")

kimuradf_list <- vector("list", length(divsumfiles))
names(kimuradf_list) <- divsum_accessions
for (i in 1:length(divsum_accessions)){
  kimuradf_list[[i]] <- kimura_df(paste0(kimuradir,divsumfiles[i]),divsum_accessions[i], fungi_metadata_assembly)
}

kimuradf_list_busco <- kimuradf_list[fungi_metadata$accession]

kimura_plot_list <- vector("list", length(kimuradf_list))
names(kimura_plot_list) <- divsum_accessions

for (i in 1:length(divsum_accessions)){
  plot_df <- kimuradf_list_busco[[divsum_accessions[i]]]
  limits <- plot_df %>% group_by(Divergence) %>% summarise(perc = sum(proportion * 100))
  ylimit <- max(limits$perc) * 1.75
  kimura_plot_list[[i]] <- ggplot(plot_df[! plot_df$classif %like% "Other",], aes(x = Divergence, y = proportion*100, fill = classif)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylimit))+theme_void()+scale_fill_manual(values = te_colours)+
    theme(legend.position = "none")
}

allkimura_void_plots <- ggarrange(plotlist = kimura_plot_list, ncol = 1)
## all samples with EarlGrey v.4.1.0 except these two which use v.5.1.0 for patch
## GCF_003290485.1,GCA_010093895.1

####
primaryfungi <- c("glPelHori1_bin.14","glPelMemb1_bin.35", "glPseNorv1_bin.13")
combinedfungi <- c("glLepBurg3_bin.203", "glLobPulm2_bin.43", "glNepLaev11_merged.0",
                   "glPelHyme1_bin.206", "glPelPrae3_bin.50", "glRicVire11_merged.0")

gene2te_dist_primaryfun <- gene2te_bulk("stringtie_final_results_v3/", "earlgrey/fungi/te_gff/") %>% drop_na(NearestTE)
saveRDS(gene2te_dist_primaryfun,"gene2te_dist_df_v3.RDS")


gene2te_dist_otherfungi <- gene2te_bulk_gff("earlgrey/fungi/metaeuk_gff/", "earlgrey/fungi/te_gff/") %>% drop_na(NearestTE) 
gene2te_dist_otherfungi <- gene2te_dist_otherfungi[gene2te_dist_otherfungi$sample_id %notin% c("glRicVire11_merged.1", "glRicVire11_merged.2", primaryfungi,combinedfungi),]

gene2te_functionquant <- readRDS("data/gene2te_annotations/GeneToTEsFunctionsQuant_v3_PrimaryFungi.RDS")
gene2te_functionquant$RNA <- (gene2te_functionquant$Rep1 + gene2te_functionquant$Rep2) / 2
gene2te_functionquant <- gene2te_functionquant[order(gene2te_functionquant$RNA, decreasing = TRUE),]
gene2te_functionquant_v2 <- gene2te_functionquant[!duplicated(gene2te_functionquant$transcript_id), ]
gene2te_functionquant_v2$expression_status <- "0"
gene2te_functionquant_v2[!is.na(gene2te_functionquant_v2$RNA) & gene2te_functionquant_v2$RNA >= 1,]$expression_status <- ">0"

te_bin_cats <-  c (0, 1, 1000, 5000, 10000, 20000, Inf)
gene2te_functionquant_v2$bin_cat <- NA
gene2te_functionquant_v2[gene2te_functionquant_v2$NearestTEDist <= 1,]$bin_cat <- "0d"
gene2te_functionquant_v2[gene2te_functionquant_v2$NearestTEDist >= 2 & gene2te_functionquant_v2$NearestTEDist <= 1000,]$bin_cat <- "1d1000"
gene2te_functionquant_v2[gene2te_functionquant_v2$NearestTEDist >= 1001 & gene2te_functionquant_v2$NearestTEDist <= 5000,]$bin_cat <- "1001d5000"
gene2te_functionquant_v2[gene2te_functionquant_v2$NearestTEDist >= 5001 & gene2te_functionquant_v2$NearestTEDist <= 10000,]$bin_cat <- "5001d10000"
gene2te_functionquant_v2[gene2te_functionquant_v2$NearestTEDist >= 10001 & gene2te_functionquant_v2$NearestTEDist <= 20000,]$bin_cat <- "10001d20000"
gene2te_functionquant_v2[gene2te_functionquant_v2$NearestTEDist > 20000,]$bin_cat <- "20000"
gene2te_functionquant_v2$species_name <-  sub('\\..*', '',gene2te_functionquant_v2$oId)
gene2te_bin_counts <- gene2te_functionquant_v2 %>% group_by(bin_cat, species_name, expression_status, NearestTEDist_class) %>% summarize(count = n())

gene2te_bin_counts$bin_cat <- factor(gene2te_bin_counts$bin_cat, 
                                           levels = c("0d", "1d1000",
                                                      "1001d5000","5001d10000",
                                                      "10001d20000","20000"))
gene2te_bin_counts$bin_cat_expression <- paste0(gene2te_bin_counts$bin_cat, "_",gene2te_bin_counts$expression_status)

gene2te_bin_counts$bin_cat_expression <- factor(gene2te_bin_counts$bin_cat_expression, 
                                     levels = c("20000_0","20000_>0",
                                                "10001d20000_0","10001d20000_>0",
                                                "5001d10000_0","5001d10000_>0",
                                                "1001d5000_0","1001d5000_>0",
                                                "1d1000_0","1d1000_>0",
                                                "0d_0","0d_>0"))

select_teclass <- c("Unknown","LINE/Tad1","DNA/TcMar-Sagan","LTR/Gypsy","LINE/L1",
                    "DNA/TcMar-Fot1","LTR/Copia","Satellite") 

gene2te_bin_counts_plot <- ggplot(gene2te_bin_counts, aes(y = species_name, x = count, fill = bin_cat_expression))+
  geom_bar(stat = "identity", position = "stack", color = "black")+
  scale_fill_manual(values = bincat_expression_colours)+theme_bw()+
  scale_x_continuous(limits = c(0, 15000), expand = c(0,0))+
  theme(axis.text.y = element_text(size = 14, face = "italic"), axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), 
        legend.title = element_text(size = 18), legend.text = element_text(size = 14))+
  labs(x = "Number of Transcripts", y = "Taxon", fill = "Bin Category Expression Status")

gene2te_bin_counts_plot_teclass <- ggplot(gene2te_bin_counts[gene2te_bin_counts$NearestTEDist_class %in% select_teclass,], 
                                          aes(y = species_name, x = count, fill = NearestTEDist_class))+
  geom_bar(stat = "identity", position = "stack", color = "black")+
  scale_fill_manual(values = tesub_colours)+theme_bw()+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.text.y = element_text(size = 14, face = "italic"), axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), 
        legend.title = element_text(size = 18), legend.text = element_text(size = 14))+
  labs(x = "Number of Transcripts", y = "Taxon", fill = "Bin Category Expression Status")+
  facet_wrap(bin_cat~expression_status, ncol = 2)

ggsave("figures/fungi/gene2te_bin_counts_plot.pdf", dpi = "retina", height = 5, width = 12)
ggsave("figures/fungi/gene2te_bin_counts_teclass_plot.pdf", gene2te_bin_counts_plot_teclass, dpi = "retina", height = 20, width = 12)

colnames(gene2te_functionquant_v2)[21] <- "sample_id"

# TE PFAM domains
te_pfam <- readLines("te_pfam.txt")
gene2te_functionquant_v3 <- gene2te_functionquant_v2
gene2te_functionquant_v3$te_pfam <- "No"
gene2te_functionquant_v3[gene2te_functionquant_v3$signature_accession %in% te_pfam,]$te_pfam <- "Yes"
gene2te_functionquant_v3$old_transcriptid <- substr(gene2te_functionquant_v3$oId, 1, nchar(gene2te_functionquant_v3$oId)-2)

tecopt_test <- gene2te_functionquant_v3 %>% group_by(old_transcriptid, te_pfam) %>% summarize(count = n()) %>% pivot_wider(names_from = te_pfam, values_from = count, values_fill = 0)
tecopt_test$copt_class <- "non_te"
tecopt_test[tecopt_test$Yes > 0 & tecopt_test$No == 0,]$copt_class <- "te_only"
tecopt_test[tecopt_test$Yes > 0 & tecopt_test$No >0,]$copt_class <- "fusion_candidate"
tecopt_test_catonly <- tecopt_test[,c("old_transcriptid","copt_class")]

gene2te_functionquant_v3 <- merge(gene2te_functionquant_v3, tecopt_test_catonly, by = "old_transcriptid")
gene2te_functionquant_v3$bin_cat <- factor(gene2te_functionquant_v3$bin_cat, 
                                           levels = c("0d", "1d1000",
                                                      "1001d5000","5001d10000",
                                                      "10001d20000","20000"))
gene2te_functionquant_v3$copt_class <- factor(gene2te_functionquant_v3$copt_class, 
                                              levels = c("non_te", "fusion_candidate", "te_only"))

# species distribution of te-copt class
species_tecopt_plot <- ggplot(gene2te_functionquant_v3, aes(y = species_name, fill = copt_class))+
  geom_bar(stat = "count", color = "black")+
  theme_bw()+scale_x_continuous(expand = c(0,0), limits = c(0,15000))+
  scale_fill_manual(values = c("#6fa8dc","#674ea7","#e06666"))+
  labs(x = "Count", y = "Species Name", fill = "Gene Class")+
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 14, face = "italic"),
        legend.text =  element_text(size = 14), legend.title = element_text(size = 18))

ggsave("figures/fungi/species_tecopt.pdf", species_tecopt_plot, 
       dpi = "retina", height = 8, width = 12)

gene2te_functionquant_v4 <- pivot_longer(gene2te_functionquant_v3, cols = c("Rep1", "Rep2"), names_to = "replicate", values_to = "rep_tpm")
tefusion_expression_reps_plot <- ggplot(gene2te_functionquant_v4, aes(y = log2(rep_tpm+1), x = replicate, fill = copt_class))+
  geom_boxplot()+theme_bw()+facet_wrap(~species_name)+
  scale_fill_manual(values =c("#6fa8dc","#674ea7","#e06666"))
ggsave("figures/fungi/tefusion_expression.pdf",
       tefusion_expression_reps_plot, dpi = "retina", height = 7, width = 10)

# bin categories of te-copt class
tefusion_only_distanceplot <- ggplot(gene2te_functionquant_v3[gene2te_functionquant_v3$copt_class == "fusion_candidate",], aes(y = species_name, fill = bin_cat))+geom_bar(stat = "count", color = "black", position = "dodge")+
  theme_bw()+scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(values = bincat_colours)+
  labs(x = "Count", y = "Species Name", fill = "TE Distance Class")+
  theme(axis.text.x = element_text(size = 14), axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 14, face = "italic"),
        legend.text =  element_text(size = 14), legend.title = element_text(size = 18))

ggsave("figures/fungi/tefusion_distanceplot.pdf", tefusion_only_distanceplot, 
       dpi = "retina", height = 8, width = 9)
tefusion_only <- gene2te_functionquant_v3[gene2te_functionquant_v3$copt_class == "fusion_candidate" & gene2te_functionquant_v3$signature_accession %notin% te_pfam,]

# te copt class GO enrichment
te_copt_classes <- unique(gene2te_functionquant_v3$copt_class)
sample_id <- unique(gene2te_functionquant_v3$sample_id.x)

sample_coptclass_genelist <- vector("list", length = length(sample_id))
names(sample_coptclass_genelist) <- sample_id

for (i in 1:length(sample_id)){
  sample_coptclass_genelist[[i]] <- vector("list", length = length(te_copt_classes))
  names(sample_coptclass_genelist[[i]]) <- te_copt_classes
  for (j in 1:length(te_copt_classes)){
    sample_coptclass_genelist[[i]][[j]] <- gene2te_functionquant_v3[gene2te_functionquant_v3$copt_class == te_copt_classes[j] & gene2te_functionquant_v3$sample_id.x == sample_id[i],]$transcript_id
  }
}

sample_coptclass_genelist_v2 <- sample_coptclass_genelist[names(sampleid_genelist)]

sampleid_topgo_tecopt_list_bp <- vector("list", length = length(sample_id))
names(sampleid_topgo_tecopt_list_bp) <- sample_id

for (i in 1:length(sampleid_topgo_tecopt_list_bp)){
  sampleid_topgo_tecopt_list_bp[[i]] <- vector("list", length = length(te_copt_classes))
  names(sampleid_topgo_tecopt_list_bp[[i]]) <- te_copt_classes
  for (j in 1:length(te_copt_classes)){
    markers <- sample_coptclass_genelist_v2[[i]][[j]]
    genelist <-  sampleid_genelist[[i]]
    sample_gene2go <- sampleid_Gene2GO_bp[[i]]
    sampleid_topgo_tecopt_list_bp[[i]][[j]] <- calc_topgo_fisher(markers, genelist, sample_gene2go, "BP", pval = 0.05)
  }
}

topgo_coopt_bp_df_list_all <- vector("list", length = length(sampleid_topgo_tecopt_list_bp))
names(topgo_coopt_bp_df_list_all) <- sample_id

topgo_coopt_bp_df_list_tmp <- vector("list", length = length(te_copt_classes))
for (i in 1:length(sampleid_topgo_tecopt_list_bp)){
  for (j in 1:length(te_copt_classes)){
    result <- sampleid_topgo_tecopt_list_bp[[i]][[j]]
    if (nrow(result) >= 1){
      bincat <- names(sampleid_topgo_tecopt_list_bp[[i]])[j]
      sample <- names(sampleid_topgo_tecopt_list_bp)[i]
      GO_term <- result$GO.ID
      description <- result$Term
      FisherPval <- as.numeric(result$FisherPval)
      logPval <- log10(FisherPval)
      Genes <- result$Genes
      topgo_df <- data.frame("sample_id" = sample, "bin_category" = bincat,
                             "GO.ID" = GO_term,"description" = description,
                             "FisherPval" = FisherPval,"logPval" = logPval, 
                             "Genes" = Genes)
      topgo_coopt_bp_df_list_tmp[[j]] <- topgo_df 
    }
  }
  topgo_coopt_bp_df_list_all[[i]] <- do.call(rbind, topgo_coopt_bp_df_list_tmp)
}

topgp_copt_bp_df_all <- do.call(rbind, topgo_coopt_bp_df_list_all)

#gene2te density plots
gene2te_dist_df <- rbind(gene2te_dist_otherfungi[,c("NearestTEDist","sample_id")], gene2te_functionquant_v2[,c("NearestTEDist","sample_id")])
gene2te_dist_df$sample_id2 <- paste0(gene2te_dist_df$sample_id, ".fa")
gene2te_0dist_density_plot <- ggplot(gene2te_dist_df, aes(x = NearestTEDist+1, y = sample_id2, group = sample_id, fill = after_stat(x)))+
  geom_density_ridges_gradient()+scale_x_log10()+
  scale_fill_viridis_c(option = "inferno", trans = "log",direction = -1)+theme_bw()+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 18))+
  labs(x = "log(Nearest TE Distance)+1", fill = "log(Nearest TE Distance)+1")
gene2te_0dist_density_tree_plot <- gene2te_0dist_density_plot %>% insert_left(fungi_mag_ggtree)
ggsave("figures/fungi/gene2te_0dist_densityplot.pdf", gene2te_0dist_density_tree_plot, 
      dpi = "retina", height = 10, width = 15)

gene2te_density_plot <- ggplot(gene2te_dist_df[gene2te_dist_df$NearestTEDist > 0,], aes(x = NearestTEDist, y = sample_id2, group = sample_id, fill = after_stat(x)))+
  geom_density_ridges_gradient()+scale_x_log10()+
  scale_fill_viridis_c(option = "inferno", trans = "log",direction = -1)+theme_bw()+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 18), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 18))+
  labs(x = "log(Nearest TE Distance)", fill = "log(Nearest TE Distance)")
gene2te_density_tree_plot <- gene2te_density_plot %>% insert_left(fungi_mag_ggtree)
ggsave("figures/fungi/gene2te_densityplot.pdf", gene2te_density_tree_plot, 
       dpi = "retina", height = 10, width = 15)

# GO term enrichemnt by TE distance bins
te_bin_categories <- unique(gene2te_functionquant_v2$bin_cat)
sample_id <- unique(gene2te_functionquant_v2$sample_id.x)

sampleid_goterm <- vector("list", length = length(sample_id))
names(sampleid_goterm) <- sample_id

for (i in 1:length(sample_id)){
  transcripts <- unique(gene2te_functionquant_v2[gene2te_functionquant_v2$sample_id.x == sample_id[i],]$transcript_id)
  sampleid_goterm[[i]] <- vector("list", length = length(transcripts))
  names(sampleid_goterm[[i]]) <- transcripts
  for (j in 1:length(sampleid_goterm[[i]])){
    transcript_id <- transcripts[j]
    go_terms <- gene2te_functionquant_v2[gene2te_functionquant_v2$transcript_id == transcript_id,]$go_term
    go_terms <- unlist(strsplit(go_terms, "|", fixed = TRUE))
    sampleid_goterm[[i]][[j]] <- go_terms[!grepl("-",go_terms)] 
  }
}

sampleid_gene2go_cc <- vector("list", length = length(sample_id))
names(sampleid_gene2go_cc) <- sample_id
sampleid_Gene2GO_cc <- vector("list", length = length(sample_id))
names(sampleid_Gene2GO_cc) <- sample_id

sampleid_gene2go_bp <- vector("list", length = length(sample_id))
names(sampleid_gene2go_bp) <- sample_id
sampleid_Gene2GO_bp <- vector("list", length = length(sample_id))
names(sampleid_Gene2GO_bp) <- sample_id

sampleid_gene2go_mf <- vector("list", length = length(sample_id))
names(sampleid_gene2go_mf) <- sample_id
sampleid_Gene2GO_mf <- vector("list", length = length(sample_id))
names(sampleid_Gene2GO_mf) <- sample_id

for (i in 1:length(sample_id)){
  sample <- sample_id[i]
  goterms <- sampleid_goterm[[i]]
  sampleid_gene2go_bp[[i]] <- annFUN.gene2GO(whichOnto = "BP", gene2GO = goterms)
  sampleid_Gene2GO_bp[[i]] <- inverseList(sampleid_gene2go_bp[[i]])
}

for (i in 1:length(sample_id)){
  sample <- sample_id[i]
  goterms <- sampleid_goterm[[i]]
  sampleid_gene2go_mf[[i]] <- annFUN.gene2GO(whichOnto = "MF", gene2GO = goterms)
  sampleid_Gene2GO_mf[[i]] <- inverseList(sampleid_gene2go_mf[[i]])
}

for (i in 1:length(sample_id)){
  sample <- sample_id[i]
  goterms <- sampleid_goterm[[i]]
  sampleid_gene2go_cc[[i]] <- annFUN.gene2GO(whichOnto = "CC", gene2GO = goterms)
  sampleid_Gene2GO_cc[[i]] <- inverseList(sampleid_gene2go_cc[[i]])
}

sample_bincat_genelist <- vector("list", length = length(sample_id))
names(sample_bincat_genelist) <- sample_id
bincategories <- unique(gene2te_functionquant_v2$bin_cat)

for (i in 1:length(sample_id)){
  sample_bincat_genelist[[i]] <- vector("list", length = length(bincategories))
  names(sample_bincat_genelist[[i]]) <- bincategories
  for (j in 1:length(bincategories)){
    sample_bincat_genelist[[i]][[j]] <- gene2te_functionquant_v2[gene2te_functionquant_v2$bin_cat == bincategories[j] & gene2te_functionquant_v2$sample_id == sample_id[i],]$transcript_id
  }
}

sampleid_genelist <- vector("list", length = length(sample_id))
names(sampleid_genelist) <- sample_id

for (i in 1:length(sample_id)){
  sampleid_genelist[[i]] <- gene2te_functionquant_v2[gene2te_functionquant_v2$sample_id == sample_id[[i]],]$transcript_id
}

sampleid_topgo_list_bp <- vector("list", length = length(sample_id))
names(sampleid_topgo_list_bp) <- sample_id

for (i in 1:length(sampleid_topgo_list_bp)){
  sampleid_topgo_list_bp[[i]] <- vector("list", length = length(bincategories))
  names(sampleid_topgo_list_bp[[i]]) <- bincategories
  for (j in 1:length(bincategories)){
    markers <- sample_bincat_genelist[[i]][[j]]
    genelist <-  sampleid_genelist[[i]]
    sample_gene2go <- sampleid_Gene2GO_bp[[i]]
    sampleid_topgo_list_bp[[i]][[j]] <- calc_topgo_fisher(markers, genelist, sample_gene2go, "BP", pval = 0.05)
  }
}

topgo_bp_df_list_all <- vector("list", length = length(sampleid_topgo_list_bp))
names(topgo_bp_df_list_all) <- sample_id

topgo_bp_df_list_tmp <- vector("list", length = length(bincategories))
for (i in 1:length(sampleid_topgo_list_bp)){
  for (j in 1:length(bincategories)){
    result <- sampleid_topgo_list_bp[[i]][[j]]
    if (nrow(result) >= 1){
      bincat <- names(sampleid_topgo_list_bp[[i]])[j]
      sample <- names(sampleid_topgo_list_bp)[i]
      GO_term <- result$GO.ID
      description <- result$Term
      FisherPval <- as.numeric(result$FisherPval)
      logPval <- log10(FisherPval)
      Genes <- result$Genes
      topgo_df <- data.frame("sample_id" = sample, "bin_category" = bincat,
                             "GO.ID" = GO_term,"description" = description,
                             "FisherPval" = FisherPval,"logPval" = logPval, 
                             "Genes" = Genes)
      topgo_bp_df_list_tmp[[j]] <- topgo_df 
    }
  }
  topgo_bp_df_list_all[[i]] <- do.call(rbind, topgo_bp_df_list_tmp)
}

topgo_bp_df_all <- do.call(rbind, topgo_bp_df_list_all)
topgo_bp_df_all$bin_category <- factor(topgo_bp_df_all$bin_category, 
                                       levels = c("20000","10001d20000","5001d10000",
                                                  "1001d5000","1d1000","0d"))
topgo_bp_df_all$species <- unname(species_mag_codes[topgo_bp_df_all$sample_id])
topgo_bp_df_all_summary <- topgo_bp_df_all %>% group_by(sample_id, species,bin_category) %>%
  summarize(count = n())

topgo_bp_enrichedterms_plot <- ggplot(topgo_bp_df_all_summary, aes(y = species, x = count, fill = bin_category))+
  geom_bar(stat = "identity", color = "black")+theme_bw()+scale_fill_manual(values = bincat_colours)+
  scale_x_continuous(expand = c(0,0), limits = c(0,600))+
  theme(axis.text.y = element_text(size = 14, face = "italic"), 
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14))+
  labs(x = "Enriched GO Term Count", y = "Taxon", fill = "TE Distance Category")

ggsave("figures/fungi/topgo_bp_enrichmentcount_plot.pdf",
       topgo_bp_enrichedterms_plot, dpi = "retina", height = 5, width = 10)

topgo_bp_df_0d <- topgo_bp_df_all[topgo_bp_df_all$bin_category == "0d",]
topgo_bp_df_0d <- topgo_bp_df_0d[order(topgo_bp_df_0d$logPval),]
topgo_bp_df_1d1000 <- topgo_bp_df_all[topgo_bp_df_all$bin_category == "1d1000",]
topgo_bp_df_1d1000 <- topgo_bp_df_1d1000[order(topgo_bp_df_1d1000$logPval),]

topgo_bp_df_20000 <- topgo_bp_df_all[topgo_bp_df_all$bin_category == "20000",]
topgo_bp_df_20000 <- topgo_bp_df_20000[order(topgo_bp_df_20000$logPval),]
ggplot(topgo_bp_df_0d[topgo_bp_df_0d$`-logPval` > 2.5,], aes(x = -logPval, y = description, fill = species))+
  geom_bar(stat = "identity")+scale_fill_viridis_d(option = "turbo")+
  facet_wrap(~species, nrow = 3, ncol = 3, scales = "free")

# bubble plot sized by # of species for significant GO terms
go_responseterms <- topgo_bp_df_0d[grepl("response",topgo_bp_df_0d$description),]$description
go_heterocycleterms <- unique(topgo_bp_df_all[grepl("heterocycle", topgo_bp_df_all$description),]$description)
go_toxinterms <- unique(topgo_bp_df_all[grepl("toxin", topgo_bp_df_all$description),]$description)
go_secondarymetabterms <- c(go_toxinterms[!(grepl("aflatoxin", go_toxinterms))], go_heterocycleterms)
go_rnaprocess <- unique(topgo_bp_df_all[grepl("RNA processing", topgo_bp_df_all$description),]$description)
go_cellcycle <- unique(topgo_bp_df_all[grepl("cell cycle", topgo_bp_df_all$description),]$description)

topgo_bp_distance_summary <- topgo_bp_df_all %>% group_by(description, bin_category) %>% summarize(count = n())
topgo_bp_distance_summary$bin_category <- factor(topgo_bp_distance_summary$bin_category, 
                                                 levels =  rev(c("20000",
                                                                     "10001d20000",
                                                                     "5001d10000",
                                                                     "1001d5000",
                                                                     "1d1000",
                                                                     "0d")))
topgo_bp_df_all_v2 <- topgo_bp_df_all
topgo_bp_df_all_v2$genecat <- "n/a"
topgo_bp_df_all_v2[topgo_bp_df_all_v2$description %in% c(go_cellcycle, go_rnaprocess),]$genecat <- "cell_cycle"
topgo_bp_df_all_v2[topgo_bp_df_all_v2$description %in% c(go_secondarymetabterms, go_responseterms),]$genecat <- "stress"
topgo_bp_df_all_v2$bin_category <- factor(topgo_bp_df_all_v2$bin_category, 
                                          levels =  rev(c("20000",
                                                          "10001d20000",
                                                          "5001d10000",
                                                          "1001d5000",
                                                          "1d1000",
                                                          "0d")))

goresponse_bin_plot <- ggplot(topgo_bp_distance_summary[topgo_bp_distance_summary$description %in% c(go_secondarymetabterms, go_responseterms),], aes(x = bin_category, y = description, color = bin_category, size = count))+
  geom_point(alpha = 0.7)+theme_bw()+scale_color_manual(values = bincat_colours)+
  theme(legend.position = "none", axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank())+
  labs(x = "Bin Distance Category", y = "GO Term Description")

gocellcycle_bin_plot <- ggplot(topgo_bp_distance_summary[topgo_bp_distance_summary$description %in% c(go_cellcycle, go_rnaprocess),], aes(x = bin_category, y = description, color = bin_category, size = count))+
  geom_point(alpha = 0.7)+theme_bw()+scale_color_manual(values = bincat_colours)+
  theme(legend.position = "none", axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank())+
  labs(x = "Bin Distance Category", y = "GO Term Description")

go_bin_plot_combined <- plot_grid(goresponse_bin_plot, gocellcycle_bin_plot, nrow = 2)

topgo_species_stress_plot <- ggplot(topgo_bp_df_all_v2[topgo_bp_df_all_v2$genecat %in% c("stress"),], aes(x = bin_category, y = description, color = bin_category, size = -logPval))+
  geom_point(alpha = 0.7)+theme_bw()+scale_color_manual(values = bincat_colours)+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), legend.title = element_text(size = 18),
        legend.text = element_text(size = 14), 
        strip.text.x = element_text(size = 14))+
  labs(x = "Bin Distance Category", y = "GO Term Description", size = "-log(Pval)")+facet_wrap(~sample_id, nrow = 1)

topgo_species_cell_plot <- ggplot(topgo_bp_df_all_v2[topgo_bp_df_all_v2$genecat %in% c("cell_cycle"),], aes(x = bin_category, y = description, color = bin_category, size = -logPval))+
  geom_point(alpha = 0.7)+theme_bw()+scale_color_manual(values = bincat_colours)+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), legend.title = element_text(size = 18),
        legend.text = element_text(size = 14), 
        strip.text.x = element_text(size = 14) 
        )+
  labs(x = "Bin Distance Category", y = "GO Term Description", size = "-log(Pval)")+facet_wrap(~sample_id, nrow = 1)

topgo_species_plot_combined <- plot_grid(topgo_species_stress_plot, topgo_species_cell_plot, nrow = 2)
ggsave("figures/fungi/enriched_goterm_bincat.pdf",
go_bin_plot_combined, dpi = "retina", height = 10, width = 8)

ggsave("fungi/topgo_species_plot_combined.pdf",
       topgo_species_plot_combined, dpi = "retina", height = 16, width = 12)
ggsave("figures/fungi/topgo_species_plot_cell.pdf",
       topgo_species_cell_plot, dpi = "retina", height = 6, width = 20)
ggsave("figures/fungi/topgo_species_plot_stress.pdf",
       topgo_species_stress_plot, dpi = "retina", height = 6, width = 20)
gene2te_functionquant_v3 <- gene2te_functionquant_v2
gene2te_functionquant_v3[is.na(gene2te_functionquant_v3$RNA),]$RNA <- 0
gene2te_functionquant_v3$logTPM <- log(gene2te_functionquant_v3$RNA + 1)
gene2te_functionquant_v3$species_name2 <- gsub("_", " ",gene2te_functionquant_v3$species_name)

gene2te_logtpm_plot <- ggplot(gene2te_functionquant_v3, aes(x = NearestTEDist, y = logTPM))+
  geom_pointdensity()+theme_bw()+scale_color_viridis_c(option = "inferno")+
  facet_wrap(~species_name2, scales = "free", nrow = 3, ncol = 3)+
  labs(x = "Distance to Nearest TE (bp)", y = "log(TPM + 1)", color = "Density")+
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18),
        strip.text = element_text(size = 18, face = "italic"))
  
ggsave("figures/fungi/gene2te_logtpm_plot.pdf", 
       gene2te_logtpm_plot, dpi = "retina", height = 10, width = 14)

# Kimura distances
kimuradir <- "earlgrey/fungi/kimura_distances/"
divsumfiles <- list.files(kimuradir)
divsum_accessions <- substr(divsumfiles, 1,nchar(divsumfiles)-7)
fungi_metadata <- read.csv("fungi_metadata.csv")
fungi_metadata$accession <- substr(fungi_metadata$accession, 1, nchar(fungi_metadata$accession)-3)
kimuradf_list <- vector("list", length(divsumfiles))
names(kimuradf_list) <- divsum_accessions

for (i in 1:length(divsum_accessions)){
  kimuradf_list[[i]] <- kimura_df(paste0(kimuradir,divsumfiles[i]),divsum_accessions[i], fungi_metadata)
}

kimuradf_list_mags <- kimuradf_list[rev(substr(fungi_mag_tree_sampleid, 1, nchar(fungi_mag_tree_sampleid)-3))]
fungimagplot_order <- rev(substr(fungi_mag_tree_sampleid, 1, nchar(fungi_mag_tree_sampleid)-3))

kimura_plot_list <- vector("list", length(fungimagplot_order))
names(kimura_plot_list) <- fungimagplot_order

for (i in 1:length(kimuradf_list_mags)){
  plot_df <- kimuradf_list_mags[[i]]
  limits <- plot_df %>% group_by(Divergence) %>% summarise(perc = sum(proportion * 100))
  ylimit <- max(limits$perc) * 1.75
  kimura_plot_list[[i]] <- ggplot(plot_df[! plot_df$classif %like% "Other",], aes(x = Divergence, y = proportion*100, fill = classif)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylimit))+theme_void()+scale_fill_manual(values = te_colours)+
    theme(legend.position = "none")
}

allkimplots <- ggarrange(plotlist = kimura_plot_list, ncol = 1)

ggsave("figures/fungi/voidkimura.pdf", height = 10, width = 7, dpi = "retina")

# CAZyme results
#gene2te_functionquant_v2[is.na(gene2te_functionquant_v2$RNA),]$RNA <- 0
gene2te_functionquant_v2$logRNA <- log(as.numeric(gene2te_functionquant_v2$RNA)+1)
dbcan_dir <- "dbcan_results/"
dbcan_files <- list.files(dbcan_dir, full.names = TRUE)

dbcan_list <- vector("list", length = length(dbcan_files))

for (i in 1:length(dbcan_files)){
  dbcan_list[[i]] <- read.csv(dbcan_files[[i]], sep = "\t")
}

dbcan_results <- do.call(rbind, dbcan_list)
dbcan_results$transcript_id <- substr(dbcan_results$Gene.ID, 1, (nchar(dbcan_results$Gene.ID)-3))
dbcan_results_3hits <- dbcan_results[dbcan_results$X.ofTools == 3,]
dbcan_results_3hits$sample <- word(dbcan_results_3hits$transcript_id, sep = "_")
dbcan_df <- data.frame("transcript_id" = dbcan_results_3hits$transcript_id, "dbcan_hit" = dbcan_results_3hits$dbCAN_sub)

dbcan_results_3hits$cazy_class <- substring(dbcan_results_3hits$dbCAN_sub, 1, 2)
for (i in 1:nrow(dbcan_results_3hits)){
  if (grepl("+", dbcan_results_3hits$dbCAN_sub[i], fixed = TRUE)){
    tmp <- strsplit(dbcan_results_3hits$dbCAN_sub[i], "+", fixed = TRUE)
    tmp_class <- substring(tmp[[1]], 1, 2) %>% unique()
    dbcan_results_3hits$cazy_class[i] <- paste(tmp_class, collapse = ",")
  }
}

#fix species names
dbcan_results_3hits[dbcan_results_3hits$sample == "PseNorv",]$sample <- "glPseNorv1_bin.13"
dbcan_results_3hits[dbcan_results_3hits$sample == "PelMemb",]$sample <- "glPelMemb1_bin.35"
dbcan_results_3hits[dbcan_results_3hits$sample == "PelHori",]$sample <- "glPelHori1_bin.14"
dbcan_results_3hits[dbcan_results_3hits$sample == "NepLaev",]$sample <- "glNepLaev11_merged.0"
dbcan_results_3hits[dbcan_results_3hits$sample == "glRicVire11",]$sample <- "glRicVire11_merged.0"
dbcan_results_3hits[dbcan_results_3hits$sample == "glPelPrae3",]$sample <- "glPelPrae3_bin.50"
dbcan_results_3hits[dbcan_results_3hits$sample == "glPelHyme1",]$sample <- "glPelHyme1_bin.206"
dbcan_results_3hits[dbcan_results_3hits$sample == "glLobPulm2",]$sample <- "glLobPulm2_bin.43"
dbcan_results_3hits[dbcan_results_3hits$sample == "glLepBurg3",]$sample <- "glLepBurg3_bin.203"
lichenmagplotorder <- fungimagplot_order[1:11] %>% rev()
dbcan_results_3hits$sample <- factor(dbcan_results_3hits$sample, levels = lichenmagplotorder)

dbcan_df_summary <- dbcan_results_3hits %>% group_by(sample, cazy_class) %>% summarize(count = n()) 


cazy_count_plot <- ggplot(dbcan_df_summary[dbcan_df_summary$cazy_class %in% c("AA","CB","GH","GT"),], aes(y = sample, x = count, fill = cazy_class))+
  geom_bar(stat = "identity", color = "black")+
  theme_bw()+scale_x_continuous(expand = c(0,0), limits = c(0,500))+
  scale_fill_manual(values = cazy12_palette)+
  labs(fill = "CAZyme Class", y = "Taxon", x = "Count")+
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14, face =  "italic"),
        axis.title.y = element_text(size = 18), axis.title.x = element_text(size = 18),
        legend.title = element_text(size = 18), legend.text = element_text(size = 14))
ggsave("figures/fungi/cazy_count_plot.pdf", 
       dpi = "retina", height =5, width = 8)

dbcan_df$transcript_id <- gsub("NepLaev", "glNepLaev11",dbcan_df$transcript_id)
gene2te_functionquant_dbcan <- merge(gene2te_functionquant_v2, dbcan_df, all.x = TRUE)
gene2te_dbcanonly <- merge(gene2te_functionquant_v2, dbcan_df)
gene2te_dbcanonly$logRNA <- log10(gene2te_dbcanonly$RNA+1)
gene2te_dbcanonly$species_name2 <- gsub("_", " ", gene2te_dbcanonly$species_name)
gene2te_dbcanonly
gene2te_dbcanonly$bin_cat <- factor(gene2te_dbcanonly$bin_cat, levels = c("20000","10001d20000","5001d10000",
                                                                          "1001d5000","1d1000","0d"))
cazyme2te_densitypoint_plot <- ggplot(gene2te_dbcanonly, aes(x = NearestTEDist, y = logRNA))+geom_pointdensity()+
  facet_wrap(~species_name,nrow = 3, ncol = 3)+theme_bw()+
  scale_color_viridis_c(option = "inferno")+
  facet_wrap(~species_name2, scales = "free", nrow = 3, ncol = 3)+
  labs(x = "Distance to Nearest TE (bp)", y = "log(TPM + 1)", color = "Density")+
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18),
        strip.text = element_text(size = 18, face = "italic"))

ggsave("figures/fungi/cazyme2te_density_plot.pdf", 
       cazyme2te_densitypoint_plot, dpi = "retina", height = 10, width = 14)

ggplot(gene2te_dbcanonly[gene2te_dbcanonly$logRNA > 2.5,], aes(y = dbcan_hit, x = logRNA, fill = bin_cat))+
  geom_bar(stat = "identity")+facet_wrap(~species_name, scales = "free", nrow = 3, ncol = 3)+
  scale_fill_viridis_d(option = "inferno")

gene2te_dbcanonly$dbcan_class <- substr(gene2te_dbcanonly$dbcan_hit,1,2)
ggplot()+geom_density(data = gene2te_functionquant_v2, aes(x = NearestTEDist+1))+
  geom_density(data = gene2te_dbcanonly[gene2te_dbcanonly$dbcan_class %in% c("GH","GT","AA"),], aes(x = NearestTEDist+1, color = dbcan_class, group = dbcan_class))+
  scale_x_log10()+facet_wrap(~species_name, scales = "free_y")+theme_bw()+
  scale_color_manual(values = cazy12_palette[c(4,8,10)])+scale_y_continuous(limits = c(0,1))

# multicellularity
multicellularity_genes <- read.csv("multicellularity.txt", sep = "\t")
gene2te_multicellularity <- gene2te_functionquant_dbcan[gene2te_functionquant_dbcan$ipro_accesion %in% multicellularity_genes$ID | gene2te_functionquant_dbcan$dbcan_hit %in% multicellularity_genes$ID,]
gene2te_multicellular_summary <- gene2te_multicellularity %>% group_by(species_name, expression_status, bin_cat) %>% summarize(count = n())

gene2te_multicellular_summary$expression_status2 <- "Not expressed"
gene2te_multicellular_summary[gene2te_multicellular_summary$expression_status == ">0",]$expression_status2 <- "Expressed"

gene2te_multicellular_summary$bin_cat2 <- "0bp (Internal)"
gene2te_multicellular_summary[gene2te_multicellular_summary$bin_cat == "10001d20000",]$bin_cat2 <- "10,001-20,000bp"
gene2te_multicellular_summary[gene2te_multicellular_summary$bin_cat == "1001d5000",]$bin_cat2 <- "1001-5000bp"
gene2te_multicellular_summary[gene2te_multicellular_summary$bin_cat == "1d1000",]$bin_cat2 <- "1-1000bp"
gene2te_multicellular_summary[gene2te_multicellular_summary$bin_cat == "20000",]$bin_cat2 <- ">20,000bp"
gene2te_multicellular_summary[gene2te_multicellular_summary$bin_cat == "5001d10000",]$bin_cat2 <- "5001-10,000bp"

gene2te_multicellular_summary$binexp2 <- paste0(gene2te_multicellular_summary$bin_cat, " - ", gene2te_multicellular_summary$expression_status)
gene2te_multicellular_summary$binexp <- paste0(gene2te_multicellular_summary$bin_cat,"_", gene2te_multicellular_summary$expression_status)

gene2te_multicellular_summary$binexp <- factor(gene2te_multicellular_summary$binexp, 
                                                levels = c("20000_0","20000_>0",
                                                           "10001d20000_0","10001d20000_>0",
                                                           "5001d10000_0","5001d10000_>0",
                                                           "1001d5000_0","1001d5000_>0",
                                                           "1d1000_0","1d1000_>0",
                                                           "0d_0","0d_>0"))

ggplot(gene2te_multicellular_summary, aes(x = count, y = species_name, fill = binexp))+
  geom_bar(stat = "identity", position = "stack", color = "black")+theme_bw()+scale_fill_manual(values = bincat_expression_colours)+
  scale_x_continuous(limits = c(0,1000), expand = c(0,1))
  
ggplot(gene2te_multicellularity, aes(x = NearestTEDist, y = logRNA))+geom_pointdensity()+
  theme_bw()+
  scale_color_viridis_c(option = "inferno")+
  facet_wrap(~species_name, scales = "free", nrow = 3, ncol = 3)+
  labs(x = "Distance to Nearest TE (bp)", y = "log(TPM + 1)", color = "Density")+
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18),
        strip.text = element_text(size = 18, face = "italic"))

# contig plots
fungi_fasta_files <- list.files("mags/fungi/")
fungi_fasta_files <-fungi_fasta_files[-7]

fungi_fasta_file_stringset_list <- vector("list", length = length(fungi_fasta_files))
names(fungi_fasta_file_stringset_list) <- fungi_fasta_files

fungi_bed_list <- vector("list", length = length(fungi_fasta_files))
names(fungi_bed_list) <- fungi_fasta_files

fungi_telomere_files <- list.files("telomeres/fungi/")

fungi_gcdf_list <- vector("list", length = length(fungi_fasta_files))
names(fungi_gcdf_list) <- fungi_fasta_files

fungi_contig_plot_df_list <- vector("list", length = length(fungi_fasta_files))
names(fungi_contig_plot_df_list) <- fungi_fasta_files

for (i in 1:length(fungi_fasta_files)){
  fungipath <- paste0("mags/fungi/", fungi_fasta_files[i])
  fungi_fasta_file_stringset_list[[i]] <- readDNAStringSet(fungipath)
  fungi_bed_list[[i]] <- data.frame("seq" = names(fungi_fasta_file_stringset_list[[i]]), "start" = 1, end = width(fungi_fasta_file_stringset_list[[i]])) %>% arrange(desc(end))
  fungi_contig_count <- length(fungi_fasta_file_stringset_list[[i]])
  fungifasta_newnames <- structure(paste0("contig_", seq(from = 1, to = fungi_contig_count)), names = fungi_bed_list[[i]]$seq)
  fungi_bed_list[[i]]$seq <- fungifasta_newnames[fungi_bed_list[[i]]$seq]
  
  if (sum(width(fungi_fasta_file_stringset_list[[i]])) > 35000000){
    fungi_gdf <- as.data.frame(get_bkg(fungi_fasta_file_stringset_list[[i]], k = 1, RC = F, window = TRUE, window.size = 10000, merge.res = FALSE))
    fungi_c_df <- fungi_gdf[fungi_gdf$klet == "C",]
    fungi_gg_df <- fungi_gdf[fungi_gdf$klet == "G",]
    fungi_gc_df <- data.frame(seq = fungi_c_df$sequence, start = fungi_c_df$start, end = fungi_c_df$stop, score = (fungi_c_df$probability + fungi_gg_df$probability)/2)
    fungi_gc_df$seq <- fungifasta_newnames[fungi_gc_df$seq]
  }else{
    fungi_gdf <- as.data.frame(get_bkg(fungi_fasta_file_stringset_list[[i]], k = 1, RC = F, window = TRUE, window.size = 1000, merge.res = FALSE))
    fungi_c_df <- fungi_gdf[fungi_gdf$klet == "C",]
    fungi_gg_df <- fungi_gdf[fungi_gdf$klet == "G",]
    fungi_gc_df <- data.frame(seq = fungi_c_df$sequence, start = fungi_c_df$start, end = fungi_c_df$stop, score = (fungi_c_df$probability + fungi_gg_df$probability)/2)
    fungi_gc_df$seq <- fungifasta_newnames[fungi_gc_df$seq]
  }

  gff_colnames <- c("seq","source","feature","start","end","score","strand","frame","attribute")
  
  fungitelomerepath <- paste0("telomeres/fungi/", fungi_fasta_files[i], "_telomeres.txt")
  if (nrow(read.delim(fungitelomerepath)) > 0 ){
    fungi_telomeres <- read.delim(fungitelomerepath, header = F)[,c(1,2)]
    fungi_telomeres$V1 <- fungifasta_newnames[fungi_telomeres$V1]
    fungi_telomeres <- fungi_telomeres %>% na.omit()
    colnames(fungi_telomeres) <- c("seq","direction")
    fungi_telomeres$position <- 0 
    fungi_telomeres$start <- fungi_bed_list[[i]][fungi_bed_list[[i]]$seq  %in% fungi_telomeres$seq,]$start
    fungi_telomeres$end <- fungi_bed_list[[i]][fungi_bed_list[[i]]$seq  %in% fungi_telomeres$seq,]$end
    fungi_telomeres[fungi_telomeres$direction == "forward",]$position <- fungi_telomeres[fungi_telomeres$direction == "forward",]$start
    fungi_telomeres[fungi_telomeres$direction == "reverse",]$position <- fungi_telomeres[fungi_telomeres$direction == "reverse",]$end
    fungi_gc_df$telomere <- "no"
    fungi_gc_df[fungi_gc_df$seq %in% fungi_telomeres$seq,]$telomere <- "yes"
    fungi_contig_plot_df <- left_join(fungi_gc_df, fungi_telomeres[,c(1,3)])  
  }else
    {
      fungi_contig_plot_df <- fungi_gc_df
      fungi_contig_plot_df$telomere <- "no"
      fungi_contig_plot_df$position <- "no"
    }
  #fungi_gc_df$telomere <- "no"
  #fungi_gc_df[fungi_gc_df$seq %in% fungi_telomeres$seq,]$telomere <- "yes"
  #fungi_contig_plot_df <-left_join(fungi_gc_df, fungi_telomeres[,c(1,5)])
  
  fungi_gcdf_list[[i]] <- fungi_gc_df
  fungi_contig_plot_df_list[[i]] <- fungi_contig_plot_df
}

for (i in 1:length(fungi_contig_plot_df_list)){
  print(names(fungi_contig_plot_df_list)[i])
  filename <- paste0("figures/fungi/", 
                     names(fungi_contig_plot_df_list)[i], "_telomereplot.pdf")
  if (length(unique(fungi_contig_plot_df_list[[i]]$telomere)) == 1){
    contigplot_fungi <- ggplot(fungi_contig_plot_df_list[[i]], aes(x = start, y = reorder(seq, end), fill = score*100))+geom_tile()+
      theme_bw()+scale_fill_viridis_c()+
      labs(x = "Contig Length (Mb)", y = "Contig ID", fill = "GC Content (%)")+
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18),
            legend.text = element_text(size = 14), legend.title = element_text(size = 18))
    
    ggsave(filename,contigplot_fungi,
           dpi = "retina", height = 8, width = 9)
  }else{
    contigplot_fungi <- ggplot(fungi_contig_plot_df_list[[i]], aes(x = start, y = reorder(seq, end), fill = score*100))+geom_tile()+
      geom_point(aes(x = position, y = seq), color = "black", pch = "*" , size = 14)+
      theme_bw()+scale_fill_viridis_c()+
      labs(x = "Contig Length (Mb)", y = "Contig ID", fill = "GC Content (%)")+
      theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18),
            legend.text = element_text(size = 14), legend.title = element_text(size = 18))
    
    ggsave(filename,contigplot_fungi,
           dpi = "retina", height = 8, width = 9)
  }
}


ggplot(fungi_contig_plot_df_list[[1]], aes(x = start, y = reorder(seq, end), fill = score*100))+geom_tile()+
  theme_bw()+scale_fill_viridis_c()+
  scale_x_continuous(breaks = c(0,1000000, 2000000, 3000000, 4000000, 5000000,6000000, 7000000,8000000),
                     labels = c("0","1", "2","3", "4", "5", "6", "7", "8"))+
  labs(x = "Contig Length (Mb)", y = "Contig ID", fill = "GC Content (%)")+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18))




symbiochloris_contig_plot <- ggplot(symbiochloris_contig_plot_df, aes(x = start, y = reorder(seq, end), fill = score*100))+geom_tile()+
  geom_point(aes(x = position, y = seq), color = "black", pch = "*" , size = 14)+
  theme_bw()+scale_fill_viridis_c()+
  scale_x_continuous(breaks = c(0,1000000, 2000000, 3000000, 4000000, 5000000,6000000, 7000000,8000000),
                     labels = c("0","1", "2","3", "4", "5", "6", "7", "8"))+
  labs(x = "Contig Length (Mb)", y = "Contig ID", fill = "GC Content (%)")+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18))






