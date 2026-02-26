library(ggtree)
library(aplot)
library(ggplot2)
library(tidyverse)
library(universalmotif)

# input files not provided in the GitHub are available upon request from authors.
# chlorophyte TE tree
chlorophyte_busco_tree <- readRDS("00-setup/RDS/chlorophyte_busco_tree.RDS")

chlorophyte_earlgrey_high_summ_table <- txt_from_folder("earlgrey/chlorophyta/summary_highlevelcounts/")
chlorophyte_earlgrey_high_summ_table$accession <- substr(chlorophyte_earlgrey_high_summ_table$sample_id, 1, 15)
chlorophyte_earlgrey_high_summ_table[chlorophyte_earlgrey_high_summ_table$sample_id == "glRicVire11_bin41.fa",]$accession <- "glRicVire11" 

chlorophyte_earlgrey_busco_table <- chlorophyte_earlgrey_high_summ_table[chlorophyte_earlgrey_high_summ_table$accession %in% chlorophyte_busco_accessions,]
chlorophyte_earlgrey_busco_table <- merge(chlorophyte_earlgrey_busco_table, chlorophyte_busco_sampleid_df)
chlorophyte_earlgrey_busco_table_summary <- chlorophyte_earlgrey_busco_table %>% group_by(sample_id) %>% summarize(prop = sum(proportion*100))
chlorophyte_te_comp <- ggplot(chlorophyte_earlgrey_busco_table, aes(x = proportion*100, y = chlorophyte_busco_sampleid, fill = tclassif))+
  geom_bar(stat = "identity", color = "black")+
  geom_vline(xintercept = median(chlorophyte_earlgrey_busco_table_summary$prop), linetype = "dotdash", colour = "firebrick2" )+
  theme_bw()+scale_x_continuous(expand = c(0,0), limits = c(0,30))+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), legend.position = "none")+
  scale_fill_manual(values = te_colours)+
  labs(fill = "TE Classification", x = "Relative Genome Content (%)")

chlorophyte_te_comp_legend <- chlorophyte_te_comp + theme(legend.position = "right")
te_class_legend <- get_legend(chlorophyte_te_comp_legend) %>% plot_grid()
ggsave("figures/trebouxia_analysis/te_class_legend.pdf", te_class_legend, dpi = "retina", height = 4, width = 3)

trebouxiophyceae_metadata <- read.csv("trebouxiophyceae_metadata.csv", header = TRUE)
trebouxiophyceae_metadata_busco <- merge(trebouxiophyceae_metadata, chlorophyte_busco_sampleid_df) %>% distinct()
trebouxiophyceae_genomesize_bar <- ggplot(trebouxiophyceae_metadata_busco, aes(x = assembly_length, y = chlorophyte_busco_sampleid, fill = busco_habitat))+
  geom_bar(stat = "identity", colour = "black")+
  geom_vline(xintercept = median(trebouxiophyceae_metadata$assembly_length), linetype = "dotdash", colour = "firebrick3")+
  theme_bw()+scale_x_continuous(expand = c(0,0), limits = c(0, 100000000))+
  scale_fill_manual(values = chlorophyte_tree_colours)+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), legend.position = "none")+
  labs(x = "Genome Size")

trebouxiophyceae_genomesize_bar_legend <- trebouxiophyceae_genomesize_bar + theme(legend.position = "right")+labs(fill = "Habitat Classification")
trebouxiophyceae_habitat <- get_legend(trebouxiophyceae_genomesize_bar_legend) %>% plot_grid()
ggsave("figures/trebouxia_analysis/trebouxia_habitat_legend.pdf", trebouxiophyceae_habitat, dpi = "retina", height = 4, width = 3)

chlorophyte_tree_tebars <- trebouxiophyceae_genomesize_bar %>% insert_left(chlorophyte_busco_tree, width = 1.75) %>% insert_right(chlorophyte_te_comp)
ggsave("figures/trebouxia_analysis/chlorophyte_tree_tebars.pdf", 
       chlorophyte_tree_tebars, dpi = "retina", height = 8, width = 15)


# chlorophyte te kimura distance
kimuradir <- "earlgrey/chlorophyta/kimura_distances/"
divsumfiles <- list.files(kimuradir)
divsum_accessions <- substr(divsumfiles, 1,15)
divsum_accessions[92] <- "glLobPulm2"
divsum_accessions[93] <- "glRicVire11"

kimuradf_list <- vector("list", length(divsumfiles))
names(kimuradf_list) <- divsum_accessions
for (i in 1:length(divsum_accessions)){
  kimuradf_list[[i]] <- kimura_df(paste0(kimuradir,divsumfiles[i]),divsum_accessions[i], trebouxiophyceae_metadata)
}

kimuradf_list_busco <- kimuradf_list[chlorophyte_busco_accessions]

kimura_plot_list <- vector("list", length(chlorophyte_busco_accessions))
names(kimura_plot_list) <- chlorophyte_busco_accessions

for (i in 1:length(chlorophyte_busco_accessions)){
  plot_df <- kimuradf_list_busco[[chlorophyte_busco_accessions[i]]]
  limits <- plot_df %>% group_by(Divergence) %>% summarise(perc = sum(proportion * 100))
  ylimit <- max(limits$perc) * 1.75
  kimura_plot_list[[i]] <- ggplot(plot_df[! plot_df$classif %like% "Other",], aes(x = Divergence, y = proportion*100, fill = classif)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylimit))+theme_void()+scale_fill_manual(values = te_colours)+
    theme(legend.position = "none")
}

allkimura_void_plots <- ggarrange(plotlist = kimura_plot_list, ncol = 1)

ggsave("figures/trebouxia_analysis/all_kimura.pdf", 
       allkimura_void_plots, dpi = "retina", height = 1, width = 1)

ricvire_kimura_plot <- ggplot(kimuradf_list_busco[["glRicVire11"]][! kimuradf_list_busco[["glRicVire11"]]$classif %like% "Other",], aes(x = Divergence, y = proportion*100, fill = classif)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 7))+
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_text(face = "italic", angle = 0),
        axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 18), axis.title.x = element_text(size = 22),
        legend.position = "none") +
  scale_fill_manual(values = te_colours)+
  labs(y = "Genome Content (%)", x = "Kimura Distance (CpG Adjusted)", fill = "")

ggsave("figures/trebouxia_analysis/ricvire11_kimura.pdf", 
       ricvire_kimura_plot, dpi = "retina", height = 3.75, width = 24)

CA_032358165.1_kimura_plot <- ggplot(kimuradf_list_busco[["GCA_032358165.1"]][! kimuradf_list_busco[["GCA_032358165.1"]]$classif %like% "Other",], aes(x = Divergence, y = proportion*100, fill = classif)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2.5))+
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_text(face = "italic", angle = 0),
        axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 18), axis.title.x = element_text(size = 22),
        legend.position = "none") +
  scale_fill_manual(values = te_colours)+
  labs(y = "Genome Content (%)", x = "Kimura Distance (CpG Adjusted)", fill = "")

ggsave("figures/trebouxia_analysis/GCA_032358165.1_kimura.pdf", 
       GCA_032358165.1_kimura_plot, dpi = "retina", height = 3.75, width = 24)

GCA_002245815.2_kimura_plot <- ggplot(kimuradf_list_busco[["GCA_002245815.2"]][! kimuradf_list_busco[["GCA_002245815.2"]]$classif %like% "Other",], aes(x = Divergence, y = proportion*100, fill = classif)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 5.5))+
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_text(face = "italic", angle = 0),
        axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 18), axis.title.x = element_text(size = 22),
        legend.position = "none") +
  scale_fill_manual(values = te_colours)+
  labs(y = "Genome Content (%)", x = "Kimura Distance (CpG Adjusted)", fill = "")

ggsave("figures/trebouxia_analysis/GCA_002245815.2_kimura.pdf", 
       GCA_002245815.2_kimura_plot, dpi = "retina", height = 3.75, width = 24)

GCA_031763795.1_kimura_plot <- ggplot(kimuradf_list_busco[["GCA_031763795.1"]][! kimuradf_list_busco[["GCA_031763795.1"]]$classif %like% "Other",], aes(x = Divergence, y = proportion*100, fill = classif)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1))+
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_text(face = "italic", angle = 0),
        axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 18), axis.title.x = element_text(size = 22),
        legend.position = "none") +
  scale_fill_manual(values = te_colours)+
  labs(y = "Genome Content (%)", x = "Kimura Distance (CpG Adjusted)", fill = "")

ggsave("figures/trebouxia_analysis/GCA_031763795.1_kimura.pdf", 
       GCA_031763795.1_kimura_plot, dpi = "retina", height = 3.75, width = 24)

GCA_023343905.1_kimura_plot <- ggplot(kimuradf_list_busco[["GCA_023343905.1"]][! kimuradf_list_busco[["GCA_023343905.1"]]$classif %like% "Other",], aes(x = Divergence, y = proportion*100, fill = classif)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 3.5))+
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_text(face = "italic", angle = 0),
        axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 18), axis.title.x = element_text(size = 22),
        legend.position = "none") +
  scale_fill_manual(values = te_colours)+
  labs(y = "Genome Content (%)", x = "Kimura Distance (CpG Adjusted)", fill = "")

ggsave("figures/trebouxia_analysis/GCA_023343905.1_kimura.pdf", 
       GCA_023343905.1_kimura_plot, dpi = "retina", height = 3.75, width = 24)

GCA_003130725.1_kimura_plot <- ggplot(kimuradf_list_busco[["GCA_003130725.1"]][! kimuradf_list_busco[["GCA_003130725.1"]]$classif %like% "Other",], aes(x = Divergence, y = proportion*100, fill = classif)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.5))+
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_text(face = "italic", angle = 0),
        axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 18), axis.title.x = element_text(size = 22),
        legend.position = "none") +
  scale_fill_manual(values = te_colours)+
  labs(y = "Genome Content (%)", x = "Kimura Distance (CpG Adjusted)", fill = "")

ggsave("figures/trebouxia_analysis/GCA_003130725.1_kimura.pdf", 
       GCA_003130725.1_kimura_plot, dpi = "retina", height = 3.75, width = 24)

GCA_963575595.1_kimura_plot <- ggplot(kimuradf_list_busco[["GCA_963575595.1"]][! kimuradf_list_busco[["GCA_963575595.1"]]$classif %like% "Other",], aes(x = Divergence, y = proportion*100, fill = classif)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.5))+
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_text(face = "italic", angle = 0),
        axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 18), axis.title.x = element_text(size = 22),
        legend.position = "none") +
  scale_fill_manual(values = te_colours)+
  labs(y = "Genome Content (%)", x = "Kimura Distance (CpG Adjusted)", fill = "")

GCA_964019345.1_kimura_plot <- ggplot(kimuradf_list_busco[["GCA_964019345.1"]][! kimuradf_list_busco[["GCA_964019345.1"]]$classif %like% "Other",], aes(x = Divergence, y = proportion*100, fill = classif)) +
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.5))+
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_text(face = "italic", angle = 0),
        axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 18), axis.title.x = element_text(size = 22),
        legend.position = "none") +
  scale_fill_manual(values = te_colours)+
  labs(y = "Genome Content (%)", x = "Kimura Distance (CpG Adjusted)", fill = "")

select_symbiont_kimuraplots <- plot_grid(ricvire_kimura_plot, GCA_032358165.1_kimura_plot,
          GCA_964019345.1_kimura_plot, GCA_963575595.1_kimura_plot, 
          GCA_002245815.2_kimura_plot, GCA_031763795.1_kimura_plot,
          GCA_023343905.1_kimura_plot, GCA_003130725.1_kimura_plot, nrow = 4, align = "v")

ggsave("figures/trebouxia_analysis/select_symbiont_kimura.pdf", select_symbiont_kimuraplots, dpi = "retina", height = 10, width = 32)

# Circos plot for Symbiochloris
svg("figures/trebouxia_analysis/circos.svg")
make_circos_plot(fasta_file = "mags/chlorophyta/glRicVire11_bin.41.fa", gff_file = "earlgrey/chlorophyta/metaeuk_gff/glRicVire11_bin.41.gff", gff_te_file = "earlgrey/chlorophyta/te_gff/glRicVire11_bin.41.filteredRepeats.gff")
dev.off()

# load the gene2te dataframes and repeat analysis like we did for fungi
trebo_gene2te_quant <- readRDS("AlgGeneToTEsFunctionsAQuant.RDS")
trebo_gene2te_quant <- trebo_gene2te_quant[!is.na(trebo_gene2te_quant$RNA),]
trebo_function_quant <- readRDS("AlgFunctionsAQuant.RDS")

trebo_gene2te_quant$RNA <- (trebo_gene2te_quant$Rep1 + trebo_gene2te_quant$Rep2) / 2
trebo_gene2te_quant <- trebo_gene2te_quant[order(trebo_gene2te_quant$RNA, decreasing = TRUE),]
trebo_gene2te_quant_v2 <- trebo_gene2te_quant[!duplicated(trebo_gene2te_quant$transcript_id), ]
trebo_gene2te_quant_v2$expression_status <- "0"
trebo_gene2te_quant_v2[!is.na(trebo_gene2te_quant_v2$RNA) & trebo_gene2te_quant_v2$RNA >= 1,]$expression_status <- ">0"

te_bin_cats <-  c (0, 1, 1000, 5000, 10000, 20000, Inf)
trebo_gene2te_quant_v2$bin_cat <- NA
trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$NearestTEDist <= 1,]$bin_cat <- "0d"
trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$NearestTEDist >= 2 & trebo_gene2te_quant_v2$NearestTEDist <= 1000,]$bin_cat <- "1d1000"
trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$NearestTEDist >= 1001 & trebo_gene2te_quant_v2$NearestTEDist <= 5000,]$bin_cat <- "1001d5000"
trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$NearestTEDist >= 5001 & trebo_gene2te_quant_v2$NearestTEDist <= 10000,]$bin_cat <- "5001d10000"
trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$NearestTEDist >= 10001 & trebo_gene2te_quant_v2$NearestTEDist <= 20000,]$bin_cat <- "10001d20000"
trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$NearestTEDist > 20000,]$bin_cat <- "20000"
trebo_gene2te_quant_v2$species_name <-  sub('\\..*', '',trebo_gene2te_quant_v2$oId)
trebo_gene2te_bin_counts <- trebo_gene2te_quant_v2[!is.na(trebo_gene2te_quant_v2$RNA),] %>% group_by(bin_cat, species_name, expression_status, NearestTEDist_class) %>% summarize(count = n())

trebo_gene2te_bin_counts$bin_cat <- factor(trebo_gene2te_bin_counts$bin_cat, 
                                     levels = c("0d", "1d1000",
                                                "1001d5000","5001d10000",
                                                "10001d20000","20000"))
trebo_gene2te_bin_counts$bin_cat_expression <- paste0(trebo_gene2te_bin_counts$bin_cat, "_",trebo_gene2te_bin_counts$expression_status)

trebo_gene2te_bin_counts$bin_cat_expression <- factor(trebo_gene2te_bin_counts$bin_cat_expression, 
                                                levels = c("20000_0","20000_>0",
                                                           "10001d20000_0","10001d20000_>0",
                                                           "5001d10000_0","5001d10000_>0",
                                                           "1001d5000_0","1001d5000_>0",
                                                           "1d1000_0","1d1000_>0",
                                                           "0d_0","0d_>0"))

trebo_gene2te_bin_counts_v2 <- trebo_gene2te_bin_counts %>% group_by(bin_cat_expression, species_name) %>% summarize(count = sum(count)) %>% distinct()
trebo_gene2te_bin_counts_plot <- ggplot(trebo_gene2te_bin_counts_v2, aes(y = species_name, x = count, fill = bin_cat_expression))+
  geom_bar(stat = "identity", position = "stack", color = "black")+
  scale_fill_manual(values = bincat_expression_colours_trebo)+theme_bw()+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.text.y = element_text(size = 14, face = "italic"), axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), 
        legend.title = element_text(size = 18), legend.text = element_text(size = 14))+
  labs(x = "Number of Transcripts", y = "Taxon", fill = "Bin Category Expression Status")
# many more genes close to TEs in the photobiont 

# choose a subset of these to plot later 
trebo_gene2te_bin_counts_plot_teclass <- ggplot(trebo_gene2te_bin_counts, 
                                          aes(y = species_name, x = count, fill = NearestTEDist_class))+
  geom_bar(stat = "identity", position = "stack", color = "black")+
  scale_fill_viridis_d(option = "turbo")+theme_bw()+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.text.y = element_text(size = 14, face = "italic"), axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), 
        legend.title = element_text(size = 18), legend.text = element_text(size = 14))+
  labs(x = "Number of Transcripts", y = "Taxon", fill = "Bin Category Expression Status")+
  facet_wrap(bin_cat~expression_status, ncol = 2)

ggsave("figures/trebouxia_analysis/gene2te_bin_counts_plot_v2.pdf",trebo_gene2te_bin_counts_plot, dpi = "retina", height = 5, width = 12)
ggsave("figures/trebouxia_analysis/gene2te_bin_counts_teclass_plot.pdf", trebo_gene2te_bin_counts_plot_teclass, dpi = "retina", height = 20, width = 12)

colnames(trebo_gene2te_quant_v2)[21] <- "sample_id"

# GO enrichment
trebo_te_bin_categories <- unique(trebo_gene2te_quant_v2$bin_cat)
trebo_sample_id <- unique(trebo_gene2te_quant_v2$species_name)

trebo_sampleid_goterm <- vector("list", length = length(trebo_sample_id))
names(trebo_sampleid_goterm) <- trebo_sample_id

for (i in 1:length(trebo_sample_id)){
  transcripts <- unique(trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$species_name == trebo_sample_id[i],]$transcript_id)
  trebo_sampleid_goterm[[i]] <- vector("list", length = length(transcripts))
  names(trebo_sampleid_goterm[[i]]) <- transcripts
  for (j in 1:length(trebo_sampleid_goterm[[i]])){
    transcript_id <- transcripts[j]
    go_terms <- trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$transcript_id == transcript_id,]$go_term
    go_terms <- unlist(strsplit(go_terms, "|", fixed = TRUE))
    trebo_sampleid_goterm[[i]][[j]] <- go_terms[!grepl("-",go_terms)] 
  }
}

trebo_sampleid_gene2go_cc <- vector("list", length = length(trebo_sample_id))
names(trebo_sampleid_gene2go_cc) <- trebo_sample_id
trebo_sampleid_Gene2GO_cc <- vector("list", length = length(trebo_sample_id))
names(trebo_sampleid_Gene2GO_cc) <- trebo_sample_id

trebo_sampleid_gene2go_bp <- vector("list", length = length(trebo_sample_id))
names(trebo_sampleid_gene2go_bp) <- trebo_sample_id
trebo_sampleid_Gene2GO_bp <- vector("list", length = length(trebo_sample_id))
names(trebo_sampleid_Gene2GO_bp) <- trebo_sample_id

trebo_sampleid_gene2go_mf <- vector("list", length = length(trebo_sample_id))
names(trebo_sampleid_gene2go_mf) <- trebo_sample_id
trebo_sampleid_Gene2GO_mf <- vector("list", length = length(trebo_sample_id))
names(trebo_sampleid_Gene2GO_mf) <- trebo_sample_id


for (i in 1:length(trebo_sample_id)){
  sample <- trebo_sample_id[i]
  goterms <- trebo_sampleid_goterm[[i]]
  trebo_sampleid_gene2go_bp[[i]] <- annFUN.gene2GO(whichOnto = "BP", gene2GO = goterms)
  trebo_sampleid_Gene2GO_bp[[i]] <- inverseList(trebo_sampleid_gene2go_bp[[i]])
}

for (i in 1:length(trebo_sample_id)){
  sample <- trebo_sample_id[i]
  goterms <- trebo_sampleid_goterm[[i]]
  trebo_sampleid_gene2go_mf[[i]] <- annFUN.gene2GO(whichOnto = "MF", gene2GO = goterms)
  trebo_sampleid_Gene2GO_mf[[i]] <- inverseList(trebo_sampleid_gene2go_mf[[i]])
}

for (i in 1:length(trebo_sample_id)){
  sample <- trebo_sample_id[i]
  goterms <- trebo_sampleid_goterm[[i]]
  trebo_sampleid_gene2go_cc[[i]] <- annFUN.gene2GO(whichOnto = "CC", gene2GO = goterms)
  trebo_sampleid_Gene2GO_cc[[i]] <- inverseList(trebo_sampleid_gene2go_cc[[i]])
}


trebo_sample_bincat_genelist <- vector("list", length = length(trebo_sample_id))
names(trebo_sample_bincat_genelist) <- trebo_sample_id
bincategories <- unique(trebo_gene2te_quant_v2$bin_cat)

for (i in 1:length(trebo_sample_id)){
  trebo_sample_bincat_genelist[[i]] <- vector("list", length = length(bincategories))
  names(trebo_sample_bincat_genelist[[i]]) <- bincategories
  for (j in 1:length(bincategories)){
    trebo_sample_bincat_genelist[[i]][[j]] <- trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$bin_cat == bincategories[j] & trebo_gene2te_quant_v2$species_name == trebo_sample_id[i],]$transcript_id
  }
}

trebo_sampleid_genelist <- vector("list", length = length(trebo_sample_id))
names(trebo_sampleid_genelist) <- trebo_sample_id

for (i in 1:length(trebo_sample_id)){
  trebo_sampleid_genelist[[i]] <- trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$species_name == trebo_sample_id[[i]],]$transcript_id
}


trebo_sampleid_topgo_list_bp <- vector("list", length = length(trebo_sample_id))
names(trebo_sampleid_topgo_list_bp) <- trebo_sample_id

for (i in 1:length(trebo_sampleid_topgo_list_bp)){
  trebo_sampleid_topgo_list_bp[[i]] <- vector("list", length = length(bincategories))
  names(trebo_sampleid_topgo_list_bp[[i]]) <- bincategories
  for (j in 1:length(bincategories)){
    markers <- trebo_sample_bincat_genelist[[i]][[j]]
    genelist <-  trebo_sampleid_genelist[[i]]
    trebo_sample_gene2go <- trebo_sampleid_Gene2GO_bp[[i]]
    trebo_sampleid_topgo_list_bp[[i]][[j]] <- calc_topgo_fisher(markers, genelist, trebo_sample_gene2go, "BP", pval = 0.05)
  }
}

trebo_topgo_bp_df_list_all <- vector("list", length = length(trebo_sampleid_topgo_list_bp))
names(trebo_topgo_bp_df_list_all) <- trebo_sample_id



trebo_topgo_bp_df_list_tmp <- vector("list", length = length(bincategories))
for (i in 1:length(trebo_sampleid_topgo_list_bp)){
  for (j in 1:length(bincategories)){
    result <- trebo_sampleid_topgo_list_bp[[i]][[j]]
    if (nrow(result) >= 1){
      bincat <- names(trebo_sampleid_topgo_list_bp[[i]])[j]
      sample <- names(trebo_sampleid_topgo_list_bp)[i]
      GO_term <- result$GO.ID
      description <- result$Term
      FisherPval <- as.numeric(result$FisherPval)
      logPval <- log10(FisherPval)
      Genes <- result$Genes
      topgo_df <- data.frame("sample_id" = sample, "bin_category" = bincat,
                             "GO.ID" = GO_term,"description" = description,
                             "FisherPval" = FisherPval,"logPval" = logPval, 
                             "Genes" = Genes)
      trebo_topgo_bp_df_list_tmp[[j]] <- topgo_df 
    }
  }
  trebo_topgo_bp_df_list_all[[i]] <- do.call(rbind, trebo_topgo_bp_df_list_tmp)
}

trebo_topgo_bp_df_all <- do.call(rbind, trebo_topgo_bp_df_list_all)
trebo_topgo_bp_df_all$bin_category <- factor(trebo_topgo_bp_df_all$bin_category,
                                             levels = (c("0d", "1d1000",
                                                            "1001d5000","5001d10000",
                                                            "10001d20000","20000")))

trebo_topgo_bp_df_summary <- trebo_topgo_bp_df_all %>% group_by(sample_id, bin_category) %>% summarize(count = n())
trebo_topgo_bp_df_summary$bin_category <- factor(trebo_topgo_bp_df_summary$bin_category,
                                                 levels = rev(c("0d", "1d1000",
                                                            "1001d5000","5001d10000",
                                                           "10001d20000","20000")))

trebogosummary_plot <- ggplot(trebo_topgo_bp_df_summary, aes(y = sample_id, x = count, fill = bin_category))+
  geom_bar(stat = "identity", color = "black")+
  scale_x_continuous(expand = c(0,0), limits = c(0, 500))+
  scale_fill_manual(values = bincat_colours_trebo)+theme_bw()+
  labs(x = "Enriched GO Term Count", y = "Lichen Species", fill = "TE Distance Category")+
  theme(axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18), axis.text.y = element_text(size = 14, face = "italic"))

ggsave("figures/trebouxia_analysis/enrichedgo_plot_v2.pdf",trebogosummary_plot, dpi = "retina", height = 5, width = 12)

vitamin_terms <- unique(trebo_topgo_bp_df_all[grepl("vitamin", trebo_topgo_bp_df_all$description, ignore.case = FALSE),]$description)
cobalamin_terms <- unique(trebo_topgo_bp_df_all[grepl("cobalamin", trebo_topgo_bp_df_all$description, ignore.case = FALSE),]$description)
thiamine_terms <- unique(trebo_topgo_bp_df_all[grepl("thiamine", trebo_topgo_bp_df_all$description, ignore.case = FALSE),]$description)
photosyn_terms <- unique(trebo_topgo_bp_df_all[grepl("photosy", trebo_topgo_bp_df_all$description, ignore.case = FALSE),]$description)
response_terms <- unique(trebo_topgo_bp_df_all[grepl("response", trebo_topgo_bp_df_all$description, ignore.case = FALSE),]$description)
communication_terms <- unique(trebo_topgo_bp_df_all[grepl("communi", trebo_topgo_bp_df_all$description, ignore.case = FALSE),]$description)
rnaprocess_terms <- unique(trebo_topgo_bp_df_all[grepl("RNA processing", trebo_topgo_bp_df_all$description, ignore.case = FALSE),]$description)
ncrna_terms <- unique(trebo_topgo_bp_df_all[grepl("ncRNA", trebo_topgo_bp_df_all$description, ignore.case = FALSE),]$description)
cellcycle_terms <- unique(trebo_topgo_bp_df_all[grepl("cell cycle", trebo_topgo_bp_df_all$description, ignore.case = FALSE),]$description)
mitosis_terms <- unique(trebo_topgo_bp_df_all[grepl("mitotic", trebo_topgo_bp_df_all$description, ignore.case = FALSE),]$description)
heterocycle_terms <- unique(trebo_topgo_bp_df_all[grepl("heterocycle", trebo_topgo_bp_df_all$description, ignore.case = FALSE),]$description)

select_terms <- c(vitamin_terms, cobalamin_terms, thiamine_terms, photosyn_terms, 
                  response_terms, rnaprocess_terms, cellcycle_terms,
                  communication_terms, mitosis_terms, ncrna_terms)
vit_selection <- c(vitamin_terms, cobalamin_terms, thiamine_terms)
lichens_selection <- c(communication_terms, ncrna_terms, heterocycle_terms, response_terms)
essential_terms <- c(photosyn_terms, rnaprocess_terms, cellcycle_terms, mitosis_terms)


trebo_topgo_bp_df_all[trebo_topgo_bp_df_all$GO.ID == "GO:0010965",]$description <- "regulation of mitotic sister chromatid separation"
trebo_topgo_bp_df_all[trebo_topgo_bp_df_all$GO.ID == "GO:0033047",]$description <- "regulation of mitotic sister chromatid segregation"
trebo_topgo_bp_df_all$term_category <- NA
trebo_topgo_bp_df_all[trebo_topgo_bp_df_all$description %in% vit_selection,]$term_category <- "vitamins"
trebo_topgo_bp_df_all[trebo_topgo_bp_df_all$description %in% lichens_selection,]$term_category <- "lichen"
trebo_topgo_bp_df_all[trebo_topgo_bp_df_all$description %in% essential_terms,]$term_category <- "essential"


ggplot(trebo_topgo_bp_df_all[trebo_topgo_bp_df_all$description %in% select_terms,], aes(x = bin_category, y = description, color = bin_category, size = -(logPval)))+
  geom_point(alpha = 0.6)+theme_bw()+facet_wrap(term_category~sample_id, scales = "free_y", ncol = 2)+
  scale_color_manual(values = bincat_colours)

vitaminplot <- ggplot(trebo_topgo_bp_df_all[trebo_topgo_bp_df_all$description %in% vit_selection,], aes(x = bin_category, y = description, color = bin_category, size = -(logPval)))+
  geom_point(alpha = 0.6)+theme_bw()+facet_wrap(~sample_id)+
  scale_color_manual(values = bincat_colours_trebo)+
  labs(x = "TE Distance Category", y = "GO Term Description", size = "-log(P-val)")+
  theme(axis.text.x = element_text(size = 14, angle = 60, hjust = 1), axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 18, face = "italic"),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18), axis.text.y = element_text(size = 14))
lichentermpgoplot <- ggplot(trebo_topgo_bp_df_all[trebo_topgo_bp_df_all$description %in% lichens_selection,], aes(x = bin_category, y = description, color = bin_category, size = -(logPval)))+
  geom_point(alpha = 0.6)+theme_bw()+facet_wrap(~sample_id)+
  scale_color_manual(values = bincat_colours_trebo)+
  labs(x = "TE Distance Category", y = "GO Term Description", size = "-log(P-val)")+
  theme(axis.text.x = element_text(size = 14, angle = 60, hjust = 1), axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 18, face = "italic"),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18), axis.text.y = element_text(size = 14))

trebo_topgo_bp_df_all[trebo_topgo_bp_df_all$GO.ID == "GO:2000816",]$description <- c("negative regulation of mitotic sister chromatid separation")
trebo_topgo_bp_df_all[trebo_topgo_bp_df_all$GO.ID == "GO:0033048",]$description <- c("negative regulation of mitotic sister chromatid segregation")

essentialgoplot <- ggplot(trebo_topgo_bp_df_all[trebo_topgo_bp_df_all$description %in% essential_terms,], aes(x = bin_category, y = description, color = bin_category, size = -(logPval)))+
  geom_point(alpha = 0.8)+theme_bw()+facet_wrap(~sample_id)+
  scale_color_manual(values = bincat_colours_trebo)+
  labs(x = "TE Distance Category", y = "GO Term Description", size = "-log(P-val)")+
  theme(axis.text.x = element_text(size = 14, angle = 60, hjust = 1), axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 18, face = "italic"),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18), axis.text.y = element_text(size = 14))
ggsave("figures/trebouxia_analysis/vitamingo_v2.pdf",
       vitaminplot, dpi = "retina", height = 6, width = 11)
ggsave("figures/trebouxia_analysis/lichengoterms_v2.pdf",
       lichentermpgoplot, dpi = "retina", height = 6, width = 11)
ggsave("figures/trebouxia_analysis/essentialgoterms_v2.1.pdf",
       essentialgoplot, dpi = "retina", height = 10, width = 11)



gluttransferase <- unique(trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$signature_accession %in% c("PF13417"),]$ipro_accesion)
rhodopsin <- unique(trebo_gene2te_quant_v2[grepl("rhodopsin",trebo_gene2te_quant_v2$ipro_description, ignore.case = FALSE),]$ipro_accesion)
tspombr <- unique(trebo_gene2te_quant_v2[grepl("tspO", trebo_gene2te_quant_v2$ipro_description, ignore.case = FALSE),]$ipro_accesion)
transporters_terms <- unique(trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$ipro_accesion %in% c("IPR003663"),]$ipro_accesion)
aquaporins <- c("IPR000425","IPR022357","IPR034294")
glycosylfamily8 <- unique(trebo_gene2te_quant_v2[grepl("Glycoside hydrolase, family 8", trebo_gene2te_quant_v2$ipro_description, fixed = TRUE),]$ipro_accesion)
glycosylfamily8 <- glycosylfamily8[-2]
chlorophytelichen_ipro <- c(gluttransferase, rhodopsin, tspombr, transporters_terms, aquaporins, glycosylfamily8)

trebo_functions_lichen <- trebo_gene2te_quant_v2[trebo_gene2te_quant_v2$ipro_accesion %in% chlorophytelichen_ipro,] %>%
  group_by(bin_cat, ipro_description, species_name) %>% summarize(count = n())
trebo_functions_lichen$bin_cat <- factor(trebo_functions_lichen$bin_cat, 
                                         levels = c("0d","1d1000",
                                                    "1001d5000","5001d10000",
                                                    "10001d20000"))
trebolicheniztion_tedist_plot <- ggplot(trebo_functions_lichen, aes(x = bin_cat, y = ipro_description, color = bin_cat, size = count))+
  geom_point(alpha = 0.6)+theme_bw()+facet_wrap(~species_name)+
  scale_color_manual(values = bincat_colours_trebo)+
  labs(x = "TE Distance Category", y = "Interpro Description", size = "Count")+
  theme(axis.text.x = element_text(size = 14, angle = 60, hjust = 1), axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 18, face = "italic"),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 18), axis.text.y = element_text(size = 14))
ggsave("figures/trebouxia_analysis/lichenization_tedist_v2.pdf",
       trebolicheniztion_tedist_plot, dpi = "retina", height = 6, width = 11)




# TE distance: genes
trebo_gene2te <- gene2te_bulk("stringtie_combo_final_separated/symbiochloris/gene/", 
             "stringtie_combo_final_separated/symbiochloris/earlgrey/") %>%
  drop_na(NearestTE)
colnames(trebo_gene2te)[2] <- "start_gff"
colnames(trebo_gene2te)[3] <- "end_gff"


lobaria_tracks <- read.gff("transcripts.gff3", GFF3 = TRUE)
lobaria_trebo_tracks <- lobaria_tracks[grepl("trebo",lobaria_tracks$seqid),]
lobaria_trebo_id <- data.frame(seqid = lobaria_trebo_tracks$seqid,
                               protid = str_extract(lobaria_trebo_tracks$attributes, "glLobPulm2_\\d{2,};"))
lobaria_trebo_id$protid <- substr(lobaria_trebo_id$protid, 1, nchar(lobaria_trebo_id$protid)-1)
lobaria_interproscan <- get_iproscan("interproscan_combination/Lobaria_pulmonaria.interproquery.faa.tsv")
lobaria_trebo_interproscan <- merge(lobaria_trebo_id,lobaria_interproscan, by.y = "Target_ID", by.x = "protid") %>% distinct()

trebo_gene2te_ipro <- merge(lobaria_trebo_interproscan, trebo_gene2te, by.y = "seqnames", by.x = "seqid") %>% distinct()

# contig plots
symbiochloris_fasta <- readDNAStringSet("chlorophyte_fasta/glRicVire11_bin41.fna")
symbiochloris_bed <- data.frame("seq" = names(symbiochloris_fasta), "start" = 1, end = width(symbiochloris_fasta)) %>% arrange(desc(end))
symbiochloris_contigcount <- length(symbiochloris_fasta)
symbiochlorisfasta_newnames <- structure(paste0("contig_", seq(from = 1, to = symbiochloris_contigcount)), names = symbiochloris_bed$seq)
symbiochloris_bed$seq <- symbiochlorisfasta_newnames[symbiochloris_bed$seq]

symbiochloris_gdf <- as.data.frame(get_bkg(symbiochloris_fasta, k = 1, RC = F, window = TRUE, window.size = 10000, merge.res = FALSE))
symbiochloris_c_df <- symbiochloris_gdf[symbiochloris_gdf$klet == "C",]
symbiochloris_gg_df <- symbiochloris_gdf[symbiochloris_gdf$klet == "G",]
symbiochloris_gc_df <- data.frame(seq = symbiochloris_c_df$sequence, start = symbiochloris_c_df$start, end = symbiochloris_c_df$stop, score = (symbiochloris_c_df$probability + symbiochloris_gg_df$probability)/2)
symbiochloris_gc_df$seq <- symbiochlorisfasta_newnames[symbiochloris_gc_df$seq]

gff_colnames <- c("seq","source","feature","start","end","score","strand","frame","attribute")
symbiochloris_tegff <- read.table("earlgrey/chlorophyta/te_gff/glRicVire11_bin.41.filteredRepeats.gff")
colnames(symbiochloris_tegff) <- gff_colnames
symbiochloris_tegff$seq <- symbiochlorisfasta_newnames[symbiochloris_tegff$seq]
symbiochloris_te_bed <- symbiochloris_tegff[,c(1,4,5,6)]

symbiochloris_telomeres <- read.delim("glRicVire11_bin.41.fa_telomeres.txt", header = F)[,c(1,2)]
symbiochloris_telomeres$V1 <- symbiochlorisfasta_newnames[symbiochloris_telomeres$V1]
symbiochloris_telomeres <- symbiochloris_telomeres %>% na.omit()
colnames(symbiochloris_telomeres) <- c("seq","direction")
symbiochloris_telomeres$position <- 0 
symbiochloris_telomeres$start <- symbiochloris_bed[symbiochloris_bed$seq  %in% symbiochloris_telomeres$seq,]$start
symbiochloris_telomeres$end <- symbiochloris_bed[symbiochloris_bed$seq  %in% symbiochloris_telomeres$seq,]$end
symbiochloris_telomeres[symbiochloris_telomeres$direction == "forward",]$position <- symbiochloris_telomeres[symbiochloris_telomeres$direction == "forward",]$start
symbiochloris_telomeres[symbiochloris_telomeres$direction == "reverse",]$position <- symbiochloris_telomeres[symbiochloris_telomeres$direction == "reverse",]$end

symbiochloris_gc_df$telomere <- "no"
symbiochloris_gc_df[symbiochloris_gc_df$seq %in% symbiochloris_telomeres$seq,]$telomere <- "yes"
symbiochloris_contig_plot_df <-left_join(symbiochloris_gc_df, symbiochloris_telomeres[,c(1,5)])

symbiochloris_contig_plot <- ggplot(symbiochloris_contig_plot_df, aes(x = start, y = reorder(seq, end), fill = score*100))+geom_tile()+
  geom_point(aes(x = position, y = seq), color = "black", pch = "*" , size = 14)+
  theme_bw()+scale_fill_viridis_c()+
  scale_x_continuous(breaks = c(0,1000000, 2000000, 3000000, 4000000, 5000000,6000000, 7000000,8000000),
                     labels = c("0","1", "2","3", "4", "5", "6", "7", "8"))+
  labs(x = "Contig Length (Mb)", y = "Contig ID", fill = "GC Content (%)")+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 18),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18))

ggsave("figures/trebouxia_analysis/symbiochloris_contigs.pdf",symbiochloris_contig_plot,
       dpi = "retina", height = 6, width = 8)
# Salmon
allsalmon <- bulktximport("salmon_repfungi_combo/", "trebouxoid_gene2tx/", txOut = TRUE)
species <- names(allsalmon)

Quant <- allsalmon
QuantM <- Quant
for (i in 1:length(species)) {
  speciesid <- species[i]
  colnames(QuantM[[i]]$counts) <- c("Rep1", "Rep2")
  QuantM[[i]] <- as.data.frame(QuantM[[i]]$counts)
  QuantM[[i]]$Transcript <- rownames(QuantM[[i]])
  QuantM[[i]] <- QuantM[[i]][, c(3, 1, 2)]
  QuantM[[i]]$species <- speciesid
}
QuantM <- do.call(rbind, QuantM)
rownames(QuantM) <- NULL

QuantM_longer <- pivot_longer(QuantM,cols = c("Rep1","Rep2"))
QuantM_longer$new_id <- paste0(QuantM_longer$species,"_",QuantM_longer$name)
QuantM_longer$species <- NULL
QuantM_longer$name <- NULL

QuantM_counts <- pivot_wider(QuantM_longer, names_from  = new_id, values_from = value, values_fill =  0) %>% as.data.frame()
rownames(QuantM_counts) <- QuantM_counts$Transcript
QuantM_counts$Transcript <- NULL
QuantM_counts[[1]] <- round(QuantM_counts[[1]])
QuantM_counts[[2]] <- round(QuantM_counts[[2]])
QuantM_counts[[3]] <- round(QuantM_counts[[3]])
QuantM_counts[[4]] <- round(QuantM_counts[[4]])

group <- factor(c("Lobaria","Lobaria","Ricasolia","Ricasolia"))

# differential expression for transcriptomic analysis
y <- DGEList(counts=QuantM_counts,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- normLibSizes(y)
design <- model.matrix(~group)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)

qlf_df <- qlf$table
qlf_df$direction <- "Non-significant"
qlf_df[qlf_df$PValue <= 0.05 & qlf_df$logFC >= 1,]$direction <- "Positive"
qlf_df[qlf_df$PValue <= 0.05 & qlf_df$logFC <= -1,]$direction <- "Negative"
qlf_df$Target_ID <- rownames(qlf_df)

trebouxoid_degenes_volcano <- ggplot(qlf_df, aes(x = logFC, y = -log(PValue, 10), color = direction))+
  geom_point()+scale_color_manual(values = c("#bf212f", "grey75","#27b376"))+
  theme_bw()+labs(color = "Significance", y = "-log(p-value)", x = "logFC")+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18), legend.title = element_text(size = 18),
        legend.text = element_text(size = 18))

ggsave("figures/trebouxia_analysis/de_volcano.pdf",
       trebouxoid_degenes_volcano, dpi = "retina", height = 5, width = 8)

trebouxoid_interproscan <- get_iproscan("Trebouxoid_symbiont.combined.interproquery.faa.tsv")
qlf_df_annotations <- merge(qlf_df, trebouxoid_interproscan)
qlf_df_annotations_sig <- qlf_df_annotations[qlf_df_annotations$PValue <= 0.05 & abs(qlf_df_annotations$logFC) >= 1,]

trebouxoid_interpro_select <- c("IPR024690","IPR004232",
                               "IPR004307",
                               "IPR000425","IPR011614",
                               "IPR001425","IPR003663")
trebouxoid_interpro_select_anno <- c("Nitrile hydratase","Nitrile hydratase",
                                     "TspO/MBR family","Major intrinsic protein",
                                     "Catalase","Rhodopsins","Sugar/inositol transporter")

trebouxoid_interpro_select_df <- data.frame("accession" = trebouxoid_interpro_select, "description" = trebouxoid_interpro_select_anno)
trebouxoid_interpro_select_geneid <- data.frame("Transcript" = trebouxoid_interproscan[trebouxoid_interproscan$ipro_accesion %in% trebouxoid_interpro_select,]$Target_ID, "accession" = trebouxoid_interproscan[trebouxoid_interproscan$ipro_accesion %in% trebouxoid_interpro_select,]$ipro_accesion)
trebouxoid_interpro_select_geneid_df <- merge(trebouxoid_interpro_select_geneid, trebouxoid_interpro_select_df)
QuantM_longer_interpro <- merge(QuantM_longer, trebouxoid_interpro_select_geneid_df)

QuantM_longer_interpro[QuantM_longer_interpro$log2 == -Inf,]$log2 <- 0
ggplot(QuantM_longer_interpro, aes(x = new_id, y = Transcript, fill = log2, color = description))+geom_tile()+
  scale_color_manual(values = alt8palette)+theme_bw()

# topGO enrichment analysis
trebouxoid_golist <- vector("list", length = length(unique(qlf_df_annotations$Target_ID)))
names(trebouxoid_golist) <- unique(qlf_df_annotations$Target_ID)

trebouxoid_golist <- vector("list", length = length(unique(qlf_df_annotations_sig$Target_ID)))
names(trebouxoid_golist) <- unique(qlf_df_annotations_sig$Target_ID)
for (i in 1:length(trebouxoid_golist)){
  transcript_id <- names(trebouxoid_golist)[i]
  go_terms <- qlf_df_annotations_sig[qlf_df_annotations_sig$Target_ID == transcript_id,]$go_term
  go_terms <- unlist(strsplit(go_terms, "|", fixed = TRUE))
  trebouxoid_golist[[i]] <- go_terms[!grepl("-",go_terms)]
}
trebouxoid_gene2go_bp <- annFUN.gene2GO(whichOnto = "BP", gene2GO = trebouxoid_golist)
trebouxoid_gene2go_cc <- annFUN.gene2GO(whichOnto = "CC", gene2GO = trebouxoid_golist)
trebouxoid_gene2go_mf <- annFUN.gene2GO(whichOnto = "MF", gene2GO = trebouxoid_golist)

trebouxoid_negativede_transcripts <- qlf_df_annotations_sig[qlf_df_annotations_sig$direction == "negative",]$PValue
trebouxoid_positivede_transcripts <- qlf_df_annotations_sig[qlf_df_annotations_sig$direction == "positive",]$PValue
names(trebouxoid_negativede_transcripts) <- qlf_df_annotations_sig[qlf_df_annotations_sig$direction == "negative",]$Target_ID
names(trebouxoid_positivede_transcripts) <- qlf_df_annotations_sig[qlf_df_annotations_sig$direction == "positive",]$Target_ID


trebouxoid_negativede_topgo_bp <- calc_topgo_fisher(trebouxoid_negativede_transcripts, "BP")
trebouxoid_positivede_topgo_bp <- calc_topgo_fisher(trebouxoid_positivede_transcripts, "BP")



                                                                             