# input data available upon request from authors for re-running code

# load data
qs50_mapping_tax_df <- read.csv("qs50_mapping_tax_df.csv")
qs50_mapping_tax_df$phy <- factor(qs50_mapping_tax_df$Phylum, levels = phylum_taxonomiclevels)
qs50_mapping_tax_df$lichen_species <- names(lichen_species)[match(qs50_mapping_tax_df$sample_id,lichen_species)]
qs50_mapping_tax_df$lichen_species <- factor(qs50_mapping_tax_df$lichen_species, levels = rev(lichen_species_fac))
qs50_mapping_tax_df$lichentree <- names(mycobiont_mag)[match(qs50_mapping_tax_df$lichen_species,mycobiont_mag)]
qs50_mapping_tax_df$lichentree <- factor(qs50_mapping_tax_df$lichentree, levels = lichen_tree_order_fa)

qs50_mapping_cyano_df <- read.csv("qs50_cyanobacteriadepth.csv")
qs50_mapping_cyano_df$phy <- factor(qs50_mapping_cyano_df$Phylum, levels = phylum_taxonomiclevels)
qs50_mapping_cyano_df$lichen_species <- names(lichen_species)[match(qs50_mapping_cyano_df$sample_id,lichen_species)]
qs50_mapping_cyano_df$lichen_species <- factor(qs50_mapping_cyano_df$lichen_species, levels = rev(lichen_species_fac))
qs50_mapping_cyano_df$new_bin_id <- factor(qs50_mapping_cyano_df$new_bin_id, levels = cyano_mag)
qs50_mapping_cyano_df$lichentree <- names(mycobiont_mag)[match(qs50_mapping_cyano_df$lichen_species,mycobiont_mag)]
qs50_mapping_cyano_df$lichentree <- factor(qs50_mapping_cyano_df$lichentree, levels = lichen_tree_order_fa)

cyanototalprops <- qs50_mapping_tax_df[qs50_mapping_tax_df$Phylum == "Cyanobacteriota",] %>% group_by(sample_id) %>% summarize(lichen_species = lichen_species, total_cyano_prop = sum(prop)) %>% distinct()
cyanototalprops$lichentree <- names(mycobiont_mag)[match(cyanototalprops$lichen_species,mycobiont_mag)]
cyanototalprops$lichentree <- factor(cyanototalprops$lichentree, levels = lichen_tree_order_fa)

qs50_mapping_summary <- qs50_mapping_tax_df %>% group_by(sample_id, phy) %>%
  summarize(count = n(),lichentree = lichentree,lichen_species = lichen_species, taxtype = taxtype) %>% distinct()
qs50_mapping_summary_cyano <- qs50_mapping_cyano_df %>% group_by(sample_id, new_bin_id) %>%
  summarize(count = n(),lichentree = lichentree,lichen_species = lichen_species, phy = phy, Family = Family)

# qs50 mapping results
qs50_summary_plot <- ggplot(qs50_mapping_summary, aes(y = lichentree, x = count, fill = phy, color = phy))+
  geom_bar(position = "stack", stat = "identity", color = "black")+scale_fill_manual(values = phylum_colours)+
  scale_color_manual(values = phylum_colours, guide = "none")+
  theme_bw()+scale_x_continuous(limits = c(0,35), expand = c(0,0))+
  labs(x = "Count", y = "Lichen Taxa", fill = "Phylum")+
  theme(axis.text.y = element_text(face = "italic", size = 16), axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20), legend.text = element_text(size = 16))

qs50_summary_plot_nolegend <- qs50_summary_plot+theme(legend.position = "none")

qs50_mapping_plot_horizontal <- ggplot(qs50_mapping_tax_propsummary, aes(y = lichentree, x = prop, fill = phy, color = phy))+
  geom_bar(stat = "identity", position = "stack", color = "black")+
  scale_color_manual(values = phylum_colours, guide = "none")+
  scale_fill_manual(values = phylum_colours)+theme_bw()+
  scale_x_continuous(expand = c(0,0), limits = c(0,103))+
  labs(x = "Relative Abundance (%)", y = "Lichen Taxa", fill = "Phylum")+
  theme(axis.text.y = element_text(face = "italic", size = 16), axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20), legend.text = element_text(size = 16))

qs50_mapping_plot_horizontal_depth <- ggplot(qs50_mapping_tax_df, aes(y = lichen_species, x = bin_depth, fill = phy, color = phy))+
  geom_bar(stat = "identity", position = "stack")+
  scale_color_manual(values = phylum_colours, guide = "none")+
  scale_fill_manual(values = phylum_colours)+theme_bw()+
  scale_x_continuous(expand = c(0,0), limits = c(0,1300))+
  labs(x = "Mean Depth", y = "Lichen Taxa", fill = "Phylum")+
  theme(axis.text.y = element_text(face = "italic", size = 16), axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20),
        legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("figures/microbial_community/qs50_mapping_horizontal.pdf", 
       qs50_mapping_plot_horizontal, 
       dpi = "retina", height = 6, width = 12)

ggsave("figures/microbial_community/qs50_mapping_horizontal_depth.pdf", 
       qs50_mapping_plot_horizontal_depth, 
       dpi = "retina", height = 6, width = 12)

ggsave("figures/microbial_community/qs50_mapping_horizontal_depth.pdf", 
       qs50_mapping_plot_horizontal_depth, 
       dpi = "retina", height = 6, width = 12)

ggsave("figures/microbial_community/qs50_summary_magcount.pdf", 
       qs50_summary_plot, 
       dpi = "retina", height = 6, width = 12)



qs50_mapping_tax_propsummary <- qs50_mapping_tax_df %>% group_by(lichen_species, phy) %>% summarize(prop = sum(prop), lichentree = lichentree) %>% distinct()
write.csv(qs50_mapping_tax_propsummary, "qs50_mapping_tax_propsummary.csv", row.names = FALSE, quote = FALSE)

qs50_bacteriamags_depthabunplot <- ggplot(qs50_mapping_tax_df[qs50_mapping_tax_df$taxtype == "bacteria" & qs50_mapping_tax_df$Phylum != "Cyanobacteriota",], aes(y = lichen_species, x = bin_id, fill = prop))+
  geom_tile()+facet_wrap(~Phylum, scales = "free")+scale_fill_gradient(low = "white", high = "navy")+theme_bw()+
  labs(x = "MAG ID", y = "Lichen Taxa", fill = "Relative Abundance (%)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16), axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, face = "italic"), axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("figures/microbial_community/qs50_bacteriamags_depthabunplot.pdf",
       qs50_bacteriamags_depthabunplot,
       dpi = "retina", height = 10, width = 15)

qs50_cyanobacteriamags_summary_plot <- ggplot(qs50_mapping_summary_cyano, aes(y = lichentree, x = new_bin_id, fill = new_bin_id), color = "black")+
  geom_tile(color = "black")+geom_point(aes(x = new_bin_id, y = lichentree, color = new_bin_id), shape = 8)+
  theme_bw()+labs(x = "MAG ID", y = "Lichen Taxa")+
  scale_color_manual(values = c("#50336b","#7b6690","#a69ab5","#d3ccda",rep("black",6)))+
  scale_fill_manual(values = cyanobacteria_bin_colours)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 16, face = "italic"),
        legend.position = "none", axis.title.x = element_text(size = 20))

qs50_cyanobacteriamags_totalprop_plot <- ggplot(cyanototalprops, aes(y = lichentree, x = total_cyano_prop))+
  geom_bar(stat = "identity", fill = "#a2ceaa", color = "black")+theme_bw()+scale_x_continuous(expand = c(0,0), limits = c(0,103))+
  labs("Relative Abundance (%)")+
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16), axis.text.y = element_blank(),
        legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.position = "none")+
  labs(x = "Relative Abundance (%)")


qs50_cyanobacteriamags_plot <- ggplot(qs50_mapping_cyano_df, aes(y = lichentree , x = prop, fill = new_bin_id))+
  geom_bar(stat = "identity", position = "stack", color = "black")+theme_bw()+
  scale_x_continuous(expand = c(0,0), limits = c(0,103))+
  scale_fill_manual(values = cyanobacteria_bin_colours)+
  labs(x = "Relative Abundance (%)", y = "Lichen Taxa", fill = "MAG ID")+
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16, face = "italic"),
        legend.title = element_text(size = 20), legend.text = element_text(size = 16))

qs50_cyanobacteriamags_plot2 <- ggplot(qs50_mapping_cyano_df, aes(y = lichentree , x = prop, fill = new_bin_id))+
  geom_bar(stat = "identity", position = "stack", color = "black")+theme_bw()+
  scale_x_continuous(expand = c(0,0), limits = c(0,103))+
  scale_fill_manual(values = cyanobacteria_bin_colours)+
  labs(x = "Relative Abundance (%)", y = "Lichen Taxa", fill = "MAG ID")+
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16, face = "italic"),
        legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.position = "none")

qs50_cyanobacteriamags_plot3 <- ggplot(qs50_mapping_cyano_df, aes(y = lichentree , x = prop, fill = new_bin_id))+
  geom_bar(stat = "identity", position = "stack", color = "black")+theme_bw()+
  scale_x_continuous(expand = c(0,0), limits = c(0,103))+
  scale_fill_manual(values = cyanobacteria_bin_colours)+
  labs(x = "Relative Abundance (%)", y = "Lichen Taxa", fill = "MAG ID")+
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16), axis.text.y = element_blank(),
        legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.position = "none")

qs50_prop_summary_cyano_plot <- plot_grid(qs50_cyanobacteriamags_plot3, qs50_cyanobacteriamags_totalprop_plot,rel_widths = c(2,1.5))
legend <- get_legend(qs50_cyanobacteriamags_plot)
qs50_prop_summary_cyano_plot_legend <- plot_grid(qs50_prop_summary_cyano_plot, legend)

plot_grid(qs50_cyanobacteriamags_totalprop_plot, legend)
qs50_cyanobacteriamags_depth_plot <- ggplot(qs50_mapping_cyano_df, aes(y = lichen_species , x = bin_depth, fill = new_bin_id2))+
  geom_bar(stat = "identity", position = "stack")+theme_bw()+
  scale_x_continuous(expand = c(0,0), limits = c(0,1050))+
  scale_fill_manual(values = cyanobacteria_bin_colours)+
  labs(x = "Mean Depth", y = "Lichen Taxa", fill = "MAG ID")+
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16, face = "italic"),
        legend.title = element_text(size = 20), legend.text = element_text(size = 16))

ggsave("figures/microbial_community/qs50_cyanomapping_horizontal.pdf", 
       qs50_cyanobacteriamags_plot, 
       dpi = "retina", height = 6, width = 10)

ggsave("figures/microbial_community/qs50_cyanomapping_depths.pdf", 
       qs50_cyanobacteriamags_depth_plot, 
       dpi = "retina", height = 6, width = 10)


lichen_primarytree <- read.tree("trees/lichen_primarytree.txt")
lichen_tree_order_fa <- lichen_primarytree$tip.label
lichen_tree_order <- gsub("[_]bin*[.]fa","",lichen_primarytree$tip.label)
lichen_photobiont_status <- data.frame(sample_id = lichen_primarytree$tip.label, status = c(rep("bipartite",9), rep("tripartite", 2)))

lichen_photobiont_status$speciesname <- c("Leptogium burgessi", "Nephroma laevigatum", "Peltigera hymenina", "Peltigera horizontalis", "Peltigera collina", "Peltigera membranacea", "Peltigera praetextata", "Pseudocyphellaria norvegica", "Sticta sylvatica", "Ricasolia virens", "Lobaria pulmonaria")

lichen_primarytree_metadata <- ggtree(lichen_primarytree, branch.length = "none") %<+% lichen_photobiont_status

lichen_primarytree_ggtree_colouredtips <- lichen_primarytree_metadata+geom_tiplab(aes(label = speciesname, color = status),align = TRUE, linetype = "dotted", linesize = 0.1, fontface = 3)+
  scale_color_manual(values = c("#00584f","#00baa6"))+theme( legend.position = "none", coord_cartesian(clip = "off"))+xlim(NA, 25)+ylim(NA,100)

qs50_summary_plot_noylab <- qs50_summary_plot+labs(y = NULL, x = "MAG Count")+theme(axis.text.y = element_blank(),legend.position = "none")
qs50_mapping_plot_noylab <- qs50_mapping_plot_horizontal+labs(y = NULL)+theme(axis.text.y = element_blank())
qs50_treeplot <- qs50_summary_plot_noylab %>% insert_left(lichen_primarytree_ggtree_colouredtips) %>% insert_right(qs50_mapping_plot_noylab)
ggsave("figures/microbial_community/qs50_mapping_tree.pdf", qs50_treeplot, dpi = "retina", height = 6, width = 15)

cyanoqs50_mag_summary_plot <- qs50_cyanobacteriamags_summary_plot+labs(y = NULL, x = "MAG Count")+theme(axis.text.y = element_blank())+labs(x = "MAG Presence")
cyanoqs50_mapping_plot_noylab <- qs50_cyanobacteriamags_plot +labs(y = NULL, x = "Relative Abundance (%)")+theme(axis.text.y = element_blank())
cyano_tree_plot <- cyanoqs50_mag_summary_plot %>% insert_left(lichen_primarytree_ggtree_colouredtips) %>% insert_right(cyanoqs50_mapping_plot_noylab) %>% insert_right(qs50_cyanobacteriamags_totalprop_plot, width = 0.6)
ggsave("figures/microbial_community/qs50_cyanomapping_tree.pdf", cyano_tree_plot, dpi = "retina", height = 6, width = 16)


# KEGG nitrogen fixation + transposase plots
all_ko <- read.table("all_mags_kofam.txt", fill = TRUE)
metamdbg_gtdbtk_tax_df <- read.csv("metamdbg_gtdbtkeuk_tax_df.csv")
all_ko <- all_ko[all_ko$V2 != "",]
colnames(all_ko) <- c("seq", "KEGG")

splitmagname <- strsplit(all_ko$seq, "_seq")
magnames_ko <- vector("character", length = length(splitmagname))

for (i in 1:length(splitmagname)){
  magnames_ko[i] <- splitmagname[[i]][1]
}
all_ko$genome <- magnames_ko
all_ko_tax <- merge(all_ko, metamdbg_gtdbtk_tax_df)

# nitrogen fixation KO id
vnf_kegg <- c("K22896","K22897","K22898","K22899")
n2_kegg <- c("K02588", "K02586", "K02591")

cyano_ko <- all_ko_tax[all_ko_tax$Phylum == "Cyanobacteriota",]
nitrogen_cyanos <- cyano_ko[cyano_ko$KEGG %in% c(vnf_kegg, n2_kegg),]
nitrogen_cyanos$style <- "Mn-Fe"
nitrogen_cyanos[nitrogen_cyanos$KEGG %in% vnf_kegg,]$style <- "Vanandium"

nitrogen_cyanos$genome <- factor(nitrogen_cyanos$genome, 
                                 levels = rev(c("glLepBurg3_bin.4","glStiSylv4_bin.1",
                                                "glLepBurg3_bin.16","glPelMemb1_bin.1",
                                                "glPelPrae3_bin.3","glPelColl1_bin.21",
                                                "glStiSylv4_bin.6","glPseNorv1_bin.12",
                                                "glLepBurg3_bin.18")))
nitrogen_cyanos$new_bin_id <- factor(nitrogen_cyanos$new_bin_id,
                                     levels = rev(c("Nostoc_3","Nostoc_4",
                                                "Nostoc_2","Nostoc_6",
                                                "Nostoc_5","Nostoc_1")))

nitrogen_cyanotype_plot <- ggplot(nitrogen_cyanos, aes(x = KEGG, y = new_bin_id, fill = style))+
  geom_tile(color = "black")+facet_wrap(~style, scales = "free_x")+theme_bw()+
  scale_fill_manual(values = c("#364F3E","#D3A577"))+
  theme(axis.text.y = element_text(size = 14, face = "italic"), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 14), legend.title = element_text(size = 18),
        strip.text = element_text(size = 14))+
  labs(fill = "Pathway Type", x = "KEGG Orthology ID")

ggsave("figures/cyanobacteria/nitrogen_fix_cyanos_plot.pdf",
       nitrogen_cyanotype_plot, dpi = "retina", height = 5, width = 9)

# transposases
transposase_ko <- c("K01152", "K07481","K07482","K07483","K07484",
                    "K07485","K07486","K07487","K07488","K07489",
                    "K07489","K07491","K07492","K07493","K07494",
                    "K07495","K07496","K07497","K07498","K07499",
                    "K18320","K23209")

transposase_ko_df <- all_ko_tax[all_ko_tax$KEGG %in% transposase_ko,]
transposase_ko_df_summary_phylum <- transposase_ko_df %>% group_by(Phylum) %>% summarize(count = n())
transposase_ko_df_summary_mag <- transposase_ko_df %>% group_by(genome) %>% summarize(count = n(), Phylum = Phylum, Family = Family) %>% distinct()

transposase_ko_df_summary_mag$Phylum <- factor(transposase_ko_df_summary_mag$Phylum, 
                                               levels = c("Pseudomonadota","Myxococcota",
                                                          "Acidobacteriota","Bacteroidota",
                                                          "Verrucomicrobiota","Planctomycetota",
                                                          "Actinomycetota","Chloroflexota","Armatimonadota",
                                                          "Cyanobacteriota"))
transposase_phylum_summplot <- ggplot(transposase_ko_df_summary_mag, aes(y = Phylum, x = count, fill = Phylum))+
  geom_boxplot()+scale_fill_manual(values = phylum_colours)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 18), legend.text = element_text(size = 14))+
  labs(y = "Phylum", x = "Transposase Count", fill = "Phylum")

cyanotransposase_phylum_summplot <- ggplot(transposase_ko_df_summary_mag[transposase_ko_df_summary_mag$Phylum == "Cyanobacteriota",], aes(y = Family, x = count, fill = Family))+
  geom_boxplot()+scale_fill_manual(values = c("#6fb899","#8175aa"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 18), legend.text = element_text(size = 14))+
  labs(y = "Family", x = "Transposase Count", fill = "Family")


ggsave("figures/cyanobacteria/transposase_count_plot.pdf",
       transposase_phylum_summplot, dpi = "retina", height = 5, width = 8)

ggsave("figures/cyanobacteria/cyanotransposase_count_plot.pdf",
       cyanotransposase_phylum_summplot, dpi = "retina", height = 5, width = 8)
######
checkm_drep <- read.csv("bin_stats_ext.tsv", sep = "\t", header = FALSE)
colnames(checkm_drep) <- c("bin_id","checkm")
checkm_drep$comp <- str_extract(checkm_drep$checkm, "'Completeness':\\s[0-9]{1,}.[0-9]{1,}") %>% str_extract("[0-9]{1,}.[0-9]{1,}")
checkm_drep$contam <- str_extract(checkm_drep$checkm, "'Contamination':\\s[0-9]{1,}.[0-9]{1,}") %>% str_extract("[0-9]{1,}.[0-9]{1,}")
checkm_drep$qs50 <- as.numeric(checkm_drep$comp) - (5*as.numeric(checkm_drep$contam))
qs50_bacteria <- unique(checkm_drep[checkm_drep$qs50 >= 50,]$bin_id)
writeLines(qs50_bacteria, "qs50_bacteriasamples.txt")






