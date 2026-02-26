library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(aplot)

`%notin%` <- Negate(`%in%`)
# lichen primary mycobiont tree
lichen_primarytree <- read.tree("data/trees/lichen_primarytree.txt")
lichen_tree_order_fa <- lichen_primarytree$tip.label
lichen_tree_order <- gsub("[_]bin*[.]fa","",lichen_primarytree$tip.label)
lichen_photobiont_status <- data.frame(sample_id = lichen_primarytree$tip.label, status = c(rep("bipartite",9), rep("tripartite", 2)))

lichen_photobiont_status$speciesname <- c("Leptogium burgessi", "Nephroma laevigatum", "Peltigera hymenina", "Peltigera horizontalis", "Peltigera collina", "Peltigera membranacea", "Peltigera praetextata", "Pseudocyphellaria norvegica", "Sticta sylvatica", "Ricasolia virens", "Lobaria pulmonaria")

lichen_primarytree_metadata <- ggtree(lichen_primarytree, branch.length = "none") %<+% lichen_photobiont_status

lichen_primarytree_ggtree_colouredtips <- lichen_primarytree_metadata+geom_tiplab(aes(label = speciesname, color = status),align = TRUE, linetype = "dotted", linesize = 0.1, fontface = 3)+
  scale_color_manual(values = c("#00584f","#00baa6"))+theme( legend.position = "none", coord_cartesian(clip = "off"))+xlim(NA, 25)+ylim(NA,100)

# chlorophyte photobiont 18S tree
chlorophyte_tree <- read.tree("data/trees/18s_chlorophytetree_midroot.txt")
chlorophyte_18s_accessions <- chlorophyte_tree$tip.label
chlorophyte_18s_accessions[4] <- ""
chlorophyte_tree_taxanames <- c("Trentepohlia arborum strain NIES-2648","Trentepohlia aurea strain NIES-2652",
                                "Symbiochloris tschermakiae SAG 46.85", "Trebouxiophyceae sp. symbiont",
                                "Symbiochloris reticulata strain SAG 53.87","Symbiochloris tropica strain CAUP H8602",
                                "Symbiochlroris irregularis strain SAG 2036", "Symbiochloris pauciautosporica strain SAG 12.86",
                                "Dictyochloropsis splendida strain SAG 244.80", "Symbiochloris handae strain CCHU 5616",
                                "Trebouxia erici strain IAM C-593", "Trebouxia sp. strain CCAP 213/3",
                                "Trebouxia sp. strain UTEX SNO74", "Trebouxia sp. strain SAG 2463",
                                "Elliptochloris bilobata strain SAG 245.80","Coccomyxa viridis strain SAG 216-14",
                                "Coccomyxa viridis strain SAG 216-1", "Coccomyxa sp. strain CCAP216/24",
                                "Coccomyxa simplex strain SAG 216-8", "Coccomyxa simplex strain SAG 216-5",
                                "Coccomyxa simplex strain SAG 216-11a")
chlorophyte_habitat <- c("Terrestrial","Terrestrial", "Lichen - Photobiont", "Me","Lichen - Photobiont",
                         "Terrestrial","Terrestrial","Lichen - Photobiont","Lichen - Photobiont","Lichen - Photobiont",
                         "Lichen - Photobiont", "Lichen - Photobiont","Terrestrial","Unspecified","Terrestrial",
                          "Lichen - Non-photobiont", "Lichen - Non-photobiont","Terrestrial", "Unspecified",
                         "Lichen - Photobiont","Lichen - Photobiont")
chlorophyte_tree$tip.label <- paste0(chlorophyte_18s_accessions, " " ,chlorophyte_tree_taxanames)
chlorophyte_tree_order <- chlorophyte_tree$tip.label
chlorophyte_tree_meta_df <- data.frame(sample_id = chlorophyte_tree_order, habitat = chlorophyte_habitat)
chlorophyte_tree_metadata <- ggtree(chlorophyte_tree) %<+% chlorophyte_tree_meta_df

chlorophyte_tree_colours <- c("Terrestrial" = "grey35","Unspecified" = "grey80", 
                              "Me" = "#bfe1a8","Lichen - Photobiont" = "#4d8769", 
                              "Lichen - Non-photobiont" = "#adcfb3",
                              "Aquatic" = "#637dbd","Pathogen" = "#836197",
                              "Endosymbiont" = "#ae95bc","Ice" = "#bac7e7",
                              "Peatland" = "#807d72")
chlorophyte_tree_ggtree_colouredtips <- chlorophyte_tree_metadata+geom_tiplab(aes(color = habitat),align = TRUE, linetype = "dotted", linesize = 0.1, fontface = 3)+
  scale_color_manual(values = chlorophyte_tree_colours)+
  theme( legend.position = "bottom")+xlim(NA,0.45)+coord_cartesian(clip = "off")+
  geom_treescale(y = 20, x = 0, width = 0.05)

ggsave("figures/trebouxia_analysis/trebouxia_18s_ggtree.pdf",
       chlorophyte_tree_ggtree_colouredtips, width = 8, height = 6, dpi = "retina")

# chlorophyte BUSCO tree
chlorophyte_busco_tree_unroot <- read.tree("data/trees/chlorophyte_busco.tree")
chlorophyte_busco_accessions <- chlorophyte_busco_tree_unroot$tip.label
chlorophyte_busco_taxa <- c("Coccomyxa subellipsoidea C-169","Micractinium conductrix strain SAG 241.80",
                            "Chlorella vulgaris strain CCAP 211/11P", "Chlorella vulgaris PKVL7422",
                            "Chlorella sorokiniana", "Chlorella sorokiniana strain UTEX 1602",
                            "Chlorella sorokiniana strain DSCG149","Chlorella sorokiniana",
                            "Chlorella sorokiniana strain DSCG147","Chlorella sorokiniana",
                            "Chlorella sorokiniana","Chlorella sorokiniana strain DSCG150",
                            "Picochlorum renovo","Picochlorum sp. BPE23",
                            "Picochlorum sp. BPE23","Picochlorum sp. 'celeri'",
                            "Picochlorum sp. SENEW3","Picochlorum sp. SENEW3",
                            "Nannochloris desiccata UTEX 2526","Nannochloris desiccata UTEX 2437",
                            "Marvania coccoides","Prototheca wickerhamii",
                            "Prototheca wickerhamii","Medakamo hakoo",
                            "Diplosphaera chodatii","Tetratostichococcus sp. P1",
                            "Deuterostichococcus epilithicus","Asterochloris sp. CNOR1",
                            "Trebouxiophyceae sp. symbiont","Coccomyxa viridis",
                            "Coccomyxa viridis strain SAG 216-4", "Coccomyxa viridis",
                            "Coccomyxa sp. Obi"
                            )
chlorophyte_busco_habitats <- c("Terrestrial","Endosymbiont", 
                               "Aquatic","Unspecified",
                               "Aquatic", "Unspecified",
                               "Aquatic","Unspecified",
                               "Aquatic","Unspecified",
                               "Unspecified","Aquatic",
                               "Unspecified","Aquatic",
                               "Aquatic","Aquatic",
                               "Aquatic","Aquatic",
                               "Aquatic","Aquatic",
                               "Aquatic","Pathogen",
                               "Pathogen","Unspecified",
                               "Lichen - Photobiont","Peatland",
                               "Ice","Lichen - Photobiont",
                               "Me","Lichen - Non-photobiont",
                               "Lichen - Photobiont","Lichen - Photobiont",
                               "Aquatic")

chlorophyte_busco_tree_meta_df <- data.frame(taxa = chlorophyte_busco_taxa, habitat = chlorophyte_busco_habitats, accessions = chlorophyte_busco_accessions)
chlorophyte_busco_tree <- read.tree("data/treeschlorophyte_busco_tree.txt")
chlorophyte_busco_tree_order <- chlorophyte_busco_tree$tip.label
rownames(chlorophyte_busco_tree_meta_df) <- chlorophyte_busco_tree_meta_df$accessions
chlorophyte_busco_tree_meta_df <- chlorophyte_busco_tree_meta_df[chlorophyte_busco_tree_order,]
chlorophyte_busco_sampleid <- paste0(chlorophyte_busco_tree_meta_df$accessions, " " ,chlorophyte_busco_tree_meta_df$taxa)
chlorophyte_busco_sampleid_df <- data.frame("accession" = chlorophyte_busco_tree_meta_df$accessions, chlorophyte_busco_sampleid, busco_habitat = chlorophyte_busco_tree_meta_df$habitat)
chlorophyte_busco_tree$tip.label <- chlorophyte_busco_sampleid
chlorophyte_busco_tree_meta_df$sample_id <- chlorophyte_busco_sampleid
chlorophyte_busco_tree_meta_df <- chlorophyte_busco_tree_meta_df[,c(4,1,2,3)]
chlorophyte_busco_tree_metadata <- ggtree(chlorophyte_busco_tree, branch.length = "none") %<+% chlorophyte_busco_tree_meta_df
chlorophyte_busco_tree_ggtree_colouredtips <- chlorophyte_busco_tree_metadata+geom_tiplab(aes(color = habitat),align = TRUE, linetype = "dotted", linesize = 0.1, fontface = 3)+
  scale_color_manual(values = chlorophyte_tree_colours)+
  theme( legend.position = "bottom")+theme( legend.position = "none", coord_cartesian(clip = "off"))+xlim(NA, 25)+ylim(NA,40)

# bacteria qs50 gtdbtk 
qs50_bacteria_tree <- read.tree("data/trees/qs50_bacteria_gtdbtk_midpoint.txt")
taxdata <- read.csv("metamdbg_gtdbtkeuk_tax_df.csv")
taxdata <- taxdata[taxdata$genome %in% qs50_bacteria_tree$tip.label,]
checkm <- read.csv("drep_genomeinfo.csv")
checkm$genome <- substr(checkm$genome, 1, (nchar(checkm$genome)-3))
checkm <- checkm[checkm$genome %in% qs50_bacteria_tree$tip.label,]
rownames(checkm) <- checkm$genome
checkm_comp <- checkm[,c(1,2)]
checkm_cont <- checkm[,c(1,3)]
checkm_comp$genome <- NULL
checkm_cont$genome <- NULL

checkm_contam <- as.data.frame(checkm[,c(2)])
ptax <- taxdata[,c(1,3)]
rownames(ptax) <- ptax$genome
ptax$genome <- NULL

qs50_mapping <- read.csv("qs50_mapping_tax_df.csv")
qs50_mapping <- qs50_mapping[qs50_mapping$genome %in% qs50_bacteria_tree$tip.label,]
qs50_mapping <- qs50_mapping[,c(1,2,10)]
qs50_mapping_wide <- pivot_wider(qs50_mapping, names_from = sample_id, values_from = Phylum, values_fill = "absent") %>% as.data.frame()
rownames(qs50_mapping_wide) <- qs50_mapping_wide$genome
qs50_mapping_wide$genome <- NULL

p <- ggtree(qs50_bacteria_tree, layout="dendrogram", open.angle=5, size = 0.15) %<+% ptax
p <- gheatmap(p, ptax, 
              colnames = FALSE, legend_title="Phylum", width = 0.1, color = "grey25")+scale_fill_manual(values = phylum_colours)
p2 <- p + new_scale_fill()
p2 <- gheatmap(p2, checkm_comp,  colnames = FALSE, legend_title = "Completeness", offset = 0.2, width = 0.1, color = "grey25")+
  scale_fill_gradient(low = "white", high = "#065535")
p3 <- p2 + new_scale_fill()
p3 <- gheatmap(p3, checkm_cont,  colnames = FALSE, legend_title = "Contamination", offset = 0.35, width = 0.1, color = "grey25")+
  scale_fill_gradient(low = "white", high = "#990000")
p4 <- p3 + new_scale_fill()
p4 <- gheatmap(p4, qs50_mapping_wide, colnames = FALSE, offset = 0.5, color = "grey25", width = 1.5)+
  scale_fill_manual(values = phylum_colours, na.value = "white")+scale_color_manual(values = c("grey30"))

ggsave("figures/microbial_community/mag_occur.pdf", p4, dpi = "retina", height = 8, width = 14)
p5 <- p4 + theme(legend.position = "none")
p4_legend <- get_legend(p4)

# cyanobacteria with refs
cyano_reftree <- read.tree("lk2JS2-72QYSiS-N7APDJg_newick.txt")
cyano_samples <- cyano_reftree$tip.label
cyano_habitats <- read.csv("cyano_habs.txt", header = FALSE, sep = ",")
colnames(cyano_habitats) <- c("genome", "colour","habitat")
rownames(cyano_habitats) <- cyano_habitats$genome
cyano_habitats$genome <- NULL
cyano_habitats$colour <- NULL
cyano_habitats$habitat <- factor(cyano_habitats$habitat, levels = c("Aquatic", "Terrestrial", "Azolla (Water Fern)",
                                                                    "Diatom Endosymbiont","Rhodophyta","Fruit Tree",
                                                                    "Cycad","Bryophyte","Lichen", "Cyanolichen (this study)"))
cyano_tax <- read.csv("itol_family_gtdbtk.csv")
rownames(cyano_tax) <- cyano_tax$genome
cyano_tax2 <- data.frame("genome" = cyano_tax$genome,"Family" = cyano_tax$Family)
rownames(cyano_tax2) <- cyano_tax2$genome
cyano_tax2$genome <- NULL

cyano_meta <- cyano_meta[cyano_meta$genome %in% cyano_samples,]
rownames(cyano_meta) <- cyano_meta$genome
cyano_meta$status <- "No"
cyano_meta[cyano_meta$habitat == "Cyanolichen (this study)",]$status <- "Yes"
cyano_meta$status_size <- 0
cyano_meta[cyano_meta$status == "Yes",]$status_size <- 3

cyanobacteria_family_colours <- c("Leptolyngbyaceae" = "#7b6690",
                               "Nostocaceae" = "#779977",
                               "Microcystaceae_A" = "#a6bea6",
                               "Thermosynechococcaceae" = "#c27ba0",
                               "LV9" = "#dbe1db",
                               "Xenococcaceae" = "#741b47",
                               "Chroococcidiopsidaceae" = "#baa6ce",
                               "Chamaesiphonaceae" = "#e6e5e7")

cyanobacteria_habitat_colours <- c("Cycad" = "#407d54",
                                  "Azolla (Water Fern)" = "#0b5394",
                                  "Aquatic" = "#b0b0b0",
                                  "Terrestrial" = "#494949",
                                  "Diatom Endosymbiont" = "#6c96bc",
                                  "Bryophyte" = "#689873",
                                  "Lichen" = "#e8b443",
                                  "Fruit Tree" = "#2e4837",
                                  "Rhodophyta" = "#c4d5e6",
                                  "Cyanolichen (this study)" = "#f9d52b")
c <- ggtree(cyano_reftree, layout="fan", open.angle=5, size = 0.15) %<+% cyano_meta
c <- c +geom_tippoint(aes(color = status, size = status_size), pch = 8)+scale_colour_manual(values = c("white","#f9d52b"))+
  scale_size_continuous(range = c(0,2))
c <- gheatmap(c, cyano_tax2, 
              colnames = FALSE, legend_title="Family", width = 0.1, color = "grey25")+scale_fill_manual(values = cyanobacteria_family_colours)
c2 <- c + new_scale_fill()
c2 <- gheatmap(c2, cyano_habitats,  colnames = FALSE, legend_title = "Habitat Classification", offset = 0.13,width = 0.1, color = "grey25")+
  scale_fill_manual(values = cyanobacteria_habitat_colours)

c3 <- c2 + theme(legend.position = "none")
ggsave("figures/cyanobacteria/cyanoreferencetree.pdf", c2, dpi = "retina", height = 12, width = 10)

# fungi - busco phylogenomics
fungi_busco_tree <- ape::read.tree("fungi_busco_phylogenomics_midpointroot_cladogram.txt")
fungi_busco_tree_tips <- fungi_busco_tree$tip.label
#fungi_busco_tree$tip.label["glNeplaev11_merged.1.fa"] <- "glNepLaev11_merged.1.fa"
fungi_metadata <- read.csv("fungi_refs_meta_accessions.csv", sep = ",")
fungi_metadata <- fungi_metadata[fungi_metadata$accession %in% fungi_busco_tree$tip.label,] %>% distinct()
fungi_busco_sampleid <- paste0(fungi_metadata$accession, " " ,fungi_metadata$organism)
fungi_meta_habitat <- fungi_metadata[,c("accession","class","order","family","life")]
fungi_meta_habitat[fungi_meta_habitat$life %in% c("not lichen (me)"),]$life <- "not lichenized (me)"
fungi_meta_habitat[fungi_meta_habitat$life %notin% c("lichen","lichen (me)", "not lichenized (me)"),]$life <- "not lichenized"
fungi_meta_habitat$sampleid <- fungi_busco_sampleid
rownames(fungi_meta_habitat) <- fungi_meta_habitat$accession
fungi_meta_habitat <- fungi_meta_habitat[(fungi_busco_tree_tips),]
write.csv("fungi_meta_habitat.csv", row.names = FALSE)
fungi_busco_tree2 <- fungi_busco_tree
fungi_busco_tree2$tip.label <- fungi_meta_habitat$sampleid
rownames(fungi_meta_habitat) <- fungi_meta_habitat$sampleid

fungi_life <- fungi_meta_habitat[,c("sampleid","life")]
fungi_life <- fungi_life[fungi_busco_tree2$tip.label,]
fungi_tree_metadata <- ggtree(fungi_busco_tree2, branch.length =  "none") %<+% fungi_life

fungi_tree_colours <- c("lichen" = "#bf9000",
                        "lichen (me)" = "#f1c232",
                        "not lichenized" = "#2f1c16",
                        "not lichenized (me)" = "#9c4b27")

fungi_tree <- fungi_tree_metadata+geom_tiplab(aes(color = life),align = TRUE, linetype = "dotted", linesize = 0.1, fontface = 3, size = 12)+
  scale_color_manual(values = fungi_tree_colours)+
  theme( legend.position = "bottom")+theme( legend.position = "none", coord_cartesian(clip = "off"))+xlim(0, 100)+ylim(NA,100)

saveRDS(fungi_tree, "00-setup/RDS/fungi_busco_ggtree.RDS")

# fungi busco - my mags only
fungi_mag_tree <- read.tree("data/trees/fungi_mag_tree.txt")
fungi_mag_tree_sampleid <- fungi_mag_tree$tip.label
fungi_mag_metadata <- fungi_metadata[fungi_metadata$accession %in% fungi_mag_tree_sampleid,] %>% distinct()

rownames(fungi_mag_metadata) <- fungi_mag_metadata$accession
fungi_mag_metadata <- fungi_mag_metadata[(fungi_mag_tree_sampleid),]
fungi_mag_life <- fungi_mag_metadata[,c("accession","life")]
fungi_mag_life <- fungi_mag_life[fungi_mag_tree$tip.label,]
fungi_mag_life[fungi_mag_life$life == "not lichen (me)",]$life <- "not lichenized (me)"
fungi_tree_metadata_mag <- ggtree(fungi_mag_tree, branch.length =  "none") %<+% fungi_mag_life

fungi_mag_ggtree <- fungi_tree_metadata_mag+geom_tiplab(aes(color = life),align = TRUE, linetype = "dotted", linesize = 0.1, fontface = 3, size = 5)+
  scale_color_manual(values = fungi_tree_colours)+
  theme(legend.position = "bottom")+theme( legend.position = "none", coord_cartesian(clip = "off"))+xlim(0, 100)

saveRDS(fungi_mag_ggtree, "00-setup/RDS/fungi_mag_ggtree.RDS")
