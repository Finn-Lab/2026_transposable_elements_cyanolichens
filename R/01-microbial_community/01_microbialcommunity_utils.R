# functions
txt_from_folder <- function(path, ext = ".tsv"){
  files <- list.files(path)
  filepaths <- paste0(path,files)
  filelist <- vector("list",length(files))
  newcolnames <- c("genome","contig_length", "contig_count", "sample_mean_depth","sample_coverage","sample_expcoverage","sample_coeffvar")
  for (i in 1:length(files)){
    table <- read.table(filepaths[i], sep = "\t", header = TRUE)
    colnames(table) <- newcolnames
    filename <- strsplit(files[i], "_qs50_")[[1]][1]
    table$sample_id <- paste0(filename)
    filelist[[i]] <- table
  }
  table_final <- do.call(rbind, filelist)
  table_final
}


# colour palettes
phylum_colours <-  c("Acidobacteriota"= "#5095d4", "Actinomycetota" = "#a6c4df", 
                      "Armatimonadota" = "#dceaf7", "Ascomycota - Primary Symbiont" = "#d17850",
                      "Ascomycota - Additional Fungi" = "#f2c181","Bacteroidota" = "#943f6a",
                      "Chloroflexota" = "#cf9ab4", "Chlorophyta" = "#638666",
                      "Cyanobacteriota" = "#a2ceaa", "Myxococcota" = "#f2e2ea",
                      "Planctomycetota" = "grey18", "Pseudomonadota" = "grey50",
                      "Verrucomicrobiota" = "snow3")

cyanobacteria_bin_colours <- c("FACHB-36_1" = "#50336b","Kovacikia_1" = "#7b6690",
                               "Stenomitos_1"= "#a69ab5","Stenomitos_2" = "#d3ccda",
                               "Nostoc_1" = "#445544","Nostoc_2" = "#779977",
                               "Nostoc_3" = "#cccc88", "Nostoc_4" = "#50682c",
                               "Nostoc_5" = "#81964e", "Nostoc_6" = "#b4c572")
# metadata lists
lichen_species_fac <- c("Ricasolia virens", 
                        "Lobaria pulmonaria", 
                        "Pseudocyphellaria norvegica",
                        "Sticta sylvatica", 
                        "Nephroma laevigatum", 
                        "Peltigera collina",
                        "Peltigera horizontalis", 
                        "Peltigera membranacea",
                        "Peltigera praetextata" ,
                        "Peltigera hymenina",
                        "Leptogium burgessi")


lichen_species <- c("Ricasolia virens" = "glRicVire11", 
                    "Lobaria pulmonaria" = "glLobPulm2", 
                    "Pseudocyphellaria norvegica" = "glPseNorv1",
                    "Sticta sylvatica" = "glStiSylv4", 
                    "Nephroma laevigatum" = "glNepLaev11", 
                    "Peltigera collina" = "glPelColl1",
                    "Peltigera horizontalis" = "glPelHori1", 
                    "Peltigera membranacea" = "glPelMemb1",
                    "Peltigera praetextata" = "glPelPrae3" ,
                    "Peltigera hymenina" = "glPelHyme1",
                    "Leptogium burgessi" = "glLepBurg3")

mycobiont_mag <- c("glLepBurg3_bin.203.fa" = "Leptogium burgessi",
                   "glPseNorv1_bin.13.fa" = "Pseudocyphellaria norvegica",
                   "glStiSylv4_bin.60.fa" = "Sticta sylvatica",
                   "glRicVire11_merged.0.fa" = "Ricasolia virens",
                   "glLobPulm2_bin.43.fa" = "Lobaria pulmonaria",
                   "glNepLaev11_merged.0.fa" = "Nephroma laevigatum",
                   "glPelHyme1_bin.206.fa" = "Peltigera hymenina",
                   "glPelHori1_bin.14.fa" = "Peltigera horizontalis",
                   "glPelColl1_bin.25.fa" = "Peltigera collina",
                   "glPelMemb1_bin.35.fa" = "Peltigera membranacea",
                   "glPelPrae3_bin.50.fa" = "Peltigera praetextata")

phylum_taxonomiclevels <- c("Ascomycota - Primary Symbiont","Ascomycota - Additional Fungi",
"Chlorophyta","Cyanobacteriota",
"Acidobacteriota","Actinomycetota",
"Armatimonadota", "Bacteroidota",
"Chloroflexota","Myxococcota",
"Planctomycetota","Pseudomonadota",
"Verrucomicrobiota")

cyano_mag <- c("FACHB-36_1","Kovacikia_1","Stenomitos_1","Stenomitos_2",
               "Nostoc_1","Nostoc_2","Nostoc_3","Nostoc_4"
               ,"Nostoc_5","Nostoc_6")

# reorder to be bacterial (3 x 3 colour pairs), fungi (1 x 2 pair), algae(1 x 2 colour pair)         
all_mapping_weighteddepth_v3$phy <- factor(all_mapping_weighteddepth_v3$Phylum, 
                                              levels = c("Ascomycota - Primary Symbiont","Ascomycota - Additional Fungi",
                                                         "Chlorophyta","Cyanobacteriota",
                                                         "Acidobacteriota","Actinomycetota",
                                                         "Armatimonadota", "Bacteroidota",
                                                         "Chloroflexota","Myxococcota",
                                                         "Planctomycetota","Pseudomonadota",
                                                         "Verrucomicrobiota"))

