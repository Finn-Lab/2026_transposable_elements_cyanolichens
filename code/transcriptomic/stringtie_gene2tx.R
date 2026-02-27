library(rtracklayer)

sp <- c(
  "glNepLaev11",
  "glPelHyme1"
)

sp_gtfs <- lapply(sp, function(x)
  import(paste0("./stringtie_combo_final_v2/", x, ".combined.gtf")))

sp_tab <- vector("list", length(sp))
for (i in seq_along(sp)) {
  sp_tab[[i]] <- data.frame(
    gene_id = sp_gtfs[[i]]$gene_id,
    tx_id = sp_gtfs[[i]]$transcript_id
  )
  sp_tab[[i]] <- sp_tab[[i]][!duplicated(sp_tab[[i]]), ]
}

for (i in seq_along(sp)) {
  readr::write_tsv(sp_tab[[i]],
    paste0("./stringtie_combo_final_v2/", sp[i], ".gene_to_tx.txt"),
    col_names = FALSE)
}
