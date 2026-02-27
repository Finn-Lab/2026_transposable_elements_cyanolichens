library(rtracklayer)

sp <- c(
  "gffcmp"
)

sp_gtfs <- lapply(sp, function(x)
  import(paste0("./combined_trebouxoid/", x, ".combined.gtf")))

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
    paste0("./combined_trebouxoid/", sp[i], ".gene_to_tx.txt"),
    col_names = FALSE)
}
