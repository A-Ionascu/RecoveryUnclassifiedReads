
packages <- c("dplyr","readr","utils")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages]) }
invisible(lapply(packages, library, character.only = TRUE))

maindir <- (paste(getwd(),sep=""))
path_unclassified <- file.path(maindir,"unclassified")
setwd(path_unclassified)

table_interest <- readr::read_tsv("blast_unique_barcodes.csv", col_names = FALSE)
table_interest <- as.data.frame(table_interest)
colnames(table_interest) <- c("Query_name", "Subject_name", "%identity", "length",
                             "mismatch", "gap", "query_start", "query_end",
                             "subject_start", "subject_end", "E_value", "Bit_score")

# eliminate repetitions
repetitie <- which(duplicated(table_interest$Subject_name))
table_no_repetitions_interest <- table_interest[!(rownames(table_interest) %in% repetitie),]

barcodes <- unique(table_no_repetitions_interest$Query_name)
for(bc in barcodes){
  barcode_subtable <-  dplyr::filter(table_no_repetitions_interest, Query_name == bc)
  id_reads_barcode <- data.frame(barcode_subtable$Subject_name)
  write.table(id_reads_barcode, file = paste(bc,".lst", sep=""), 
              sep="\n", qmethod = "escape",
              row.names = FALSE, col.names = FALSE)
}

# End of script