
packages <- c("dplyr","readr","utils")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages]) }
invisible(lapply(packages, library, character.only = TRUE))

maindir <- (paste(getwd(),sep=""))
path_unclassified <- file.path(maindir,"unclassified")
setwd(path_unclassified)

table <- readr::read_tsv("blast_all_barcodes.csv", col_names = FALSE)
table <- as.data.frame(table)
colnames(table) <- c("Query_name", "Subject_name", "%identity", "length",
                     "mismatch", "gap", "query_start", "query_end",
                     "subject_start", "subject_end", "E_value", "Bit_score")

# eliminate repetitions
repetitie <- which(duplicated(table$Subject_name))
table_no_repetitions <- table[!(rownames(table) %in% repetitie),]

# prepare unique barcoded reads
unique_barcodes <- data.frame(table_no_repetitions$Subject_name)
colnames(unique_barcodes) <- c("read_id")
unique_barcodes$read_id <- paste0(">", unique_barcodes$read_id)

write.table(unique_barcodes, file = "unique_all_barcodes.txt", sep="\n", qmethod="escape",
            row.names = FALSE, col.names = FALSE)

# End of script