
# List packages to check installation
packages <- c("tibble","stringr")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
suppressWarnings(suppressMessages(invisible(lapply(packages, library, character.only = TRUE))))



mainDir <- (paste(getwd(),sep=""))

setwd(mainDir)
directories <- list.dirs(recursive=FALSE)

table <- data.frame(matrix(NA, nrow=13, ncol=1))
colnames(table) <- c("t")
rownames(table) <- c("reads","bases","mean_length","median_length","length_stdev",
                    "N50","mean_quality","median_quality","Q>5","Q>7","Q>10","Q>12","Q>15")

for(d in directories){
  dir <- gsub("./","",d)
  path <- file.path(mainDir,dir)
  setwd(path)

  nanoplot <- readr::read_tsv("NanoStats.csv", col_names = TRUE, show_col_types = FALSE)
  nanoplot <- as.data.frame(nanoplot)


  data <- data.frame(matrix(NA, nrow=13, ncol=1))
  colnames(data) <- c(substring(dir,23))
  rownames(data) <- c("reads","bases","mean_length","median_length","length_stdev",
                    "N50","mean_quality","median_quality","Q>5","Q>7","Q>10","Q>12","Q>15")

  data[1,1] <- as.numeric(nanoplot[1,2])
  data[2,1] <- as.numeric(nanoplot[2,2])
  data[3,1] <- as.numeric(nanoplot[3,2])
  data[4,1] <- as.numeric(nanoplot[4,2])
  data[5,1] <- as.numeric(nanoplot[5,2])
  data[6,1] <- as.numeric(nanoplot[6,2])
  data[7,1] <- as.numeric(nanoplot[7,2])
  data[8,1] <- as.numeric(nanoplot[8,2])

  data[9,1] <- stringr::str_extract(string = nanoplot[19,2],
                                   pattern = "(?<=\\().*(?=\\%)")
  data[10,1] <- stringr::str_extract(string = nanoplot[20,2],
                                  pattern = "(?<=\\().*(?=\\%)")
  data[11,1] <- stringr::str_extract(string = nanoplot[21,2],
                                  pattern = "(?<=\\().*(?=\\%)")
  data[12,1] <- stringr::str_extract(string = nanoplot[22,2],
                                  pattern = "(?<=\\().*(?=\\%)")
  data[13,1] <- stringr::str_extract(string = nanoplot[23,2],
                                  pattern = "(?<=\\().*(?=\\%)")

  table <- cbind(table, data)
  
  setwd(mainDir)
}

table <- table[-1]
write.csv(table, "nanoplot_stats.csv")



