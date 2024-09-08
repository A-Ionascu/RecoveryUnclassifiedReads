
# Load required libries 

#List packages to check installation
packages <- c("dplyr","readr","utils","ggplot2","forcats","nortest","tseries",
              "car","gridExtra","grid","vioplot","patchwork","rstatix","utils",
              "ggthemes","data.table","FactoMineR","factoextra","tidyverse",
              "corrplot")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages]) }
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


###############################################################################

setwd("S:/doc/etapa3_DNAseq/data_16S_sequencing/experimente_recovery_unclassified")

mainDir <- getwd()

folders <- list.dirs(recursive = FALSE)
folders <- gsub("./","",folders)
folders <- folders[-length(folders)]



for(dir in folders){
  
  # Enter directory
  path <- file.path(mainDir,dir,"epi2me","wf-16s_epi2me","output")
  setwd(path)
  
  
  
  ### Create directories for all levels of data
  ifelse(!dir.exists(file.path(path,"data_analysis")), 
         dir.create(file.path(path,"data_analysis")), FALSE)
  setwd(file.path(path,"data_analysis"))
  
  ifelse(!dir.exists(file.path(path,"data_analysis","genus")), 
         dir.create(file.path(path,"data_analysis","genus")), FALSE)
  ifelse(!dir.exists(file.path(path,"data_analysis","family")), 
         dir.create(file.path(path,"data_analysis","family")), FALSE)
  ifelse(!dir.exists(file.path(path,"data_analysis","order")), 
         dir.create(file.path(path,"data_analysis","order")), FALSE)
  ifelse(!dir.exists(file.path(path,"data_analysis","class")), 
         dir.create(file.path(path,"data_analysis","class")), FALSE)
  ifelse(!dir.exists(file.path(path,"data_analysis","phylum")), 
         dir.create(file.path(path,"data_analysis","phylum")), FALSE)
  

  ### Prepare and edit raw tables
  
  # GENUS
  #############################################################################
  
  setwd(file.path(path))
  
  data_genus <- read.csv("wf-16s-counts-genus.csv")
  data_genus <- data.frame(data_genus)
  
  setwd(file.path(path,"data_analysis","genus"))
  
  if(dir == "exp1_Abaumannii"){
    colnames(data_genus) <- c("genus","all_barcodes",
                              "Or_AB_51h","m057_AB_51h","m14_AB_51h",
                              "Or_AB+COL_69h","Or_AB+MEM_69h","Or_AB_69h",
                              "m057_AB_69h","m057_AB+COL_69h","m057_AB+MEM_69h",
                              "m14_AB_69h","m14_AB+COL_69h","m14_AB+MEM_69h",
                              "Or_control","m057_control","m14_control",
                              "total","superkingdom","kingdom","plylum","class","order","family","tax") }
  if(dir == "exp2_antibiotice"){
    colnames(data_genus) <- c("genus","all_barcodes",
                              "Or_COL_4h","Or_MEM_4h","Or_COL+MEM_4h",
                              "m057_COL_4h","m057_MEM_4h","m057_COL+MEM_4h",
                              "m14_COL_4h","m14_MEM_4h","m14_COL+MEM_4h",
                              "total","superkingdom","kingdom","plylum","class","order","family","tax") }
  if(dir == "exp3_Ekobei"){
    colnames(data_genus) <- c("genus","all_barcodes",
                              "Or_EK_51h","Or_EK_69h","Or_EK+COL_69h","Or_EK+MEM_69h",
                              "m14_EK_51h","m14_EK+COL_69h","m14_EK+MEM_69h",
                              "total","superkingdom","kingdom","plylum","class","order","family","tax") }
  if(dir == "exp4_Kpneumoniae"){
    colnames(data_genus) <- c("genus","all_barcodes",
                              "Or_KP_51h","Or_KP_69h","Or_KP+COL_69h","Or_KP+MEM_69h",
                              "m14_KP_51h","m14_KP_69h","m14_KP+COL_69h","m14_KP+MEM_69h",
                              "total","superkingdom","kingdom","plylum","class","order","family","tax") }

  
  data_genus <- data_genus %>%
    filter(genus != "Unknown") %>%
    filter(total > 100)
  rownames(data_genus) <- data_genus[,1]
  
  if(dir == "exp1_Abaumannii"){ data_genus_selection <- data_genus[,3:17] }
  if(dir == "exp2_antibiotice"){ data_genus_selection <- data_genus[,3:11] }
  if(dir == "exp3_Ekobei"){ data_genus_selection <- data_genus[,3:9] }
  if(dir == "exp4_Kpneumoniae"){ data_genus_selection <- data_genus[,3:10] }
  
  write.csv(x=data_genus, file="wf-16s-counts-genus_noUnK.csv")
  write.csv(x=data_genus_selection, file="wf-16s-counts-genus_selection.csv")
  #############################################################################
  
  
  # FAMILY
  #############################################################################
  
  setwd(file.path(path))
  
  data_family <- read.csv("wf-16s-counts-family.csv")
  data_family <- data.frame(data_family)
  
  setwd(file.path(path,"data_analysis","family"))
  
  if(dir == "exp1_Abaumannii"){
    colnames(data_family) <- c("family","all_barcodes",
                              "Or_AB_51h","m057_AB_51h","m14_AB_51h",
                              "Or_AB+COL_69h","Or_AB+MEM_69h","Or_AB_69h",
                              "m057_AB_69h","m057_AB+COL_69h","m057_AB+MEM_69h",
                              "m14_AB_69h","m14_AB+COL_69h","m14_AB+MEM_69h",
                              "Or_control","m057_control","m14_control",
                              "total","superkingdom","kingdom","plylum","class","order","tax") }
  if(dir == "exp2_antibiotice"){
    colnames(data_family) <- c("family","all_barcodes",
                              "Or_COL_4h","Or_MEM_4h","Or_COL+MEM_4h",
                              "m057_COL_4h","m057_MEM_4h","m057_COL+MEM_4h",
                              "m14_COL_4h","m14_MEM_4h","m14_COL+MEM_4h",
                              "total","superkingdom","kingdom","plylum","class","order","tax") }
  if(dir == "exp3_Ekobei"){
    colnames(data_family) <- c("family","all_barcodes",
                              "Or_EK_51h","Or_EK_69h","Or_EK+COL_69h","Or_EK+MEM_69h",
                              "m14_EK_51h","m14_EK+COL_69h","m14_EK+MEM_69h",
                              "total","superkingdom","kingdom","plylum","class","order","tax") }
  if(dir == "exp4_Kpneumoniae"){
    colnames(data_family) <- c("family","all_barcodes",
                              "Or_KP_51h","Or_KP_69h","Or_KP+COL_69h","Or_KP+MEM_69h",
                              "m14_KP_51h","m14_KP_69h","m14_KP+COL_69h","m14_KP+MEM_69h",
                              "total","superkingdom","kingdom","plylum","class","order","tax") }
  
  
  data_family <- data_family %>%
    filter(family != "Unknown") %>%
    filter(total > 100)
  rownames(data_family) <- data_family[,1]
  
  if(dir == "exp1_Abaumannii"){ data_family_selection <- data_family[,3:17] }
  if(dir == "exp2_antibiotice"){ data_family_selection <- data_family[,3:11] }
  if(dir == "exp3_Ekobei"){ data_family_selection <- data_family[,3:9] }
  if(dir == "exp4_Kpneumoniae"){ data_family_selection <- data_family[,3:10] }
  
  write.csv(x=data_family, file="wf-16s-counts-family_noUnK.csv")
  write.csv(x=data_family_selection, file="wf-16s-counts-family_selection.csv")
  #############################################################################
  
  
  # ORDER
  #############################################################################
  
  setwd(file.path(path))
  
  data_order <- read.csv("wf-16s-counts-order.csv")
  data_order <- data.frame(data_order)
  
  setwd(file.path(path,"data_analysis","order"))
  
  if(dir == "exp1_Abaumannii"){
    colnames(data_order) <- c("order","all_barcodes",
                               "Or_AB_51h","m057_AB_51h","m14_AB_51h",
                               "Or_AB+COL_69h","Or_AB+MEM_69h","Or_AB_69h",
                               "m057_AB_69h","m057_AB+COL_69h","m057_AB+MEM_69h",
                               "m14_AB_69h","m14_AB+COL_69h","m14_AB+MEM_69h",
                               "Or_control","m057_control","m14_control",
                               "total","superkingdom","kingdom","plylum","class","tax") }
  if(dir == "exp2_antibiotice"){
    colnames(data_order) <- c("order","all_barcodes",
                               "Or_COL_4h","Or_MEM_4h","Or_COL+MEM_4h",
                               "m057_COL_4h","m057_MEM_4h","m057_COL+MEM_4h",
                               "m14_COL_4h","m14_MEM_4h","m14_COL+MEM_4h",
                               "total","superkingdom","kingdom","plylum","class","tax") }
  if(dir == "exp3_Ekobei"){
    colnames(data_order) <- c("order","all_barcodes",
                               "Or_EK_51h","Or_EK_69h","Or_EK+COL_69h","Or_EK+MEM_69h",
                               "m14_EK_51h","m14_EK+COL_69h","m14_EK+MEM_69h",
                               "total","superkingdom","kingdom","plylum","class","tax") }
  if(dir == "exp4_Kpneumoniae"){
    colnames(data_order) <- c("order","all_barcodes",
                               "Or_KP_51h","Or_KP_69h","Or_KP+COL_69h","Or_KP+MEM_69h",
                               "m14_KP_51h","m14_KP_69h","m14_KP+COL_69h","m14_KP+MEM_69h",
                               "total","superkingdom","kingdom","plylum","class","tax") }
  
  
  data_order <- data_order %>%
    filter(order != "Unknown") %>%
    filter(total > 100)
  rownames(data_order) <- data_order[,1]
  
  if(dir == "exp1_Abaumannii"){ data_order_selection <- data_order[,3:17] }
  if(dir == "exp2_antibiotice"){ data_order_selection <- data_order[,3:11] }
  if(dir == "exp3_Ekobei"){ data_order_selection <- data_order[,3:9] }
  if(dir == "exp4_Kpneumoniae"){ data_order_selection <- data_order[,3:10] }
  
  write.csv(x=data_order, file="wf-16s-counts-order_noUnK.csv")
  write.csv(x=data_order_selection, file="wf-16s-counts-order_selection.csv")
  ############################################################################
  
  
  
  # CLASS
  #############################################################################
  
  setwd(file.path(path))
  
  data_class <- read.csv("wf-16s-counts-class.csv")
  data_class <- data.frame(data_class)
  
  setwd(file.path(path,"data_analysis","class"))
  
  if(dir == "exp1_Abaumannii"){
    colnames(data_class) <- c("class","all_barcodes",
                              "Or_AB_51h","m057_AB_51h","m14_AB_51h",
                              "Or_AB+COL_69h","Or_AB+MEM_69h","Or_AB_69h",
                              "m057_AB_69h","m057_AB+COL_69h","m057_AB+MEM_69h",
                              "m14_AB_69h","m14_AB+COL_69h","m14_AB+MEM_69h",
                              "Or_control","m057_control","m14_control",
                              "total","superkingdom","kingdom","plylum","tax") }
  if(dir == "exp2_antibiotice"){
    colnames(data_class) <- c("class","all_barcodes",
                              "Or_COL_4h","Or_MEM_4h","Or_COL+MEM_4h",
                              "m057_COL_4h","m057_MEM_4h","m057_COL+MEM_4h",
                              "m14_COL_4h","m14_MEM_4h","m14_COL+MEM_4h",
                              "total","superkingdom","kingdom","plylum","tax") }
  if(dir == "exp3_Ekobei"){
    colnames(data_class) <- c("class","all_barcodes",
                              "Or_EK_51h","Or_EK_69h","Or_EK+COL_69h","Or_EK+MEM_69h",
                              "m14_EK_51h","m14_EK+COL_69h","m14_EK+MEM_69h",
                              "total","superkingdom","kingdom","plylum","tax") }
  if(dir == "exp4_Kpneumoniae"){
    colnames(data_class) <- c("class","all_barcodes",
                              "Or_KP_51h","Or_KP_69h","Or_KP+COL_69h","Or_KP+MEM_69h",
                              "m14_KP_51h","m14_KP_69h","m14_KP+COL_69h","m14_KP+MEM_69h",
                              "total","superkingdom","kingdom","plylum","tax") }
  
  
  data_class <- data_class %>%
    filter(class != "Unknown") %>%
    filter(total > 100)
  rownames(data_class) <- data_class[,1]
  
  if(dir == "exp1_Abaumannii"){ data_class_selection <- data_class[,3:17] }
  if(dir == "exp2_antibiotice"){ data_class_selection <- data_class[,3:11] }
  if(dir == "exp3_Ekobei"){ data_class_selection <- data_class[,3:9] }
  if(dir == "exp4_Kpneumoniae"){ data_class_selection <- data_class[,3:10] }
  
  write.csv(x=data_class, file="wf-16s-counts-class_noUnK.csv")
  write.csv(x=data_class_selection, file="wf-16s-counts-class_selection.csv")
  ############################################################################
  
  
  
  # PHYLUM
  #############################################################################
  
  setwd(file.path(path))
  
  data_phylum <- read.csv("wf-16s-counts-phylum.csv")
  data_phylum <- data.frame(data_phylum)
  
  setwd(file.path(path,"data_analysis","phylum"))
  
  if(dir == "exp1_Abaumannii"){
    colnames(data_phylum) <- c("phylum","all_barcodes",
                              "Or_AB_51h","m057_AB_51h","m14_AB_51h",
                              "Or_AB+COL_69h","Or_AB+MEM_69h","Or_AB_69h",
                              "m057_AB_69h","m057_AB+COL_69h","m057_AB+MEM_69h",
                              "m14_AB_69h","m14_AB+COL_69h","m14_AB+MEM_69h",
                              "Or_control","m057_control","m14_control",
                              "total","superkingdom","kingdom","tax") }
  if(dir == "exp2_antibiotice"){
    colnames(data_phylum) <- c("phylum","all_barcodes",
                              "Or_COL_4h","Or_MEM_4h","Or_COL+MEM_4h",
                              "m057_COL_4h","m057_MEM_4h","m057_COL+MEM_4h",
                              "m14_COL_4h","m14_MEM_4h","m14_COL+MEM_4h",
                              "total","superkingdom","kingdom","tax") }
  if(dir == "exp3_Ekobei"){
    colnames(data_phylum) <- c("phylum","all_barcodes",
                              "Or_EK_51h","Or_EK_69h","Or_EK+COL_69h","Or_EK+MEM_69h",
                              "m14_EK_51h","m14_EK+COL_69h","m14_EK+MEM_69h",
                              "total","superkingdom","kingdom","tax") }
  if(dir == "exp4_Kpneumoniae"){
    colnames(data_phylum) <- c("phylum","all_barcodes",
                              "Or_KP_51h","Or_KP_69h","Or_KP+COL_69h","Or_KP+MEM_69h",
                              "m14_KP_51h","m14_KP_69h","m14_KP+COL_69h","m14_KP+MEM_69h",
                              "total","superkingdom","kingdom","tax") }
  
  
  data_phylum <- data_phylum %>%
    filter(phylum != "Unknown") %>%
    filter(total > 100)
  rownames(data_phylum) <- data_phylum[,1]
  
  if(dir == "exp1_Abaumannii"){ data_phylum_selection <- data_phylum[,3:17] }
  if(dir == "exp2_antibiotice"){ data_phylum_selection <- data_phylum[,3:11] }
  if(dir == "exp3_Ekobei"){ data_phylum_selection <- data_phylum[,3:9] }
  if(dir == "exp4_Kpneumoniae"){ data_phylum_selection <- data_phylum[,3:10] }
  
  write.csv(x=data_phylum, file="wf-16s-counts-phylum_noUnK.csv")
  write.csv(x=data_phylum_selection, file="wf-16s-counts-phylum_selection.csv")
  ############################################################################
  
  setwd(mainDir)
}



###############################################################################
### PCA for all types of data
###############################################################################


for(dir in folders){
  
  # Enter directory
  path <- file.path(mainDir,dir,"epi2me","wf-16s_epi2me","output")
  setwd(file.path(path,"data_analysis"))
  folders2 <- list.dirs(recursive = FALSE)
  folders2 <- gsub("./","",folders2)

for(dir2 in folders2){
  setwd(file.path(path,"data_analysis",dir2))
  files <- list.files()
  table_specific <- read.csv(files[2])
  table_specific <- data.frame(table_specific)
  rownames(table_specific) <- table_specific[,1]
  table_specific <- table_specific[,-1]
  
  
  PCA_analysis <- PCA(table_specific, graph=TRUE)
  ggsave("generic_PCA.jpg", plot=last_plot() )
  
  if(nrow(table_specific) <= 2){
    next
  }
  
  
  #########
  # Default PCA graph
  fviz_pca_ind(PCA_analysis, 
               geom.ind = "text", # show points only (nbut not "text")
               col.ind = rownames(table_specific), # color by groups
               palette = rainbow(nrow(table_specific)),
               #addEllipses = TRUE, # Concentration ellipses
               legend.title = "Reads for:",
               labelsize=6) +
    theme_bw() +
    theme(
      plot.title = element_text(size=20, face="bold"),
      axis.text.x = element_text(size=20, face="bold", angle=0, hjust=1, vjust=0.5),
      axis.title.x = element_text(size=20, face="bold"),
      axis.title.y = element_text(size=20, face="bold"),
      axis.text.y = element_text(size=20, face="bold"),
      legend.title = element_text(size=15),
      legend.text = element_text(size=15) ) +
    ggtitle(paste("PCA individuals"))
  ggsave("PCA_a_individuals.jpg", plot=last_plot() )
  
  
  ######################
  # Bacterial genera contributions in PC1
  fviz_contrib(PCA_analysis, choice="ind", sort.val = "desc", axes=1,
               color = "black", addlabels=TRUE) +
    geom_text(label = round(PCA_analysis$ind$contrib[,1],2), vjust=-0.5, hjust=0.4, size = 8) +
    ylim(0, max(PCA_analysis$ind$contrib[,1])+5) +
    theme_bw() +
    theme(
      plot.title = element_text(size=20, face="bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=20, face="bold", angle=90, hjust=1, vjust=0.5),
      axis.title.y = element_text(size=20, face="bold"),
      axis.text.y = element_text(size=20, face="bold") ) +
    ggtitle(paste("Contribution of individuals to PC1"))
  ggsave("PCA_b_contribution_individuals_PC1.jpg", plot=last_plot() )
  
  ######################
  # Bacterial genera contributions in PC2
  fviz_contrib(PCA_analysis, choice="ind", sort.val = "desc", axes=2,
               color = "black", addlabels=TRUE) +
    geom_text(label = round(PCA_analysis$ind$contrib[,2],2), vjust=-0.5, hjust=0.4, size = 8) +
    ylim(0, max(PCA_analysis$ind$contrib[,2])+5) +
    theme_bw() +
    theme(
      plot.title = element_text(size=20, face="bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=20, face="bold", angle=90, hjust=1, vjust=0.5),
      axis.title.y = element_text(size=20, face="bold"),
      axis.text.y = element_text(size=20, face="bold") ) +
    ggtitle(paste("Contribution of individuals to PC2"))
  ggsave("PCA_b_contribution_individuals_PC2.jpg", plot=last_plot() )
  
  ######################
  # D.mel lines contributions in PC1
  fviz_contrib(PCA_analysis, choice="var", sort.val = "none", axes=1,
               color = "black", addlabels=TRUE) +
    geom_text(label = round(PCA_analysis$var$contrib[,1],2), vjust=-0.5, hjust=0.4, size = 8) +
    ylim(0, max(PCA_analysis$var$contrib[,1])+5) +
    theme_bw() +
    theme(
      plot.title = element_text(size=20, face="bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=20, face="bold", angle=90, hjust=1, vjust=0.5),
      axis.title.y = element_text(size=20, face="bold"),
      axis.text.y = element_text(size=20, face="bold") ) +
    ggtitle(paste("Contribution of variables to PC1"))
  ggsave("PCA_c_contribution_variables_PC1.jpg", plot=last_plot() )
  
  ######################
  # D.mel lines contributions in PC2
  fviz_contrib(PCA_analysis, choice="var", sort.val = "none", axes=2,
               color = "black", addlabels=TRUE) +
    geom_text(label = round(PCA_analysis$var$contrib[,2],2), vjust=-0.5, hjust=0.4, size = 8) +
    ylim(0, max(PCA_analysis$var$contrib[,2])+5) +
    theme_bw() +
    theme(
      plot.title = element_text(size=20, face="bold"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=20, face="bold", angle=90, hjust=1, vjust=0.5),
      axis.title.y = element_text(size=20, face="bold"),
      axis.text.y = element_text(size=20, face="bold") ) +
    ggtitle(paste("Contribution of variables to PC2"))
  ggsave("PCA_c_contribution_variables_PC2.jpg", plot=last_plot() )
  
  
  ######################
  # PCA dimensions variance
  fviz_eig(PCA_analysis, choice="variance", addlabels = FALSE) +
    geom_text(label = round(get_eig(PCA_analysis)[,2],2), vjust=-0.5, hjust=0.4, size = 8) +
    ylim(0, 100) +
    theme_bw() +
    theme(
      plot.title = element_text(size=20, face="bold"),
      axis.title.x = element_text(size=20, face="bold"),
      axis.text.x = element_text(size=20, face="bold"),
      axis.title.y = element_text(size=20, face="bold"),
      axis.text.y = element_text(size=20, face="bold") ) +
    ggtitle(paste("Explained percentage of variance in dimensions"))
  ggsave("PCA_d_variance.jpg", plot=last_plot() )
  
  
  ######################
  # PCA dimensions variance
  fviz_eig(PCA_analysis, choice="eigenvalue", addlabels = FALSE) +
    geom_text(label = round(get_eig(PCA_analysis)[,1],2), vjust=-0.5, hjust=0.4, size = 8) +
    ylim(0, max(get_eig(PCA_analysis)[,1])+1) +
    theme_bw() +
    theme(
      plot.title = element_text(size=20, face="bold"),
      axis.title.x = element_text(size=20, face="bold"),
      axis.text.x = element_text(size=20, face="bold"),
      axis.title.y = element_text(size=20, face="bold"),
      axis.text.y = element_text(size=20, face="bold") ) +
    ggtitle(paste("Eigenvalues in dimensions"))
  ggsave("PCA_e_eigenvalue.jpg", plot=last_plot() )
  
  
  ######################
  # PCA bi-plot variance
  fviz_pca_biplot(PCA_analysis, 
                  palette = rep(c("black"),9),
                  #col.ind = rownames(PCA_analysis),
                  pointsize=2, labelsize=6,
                  repel = TRUE,
                  col.var = "blue", arrowsize = 1) +
    theme_bw() +
    theme(
      plot.title = element_text(size=20, face="bold"),
      axis.text.x = element_text(size=20, face="bold"),
      axis.title.x = element_text(size=20, face="bold"),
      axis.title.y = element_text(size=20, face="bold"),
      axis.text.y = element_text(size=20, face="bold") ) +
    ggtitle(paste("PCA biplot"))
  ggsave("PCA_f_biplot_1.jpg", plot=last_plot() )
  
  
  fviz_pca_biplot(PCA_analysis, 
                  #col.ind = colnames(all_reads), 
                  #palette = rainbow(15), 
                  #habillage=colnames(all_reads),
                  #addEllipses = TRUE, label = "var", ellipse.type="confidence",
                  pointsize = 2, labelsize = 6,
                  col.var = "grey40", 
                  alpha.var = 0.75,
                  repel = TRUE,
                  legend.title = "Strain") +
    theme_bw() +
    theme(
      #legend.position = "none",
      plot.title = element_text(size=20, face="bold"),
      axis.text.x = element_text(size=20, face="bold"),
      axis.title.x = element_text(size=20, face="bold"),
      axis.title.y = element_text(size=20, face="bold"),
      axis.text.y = element_text(size=20, face="bold"),
      legend.title = element_text(size=20, face="bold"), 
      legend.text = element_text(size=20),
      legend.spacing.y = unit(0.5, "cm")
    ) +
    guides(colour = guide_legend(byrow = TRUE)) +
    ggtitle(paste("PCA applied on all read counts.")) 
  ggsave("PCA_f_biplot_2.jpg", plot=last_plot() )
  
  

  fviz_pca_biplot(PCA_analysis, 
                  #col.ind = colnames(all_reads), 
                  #palette = rainbow(15), 
                  #habillage=colnames(all_reads),
                  #addEllipses = TRUE, label = "var", ellipse.type="confidence",
                  pointsize = 2, labelsize = 6,
                  col.var = "blue", 
                  alpha.var = 0.75,
                  repel = TRUE,
                  legend.title = "Strain") +
    theme_bw() +
    theme(
      #legend.position = "none",
      plot.title = element_text(size=20, face="bold"),
      axis.text.x = element_text(size=20, face="bold"),
      axis.title.x = element_text(size=20, face="bold"),
      axis.title.y = element_text(size=20, face="bold"),
      axis.text.y = element_text(size=20, face="bold"),
      legend.title = element_text(size=20, face="bold"), 
      legend.text = element_text(size=20),
      legend.spacing.y = unit(0.5, "cm")
    ) +
    guides(colour = guide_legend(byrow = TRUE)) +
    ggtitle(paste("PCA applied on all read counts.")) 
  ggsave("PCA_f_biplot_3.jpg", plot=last_plot() )
  

  
  
  
  
}
  }













