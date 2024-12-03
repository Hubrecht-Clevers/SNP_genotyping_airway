#####################################################################################################
#                                                                                                   #
#                            SNP profiling between DNA and scRNA samples                            #
#                                          Gijs van Son                                             #
#                                      December 3rd 2024                                            #
#                               g.j.f.vanson@prinsesmaximacentrum.nl                                #
#                                                                                                   #
#####################################################################################################

# set up the environment
library("VariantAnnotation")
library("RColorBrewer")
library("ggplot2")
library("ggpubr")
library("cowplot")

setwd("#############")
myColors <- c("blue","green","red","white")
names(myColors) <- c("0/0","0/1","1/1",".")

myValues <- c(-1,0,1,NA)
names(myValues) <- c("0/0","0/1","1/1",".")

matrixcols2 <-  colorRampPalette(c("#377eb8","#4daf4a","#e41a1c"))
matrixcols3 <-  colorRampPalette(c("#377eb8","#4daf4a","#e41a1c"))

######## functions ############
dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  return(d)
} 

get_info_from_file <- function(inputfile){
  # have a look at what the file looks like, searching for the starting rows off all variable fields
  # we currently only use the Median values, but the rest is stored for easy access. 
  results_fileinfo <- read.table(inputfile, sep = "]", header = F, blank.lines.skip = F)
  Median_row <- grep("DataType:,Median", results_fileinfo[,1])
  NetMFI_row <- grep("DataType:,Net MFI", results_fileinfo[,1])
  Count_row <- grep("DataType:,Count", results_fileinfo[,1])
  AvgNetMFI_row <- grep("DataType:,Avg Net MFI", results_fileinfo[,1]) 
  Units_row <- grep("DataType:,Units", results_fileinfo[,1])
  PerBeadCount_row <- grep("DataType:,Per Bead Count", results_fileinfo[,1])
  DilutionFactor_row <- grep("DataType:,Dilution Factor", results_fileinfo[,1])
  AnalysisTypes_row <- grep("DataType:,Analysis Types", results_fileinfo[,1])
  calinfo_row <- grep("CALInfo:", results_fileinfo[,1])
  # Reading the file into different tables with information on the run
  calibrator_info <- read.table(inputfile, skip = calinfo_row+1, nrows = 1, sep = ",", header = T, blank.lines.skip = F)
  calibrator_info <- calibrator_info[,c(1:17)]
  results_median <- read.table(inputfile, skip = Median_row, nrows = 5, sep = ",", header = T, blank.lines.skip = F)
  results_netMFI <-  read.table(inputfile, skip = NetMFI_row, nrows = 5, sep = ",", header = T, blank.lines.skip = F)
  results_counts <-  read.table(inputfile, skip = Count_row, nrows = 5, sep = ",", header = T, blank.lines.skip = F)
  results_AVG_netMFI <- read.table(inputfile, skip = AvgNetMFI_row, nrows = 5, sep = ",", header = T, blank.lines.skip = F)
  results_Dilution_factor <- read.table(inputfile, skip = DilutionFactor_row, nrows = 5, sep = ",", header = T, blank.lines.skip = F)
  results_Dilution_factor <- results_Dilution_factor[,c(1,2,3)]
  # Storing the different data frames in variables to use in the script.
  assign("calibrator_info", calibrator_info, envir = .GlobalEnv)
  assign("results_median", results_median, envir = .GlobalEnv)
  assign("results_netMFI", results_netMFI, envir = .GlobalEnv)
  assign("results_counts", results_counts, envir = .GlobalEnv)
  assign("results_AVG_netMFI", results_AVG_netMFI, envir = .GlobalEnv)
  assign("results_Dilution_factor", results_Dilution_factor, envir = .GlobalEnv)
  
}

genotype_classification <- function(sample_in, design,  callibrator_file_name, method = "euclidian"){
  # function to assign genotypes to a sample based on a callibrated LUMINEX assay. 
  genotype_frame <- data.frame(matrix(NA, nrow = length(sample_in$Sample), ncol = length(design$ID)))
  group_frame <- data.frame(matrix(NA, nrow = length(sample_in$Sample), ncol = length(design$ID)))
  rownames(genotype_frame) <- sample_in$Sample
  colnames(genotype_frame) <- design$ID
  rownames(group_frame) <- sample_in$Sample
  colnames(group_frame) <- design$ID
  # Calculate the distance to the calibrator genotypes in X-Y space for all the location (2 probes per location)
  for (location in design$ID){
    calibrator_file <- read.csv(callibrator_file_name, row.names=1)
    calibrator_file_temp <- calibrator_file[,c(design[which(design$ID == location), "PROBE_ALT"], design[which(design$ID == location), "PROBE_REF"])]
    sample_list_location <- sample_in[,c("Sample",design[which(design$ID == location), "PROBE_ALT"], design[which(design$ID == location), "PROBE_REF"])]
    rownames(sample_list_location) <- sample_list_location$Sample
    sample_list_location <- sample_list_location[,-1]
    # Do this for every sample
    for (sample  in rownames(sample_list_location)) {
      x <- sample
      temp1 <- rbind(calibrator_file_temp, sample_list_location[x,])
      temp1_dist <- as.matrix(dist(temp1[,c(1,2)], method = method))
      sample_loc <- sample_list_location[x,]
      group <- c(
        mean(unlist(lapply(which(colnames(temp1_dist) %in% 
                                   rownames(calibrator_file[which(calibrator_file[[grep(paste0("chr",design$CHR[which(design$ID == location)],"_"),
                                                                                        colnames(calibrator_file), value = T)]] == "0/0"),])), 
                           function(y){
                             cal <- calibrator_file_temp[y,]
                             dist2d(as.numeric(sample_loc),c(0,0),as.numeric(cal))}))),
        mean(unlist(lapply(which(colnames(temp1_dist) %in% 
                                   rownames(calibrator_file[which(calibrator_file[[grep(paste0("chr",design$CHR[which(design$ID == location)],"_"),
                                                                                        colnames(calibrator_file), value = T)]] == "0/1"),])), 
                           function(y){
                             cal <- calibrator_file_temp[y,]
                             dist2d(as.numeric(sample_loc),c(0,0),as.numeric(cal))}))),
        mean(unlist(lapply(which(colnames(temp1_dist) %in% 
                                   rownames(calibrator_file[which(calibrator_file[[grep(paste0("chr",design$CHR[which(design$ID == location)],"_"),
                                                                                        colnames(calibrator_file), value = T)]] == "1/1"),])), 
                           function(y){
                             cal <- calibrator_file_temp[y,]
                             dist2d(as.numeric(sample_loc),c(0,0),as.numeric(cal))})))
      )
      call <- order(group)[1]
      difference_calls <- group[order(group)[1]]/group[order(group)[2]]
      genotype <- c("0/0", "0/1", "1/1")[call]
      if(difference_calls > 0.95){
        genotype <- NA
        call <- NA}
      genotype_frame[sample, location] <- genotype
      group_frame[sample, location] <- call
    }}
  return(list(genotype_frame, group_frame))
}

####### Load data ###########
get_info_from_file(inputfile = "20240215_CleversXilis_SNP_13_Nella.csv")
design <- read.csv("design.txt")
####### classify genotypes ##########
new_sample_file_table <- genotype_classification(results_median, design = design, 
                                                 callibrator_file_name = "/Users/gijsvanson/OneDrive - Prinses Maxima Centrum/General/SNP_Profiling/my_app/Calibrator/2022_03_17_calibrator.csv", method = "euclidean" )

dev.off()
####### print genotypes in a heatmap #######
pdf(file= "nella_genotypes.pdf", width=30, height=15, pointsize=15)
gplots::heatmap.2(as.matrix(new_sample_file_table[[2]][grep("Nella", rownames(new_sample_file_table[[2]])),]), trace="none", scale="none", 
          main="Genotype Nella", margins=c(10, 15), col=matrixcols3, 
          cexCol=1.2, dendrogram="row", Colv=F, symkey=T, 
          sepwidth=c(0.02,0.02), sepcolor="lightgray", colsep=1:ncol(new_sample_file_table[[2]]), 
          rowsep=1:nrow(new_sample_file_table[[2]]), keysize=0.5)
dev.off()
##### write genotypes to file #########
write.table(new_sample_file_table[[1]], "Nella_genotypes.csv")
##### view genotype from a vcf-file (from scRNA seq data) ##########
genotype_LX469 <- system("bcftools query -r 3:42209771,5:140668624,6:87258847,7:139047751,8:51820490,9:113188426,10:7964835,11:74074281,12:67312686,13:49534116,14:55138318,15:51622465,16:4701044,17:3026098,18:62574499,19:10489766,21:42904255 -f '%CHROM %POS %AF %GT\n' run.vcf")
##### make a new dataframe with all genotypes to perform clustering
results_sc_snp <- data.frame("chromosome" = c(1:21), "0.bam_LX469" = c(NA, NA, 2, 1, NA, 1, NA, 
                                                                       2, NA, 2, NA, 2, NA, 1, 2,
                                                                       NA, NA, NA, NA, 2,2),
                             "1.bam_LX469" = c(NA, NA, NA, 1, 1, 1, NA, 1, NA , 1, NA, NA, NA, 1, NA, NA, NA, NA, NA, 1, NA),
                             "2.bam_LX469" = c(NA, NA, NA, 2, NA, 2, NA, NA, NA, 2, NA, 1, 1, 2, NA, NA, NA, NA, NA, 1, 1),
                             "Nella1" = c(NA, NA, 1,2,1,2,1,1,2,2,1,1,1,2,0,0,1,2,1,NA,1),
                             "Nella2" = c(NA, NA, 2,1,1,1,1,1,2,1,1,1,1,1,1,0,1,1,2,NA,0),
                             "Nella3" = c(NA, NA, 1,2,0,1,2,2,1,2,2,2,1,1,2,1,1,0,0,NA,2),
                             "scale" = rep(1, 21))
pheatmap(results_sc_snp[,-1], cluster_cols = T, cluster_rows = T)
