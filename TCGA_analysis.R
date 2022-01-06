# uses TCGA ABSOLUTE copy number calls to get 1Mb sized copy number data for samples
# calculates ploidy and PGA
# plots heatmaps of duodenal lpWGS absolute copy number calls
# plots heatmaps of TCGA COAT, PRAAD, OES, READ, STOM and duodenals
# compares ploidy and PGA for all sample groups
# calcultes PIC score of population-wise diversity across the genome and compares how heterogeneous copy number calls are
# compares number of samples that are GD

#Load Libraries
library(data.table)

#############################################################
######## get 1Mb copy number calls from TCGA dataset ########
#############################################################

tcga = read.table("/Users/nowins01/Documents/ItsaCIN/publicly_avail/tcga/absolute_segs/TCGA_mastercalls.abs_segtabs.fixed.txt", header=T, sep="\t")

data<-as.data.frame(fread("/Users/nowins01/Documents/ItsaCIN/publicly_avail/tcga/clincal/clinical.cases_selection.2021-02-17/clinical.tsv"))

# all the colorectal data
coad_data = data[data$project_id=="TCGA-COAD",]

# just the left side
LSCC = coad_data[coad_data$tissue_or_organ_of_origin == "Cecum" | coad_data$tissue_or_organ_of_origin == "Ascending colon" | coad_data$tissue_or_organ_of_origin == "Hepatic flexure of colon",]

# just the right side
RSCC = coad_data[coad_data$tissue_or_organ_of_origin == "Descending colon" | coad_data$tissue_or_organ_of_origin == "Sigmoid colon" | coad_data$tissue_or_organ_of_origin == "Splenic flexure of colon" | coad_data$tissue_or_organ_of_origin == "Rectosigmoid junction",]

# get the right details for the COAD
tcga$TCGA.D5.6536Sample <- gsub("-01", "", tcga$TCGA.D5.6536Sample)

joined_df <- merge(coad_data, tcga, by.x = "case_submitter_id", 
             by.y = "TCGA.D5.6536Sample", all.x = TRUE, all.y = FALSE)


# load premade 1Mb TCGA samples and duodenals
load("/Users/nowins01/Documents/Duodenal_project/Analysis/perf_cn/Merged_TCGA_DUODENALS_perfect.RData")
Merge5$Chromosome.y <- NULL
Merge5$Start.y <- NULL
Merge5$End.y <- NULL

############################################
################# PIC SCORE ################
############################################



############################################
################ GD/ploidy #################
############################################

gd = as.data.frame(names(Merge5)[4:ncol(Merge5)])
names(gd)[1] <- "sample"
gd$GD = 0

for (name in names(Merge5)[4:ncol(Merge5)]) {
	gd$GD[gd$sample==name] <- sum(Merge5[,name])/length(Merge5[,name])
}
	















 
