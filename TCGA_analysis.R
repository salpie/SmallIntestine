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

# get the start and stop locations
chrom_start_stop = Merge5[,c(1:3)]

# add new column names to the chrom start and stop
chrom_start_stop[,unique(tcga[,1])] <- NA

for (name in unique(tcga[,1])) {
	sub = tcga[tcga[,1]==name,]
	for (i in 1:nrow(chrom_start_stop)) {
		CN_rows = sub[(sub$Chromosome == chrom_start_stop$Chromosome[i]) & (sub$Start <= chrom_start_stop$Start[i]) & (sub$End >= chrom_start_stop$End[i]) | (sub$Chromosome == chrom_start_stop$Chromosome[i]) & (sub$Start <= chrom_start_stop$Start[i]) & (sub$End >= chrom_start_stop$Start[i]) | (sub$Chromosome == chrom_start_stop$Chromosome[i]) & (sub$Start <= chrom_start_stop$End[i]) & (sub$End >= chrom_start_stop$End[i]) | (sub$Chromosome == chrom_start_stop$Chromosome[i]) & (sub$Start >= chrom_start_stop$Start[i]) & (sub$End <= chrom_start_stop$End[i]),]
		CN = CN_rows$Modal_Total_CN
		if (dim(table(CN)) ==1) {
			chrom_start_stop[i,name] <- CN[1]
		}
		if (dim(table(CN)) > 1) {
			weights=vector()
			for (r in 1:nrow(CN_rows)) {
				weights=c(weights, sum(CN_rows$Start[r]:CN_rows$End[r] %in% chrom_start_stop$Start[i]:chrom_start_stop$End[i])/500000)
			}
			chrom_start_stop[i,name] <- weighted.mean(CN, weights)
		}
	}
}



# load premade 1Mb TCGA samples and duodenals
load("/Users/nowins01/Documents/Duodenal_project/Analysis/perf_cn/Merged_TCGA_DUODENALS_perfect.RData")
Merge5$Chromosome.y <- NULL
Merge5$Start.y <- NULL
Merge5$End.y <- NULL

############################################
################# PIC SCORE ################
############################################

### convert tcga into ga
copy = Merge5[,c(4:ncol(Merge5))]
copy[copy >= 3] <- 3
copy[copy < 2] <- 1

## calculate diversity per bin - PIC score

## number of samples that are loss at each bin

get_diversity <- function(x) {
	(1- ((sum(x == 1)/length(x))^2 + (sum(x == 2)/length(x))^2 + (sum(x == 3)/length(x))^2))
}

PIC_SCORE = as.data.frame(apply(copy, 2, get_diversity))
names(PIC_SCORE) <- "PIC"

############################################
################ GD/ploidy #################
############################################

gd = as.data.frame(names(Merge5)[4:ncol(Merge5)])
names(gd)[1] <- "sample"
gd$GD = 0

for (name in names(Merge5)[4:ncol(Merge5)]) {
	gd$GD[gd$sample==name] <- sum(Merge5[,name])/length(Merge5[,name])
}
	
gd$type <- "cohort"
gd$type_RL <- "cohort"

gd$clinical














 
