# uses TCGA ABSOLUTE copy number calls to get 1Mb sized copy number data for samples
# calculates ploidy and PGA
# plots heatmaps of duodenal lpWGS absolute copy number calls
# plots heatmaps of TCGA COAT, PRAAD, OES, READ, STOM and duodenals
# compares ploidy and PGA for all sample groups
# calcultes PIC score of population-wise diversity across the genome and compares how heterogeneous copy number calls are
# compares number of samples that are GD

#Load Libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(viridis)

# functions
'%!in%' <- function(x,y)!('%in%'(x,y))

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

for_analysis = unique(tcga$TCGA.D5.6536Sample)[unique(tcga$TCGA.D5.6536Sample) %in% names(Merge5)]
tcga = tcga[tcga$TCGA.D5.6536Sample %in% for_analysis,]

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

save(chrom_start_stop, file="chrom_start_stop.RData")

# load duodenals
load("/Users/nowins01/Documents/Duodenal_project/lpWGS/jupyter/Complete_Duodenal_Samples.RData")
comp$local = paste(comp$Chromosome, comp$Start, comp$End, sep=":")
chrom_start_stop$local = paste(chrom_start_stop$Chromosome, chrom_start_stop$Start, chrom_start_stop$End, sep=":")
Merge5 = merge(comp, chrom_start_stop, by="local")

Merge5$local <- NULL
Merge5$Start.y <-NULL
Merge5$End.y <-NULL
Merge5$Chromosome.y <-NULL

names(Merge5)[1] <-"Chromosome"
names(Merge5)[2] <- "Start"
names(Merge5)[3] <- "End"

# load premade 1Mb TCGA samples and duodenals
#load("/Users/nowins01/Documents/Duodenal_project/Analysis/perf_cn/Merged_TCGA_DUODENALS_perfect.RData")
#Merge5$Start.y <-NULL
#Merge5$End.y <-NULL
#Merge5$Chromosome.y <-NULL
#Merge5$Start.x <-NULL
#Merge5$End.x <-NULL



############################################
################ GD/ploidy #################
############################################

gd = as.data.frame(names(Merge5)[4:ncol(Merge5)])
names(gd)[1] <- "sample"
gd$GD = 0

New_Merge =Merge5
New_Merge = round(New_Merge)

for (name in names(Merge5)[4:ncol(Merge5)]) {
        gd$GD[gd$sample==name] <- sum(Merge5[,name])/length(Merge5[,name])
        New_Merge[,name] = (New_Merge[,name] - abs(round(gd$GD[gd$sample==name])- 2))
}




############################################
################# PIC SCORE ################
############################################

### convert tcga into ga
copy = New_Merge[,c(4:ncol(New_Merge))]
copy[copy >= 3] <- 3
copy[copy < 2] <- 1

## calculate diversity per bin - PIC score

## number of samples that are loss at each bin

get_diversity <- function(x) {
	(1- ((sum(x == 1)/length(x))^2 + (sum(x == 2)/length(x))^2 + (sum(x == 3)/length(x))^2))
}

PIC_SCORE = as.data.frame(apply(copy, 2, get_diversity))
names(PIC_SCORE) <- "PIC"

PIC_SCORE$sample <- rownames(PIC_SCORE)

##################
#################  Merge pic and ploidy etc

# merge PIC and ploidy and PGA
diversity_ploidy = merge(PIC_SCORE, gd, by="sample")

# and clinical
clin = data[,c(2,3)]
diversity_ploidy_clin = merge(diversity_ploidy, clin, by.x="sample", by.y="case_submitter_id")

# and duodenals
DUO_clin = diversity_ploidy[diversity_ploidy$sample %!in% diversity_ploidy_clin$sample,]
DUO_clin$project_id <- "Duodenum"

total = rbind(diversity_ploidy_clin, DUO_clin)

####### get PGA

pga = as.data.frame(names(Merge5)[4:ncol(Merge5)])
names(pga)[1] <- "sample"
pga$PGA = 0

for (name in names(Merge5)[4:ncol(Merge5)]) {
        pga$PGA[gd$sample==name] <- sum(Merge5[,name] != 2)/length(Merge5[,name])
}

# merge PIC and ploidy and PGA
diversity_ploidy_pga = merge(total, pga, by="sample")
diversity_ploidy_pga = unique(diversity_ploidy_pga)

##################### PLOT

ggplot(total, aes(x=project_id, y=PIC, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#9986A5", "#79402E", "#CCBA72", "#0F0D0E", "#D9D0D3", "#8D8680"))+
geom_boxplot(width=0.1, fill="white")+ylim(c(0,1))

pdf("Figure2.pdf")
ggplot(diversity_ploidy_pga, aes(x=project_id, y=GD, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#9986A5", "#79402E", "#CCBA72", "#0F0D0E", "#D9D0D3", "#8D8680"))+
geom_boxplot(width=0.1, fill="white")+ylim(c(0,8))
dev.off()

pdf("Figure3.pdf")
ggplot(diversity_ploidy_pga, aes(x=project_id, y=PGA, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#9986A5", "#79402E", "#CCBA72", "#0F0D0E", "#D9D0D3", "#8D8680"))+
geom_boxplot(width=0.1, fill="white")+ylim(c(0,1))
dev.off()


# Plot along genome
total %>%
  ggplot( aes(x=year, y=n, group=name, color=name)) +
    geom_line() +
    scale_color_viridis(discrete = TRUE) +
    ggtitle("Population-wde diversity") +
    ylab("Diversity along genome")
	
gd$type <- "cohort"
gd$type_RL <- "cohort"

gd$clinical

###################### 

# PIC score along genome
PIC_SCORE_bin = as.data.frame(apply(copy, 1, get_diversity))

names = total[,c(1,4)]
panc = names$sample[names$project_id == "TCGA-PAAD"]
panc_pic = as.data.frame(apply(copy[,names(copy) %in% panc], 1, get_diversity))
names(panc_pic)[1] <- "PIC_bin"
panc_pic$type <- "Pancreas"
panc_pic$num <- rownames(panc_pic)

coad = names$sample[names$project_id == "TCGA-COAD"]
coad_pic = as.data.frame(apply(copy[,names(copy) %in% panc], 1, get_diversity))
names(coad_pic)[1] <- "PIC_bin"
coad_pic$type <- "Colon"
coad_pic$num <- rownames(coad_pic)

oes = names$sample[names$project_id == "TCGA-ESCA"]
oes_pic = as.data.frame(apply(copy[,names(copy) %in% panc], 1, get_diversity))
names(oes_pic)[1] <- "PIC_bin"
oes_pic$type <- "Oesophagus"
oes_pic$num <- rownames(oes_pic)

stom = names$sample[names$project_id == "TCGA-STAD"]
stom_pic = as.data.frame(apply(copy[,names(copy) %in% panc], 1, get_diversity))
names(stom_pic)[1] <- "PIC_bin"
stom_pic$type <- "Stomach"
stom_pic$num <- rownames(stom_pic)

rec = names$sample[names$project_id == "TCGA-READ"]
rec_pic = as.data.frame(apply(copy[,names(copy) %in% panc], 1, get_diversity))
names(rec_pic)[1] <- "PIC_bin"
rec_pic$type <- "Rectum"
rec_pic$num <- rownames(rec_pic)

duo = names$sample[names$project_id == "Duodenum"]
duo_pic = as.data.frame(apply(copy[,names(copy) %in% panc], 1, get_diversity))
names(duo_pic)[1] <- "PIC_bin"
duo_pic$type <- "Duodenum"
duo_pic$num <- rownames(duo_pic)

total_pic_bin = rbind(oes_pic, stom_pic, panc_pic, duo_pic, coad_pic, rec_pic)

ggplot(total_pic_bin, aes(x = num, y = PIC_bin, colour = type)) +
  geom_line()+ylim=c(0,1)

ggplot(oes_pic, aes(x = num, y = PIC_bin, colour = type)) +
  geom_line()+ylim=c(0,1)

# Plot along genome
total_pic_bin %>%
  ggplot( aes(x=num, y=PIC_bin, group=type, color=type)) +
    geom_line() +
    scale_color_viridis(discrete = TRUE) +
    ggtitle("Population-wde diversity") +
    ylab("Diversity along genome")

ggplot(total_pic_bin, aes(fill=type, y=PIC_bin, x=num)) + 
    geom_bar(stat="identity", position="fill") +
    geom_line(data=total_pic_bin, aes(x=num,y=PIC_bin, group=type), color="black") +
    scale_fill_manual(values=c("#9986A5", "#79402E", "#CCBA72", "#0F0D0E", "#D9D0D3", "#8D8680")) +
    ggtitle("bins_PIC")


ggplot(total_pic_bin, aes(y=PIC_bin, x=num)) +
  geom_line(aes(linetype = type))


ggplot(total_pic_bin, aes(colour=type, y=PIC_bin, x=num, group=type)) + geom_line()   
    geom_bar(stat="identity", position="fill") +
    geom_line(data=total_pic_bin, aes(x=num,y=PIC_bin, group=type))

save(total_pic_bin, file="total_pic_bin.RData")
















for (name in unique(tcga[,1])) {
        sub = tcga[tcga[,1]==name,]
	print(name)
        for (i in 1:nrow(chrom_start_stop)) {
                CN_rows = sub[(sub$Chromosome == chrom_start_stop$Chromosome[i]) & (sub$Start <= chrom_start_stop$Start[i]) & (sub$End >= chrom_start_stop$End[i]) | (sub$Chromosome == chrom_start_stop$Chromosome[i]) & (sub$Start <= chrom_start_stop$Start[i]) & (sub$End >= chrom_start_stop$Start[i]) | (sub$Chromosome == chrom_start_stop$Chromosome[i]) & (sub$Start <= chrom_start_stop$End[i]) & (sub$End >= chrom_start_stop$End[i]) | (sub$Chromosome == chrom_start_stop$Chromosome[i]) & (sub$Start >= chrom_start_stop$Start[i]) & (sub$End <= chrom_start_stop$End[i]),]
                chrom_start_stop[i,name] = mean(CN_rows$Modal_Total_CN)
	}
}
		
save(chrom_start_stop, file="chrom_start_stop_copynumbermean.RData")

chrom_start_stop = round(chrom_start_stop)




 
