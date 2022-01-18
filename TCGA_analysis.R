# uses TCGA ABSOLUTE copy znumber calls to get 1Mb sized copy number data for samples
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
library(cowplot)
library(ggpubr)

# functions
'%!in%' <- function(x,y)!('%in%'(x,y))

#### this finds chromothriptic event as described as a oscillation of varying copynumber
findChromothripsis <- function(Cn) {
  if (length(Cn[Cn >= 5]) > 2) {
    y1 <- as.numeric(Cn > 7 & lead(Cn) %in% 1:2)
    as.numeric(Cn > 7 & lead(Cn) %in% 1:2)
  }
  else {
    as.numeric(!(Cn))
  }
}


### this finds shattering events as the number of segments in a chromosome as higher than 5
findshattering <- function(x) {
    return(length(rle(sub)$lengths))
}



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

Merge5 = round(Merge5)

Merge5 = Merge5[
  with(Merge5, order(Chromosome, Start)),
]

# load premade 1Mb TCGA samples and duodenals
#load("/Users/nowins01/Documents/Duodenal_project/Analysis/perf_cn/Merged_TCGA_DUODENALS_perfect.RData")
#Merge5$Start.y <-NULL
#Merge5$End.y <-NULL
#Merge5$Chromosome.y <-NULL
#Merge5$Start.x <-NULL
#Merge5$End.x <-NULL

### plot copynumber all samples, tiny heatmap

#Chromosome separation positions
chr.ends = cumsum(table(Merge5[,"Chromosome"]))[1:22]

filtered <- Merge5
filtered$Chromosome <- NULL
filtered$Start <- NULL
filtered$End <- NULL

filtered[filtered>7] <- 7
filtered[filtered<1] <- 1

#put them in the order you want
new = cbind(filtered[, names(filtered) %in% oes],filtered[, names(filtered) %in% stom],filtered[, names(filtered) %in% panc],filtered[, names(filtered) %in% duo], filtered[, names(filtered) %in% coad], filtered[, names(filtered) %in% rec])


filtered = new

### make the heatmap
#count number of segments with most number of aberrant bins
order_for_plot <- function(samples) {
    convert <- samples
    samples <- sapply(samples, as.numeric)
    samples[samples=="2"] <- 0
    samples[samples=="3"] <- 1
    convert <- convert[,order(colSums(samples))]
 return(convert)
}



###### make a heatmap first of just the duodenal samples

Duodenals = filtered[,names(filtered) %in% duo]
Duodenals <- order_for_plot(Duodenals[,c(1:125)])

pdf("Duodenum.pdf", width=10, height=10)

new.positions <- (diff(chr.ends))/2
n.pos <- c(113, new.positions)
n.pos <- (chr.ends - n.pos)

bot.column.anno = columnAnnotation(link = column_anno_link(at = as.vector(chr.ends) - 30,
                                                           labels = 1:22,
                                                           side = "bottom",
                                                           link_width = unit(0.25, "cm")),
                                   width = unit(2, "cm"))

samps = rep("Duodenum", sum(names(filtered) %in% duo))


df = as.data.frame(samps)
names(df) = "location"

#location = c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D")
location = "#BDD7E7"
names(location) = unique(samps)
colors = list(location = location)

CopyNumber <- Heatmap(t(Duodenals),
        name = "CopyNumber",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        # row_split = as.character(df$type),
        col = colorRamp2(c(1, 2, 3, 4,5,6,7), c("blue", "white", "pink","red", "darkred", "orange", "yellow")),
        heatmap_legend_param = list(at = c(1, 2, 3, 4,5,6,7), labels = c("CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7")))


ha <- HeatmapAnnotation(df = df, which = "row", col=colors, width = unit(1, "cm"))

draw(ha + CopyNumber, gap = unit(0.5, "cm"))


#Add lines
for(boundary in chr.ends / ncol(t(Duodenals))) {
  
  #Add the lines
  decorate_heatmap_body("CopyNumber", {
    grid.lines(c(boundary, boundary), c(0, 1), gp = gpar(lty = "dotted", lwd = 1))
  })
  
}

row_lines <- c(seq(0,nrow(t(filtered)),1))

dev.off()












pdf("heatmap_duodenal_tcga_samples_ordered_clustered.pdf", width=10, height=10)

new.positions <- (diff(chr.ends))/2
n.pos <- c(113, new.positions)
n.pos <- (chr.ends - n.pos)

bot.column.anno = columnAnnotation(link = column_anno_link(at = as.vector(chr.ends) - 30, 
                                                           labels = 1:22, 
                                                           side = "bottom",
                                                           link_width = unit(0.25, "cm")), 
                                   width = unit(2, "cm"))


samps = (c(rep("Oesophagus", sum(names(filtered) %in% oes)), rep("Stomach", sum(names(filtered) %in% stom)), rep("Pancreas", sum(names(filtered) %in% panc)),rep("Duodenum", sum(names(filtered) %in% duo)),rep("Colon", sum(names(filtered) %in% coad)), rep("Rectum", sum(names(filtered) %in% rec))))

df = as.data.frame(samps)
names(df) = "location"

#location = c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D")
location = c("#F0F0F0", "#76000D", "#FCAE91", "#BDD7E7", "#6283A9", "#B9574E")
names(location) = unique(samps)
colors = list(location = location)

CopyNumber <- Heatmap(t(filtered),
        name = "CopyNumber",
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        bottom_annotation = bot.column.anno,
        # row_split = as.character(df$type),
        col = colorRamp2(c(1, 2, 3, 4,5,6,7), c("blue", "white", "pink","red", "darkred", "orange", "yellow")),
        heatmap_legend_param = list(at = c(1, 2, 3, 4,5,6,7), labels = c("CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7")))


ha <- HeatmapAnnotation(df = df, which = "row", col=colors, width = unit(1, "cm"))

draw(ha + CopyNumber, gap = unit(0.5, "cm"))


#Add lines
for(boundary in chr.ends / ncol(t(filtered))) {
  
  #Add the lines
  decorate_heatmap_body("CopyNumber", {
    grid.lines(c(boundary, boundary), c(0, 1), gp = gpar(lty = "dotted", lwd = 1))
  })
  
}

row_lines <- c(seq(0,nrow(t(filtered)),1))


for(lines in row_lines / nrow(t(filtered))) {
decorate_heatmap_body("CopyNumber", {
    grid.lines(c(0, 1), unit(c(lines, lines), "native"), gp = gpar(col = "black"))
})
}


dev.off()




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
copy[copy >= 5] <- 5
copy[copy < 2] <- 1

## calculate diversity per bin - PIC score

## number of samples that are loss at each bin

get_diversity <- function(x) {
	(1- ((sum(x == 1)/length(x))^2 + (sum(x == 2)/length(x))^2 + (sum(x == 3)/length(x))^2 + (sum(x == 4)/length(x))^2 + (sum(x == 5)/length(x))^2))
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

diversity_ploidy_pga$project_id[diversity_ploidy_pga$project_id=="TCGA-COAD"] <- "Colon"
diversity_ploidy_pga$project_id[diversity_ploidy_pga$project_id=="TCGA-ESCA"] <- "Oesophagus"
diversity_ploidy_pga$project_id[diversity_ploidy_pga$project_id=="TCGA-PAAD"] <- "Pancreas"
diversity_ploidy_pga$project_id[diversity_ploidy_pga$project_id=="TCGA-READ"] <- "Rectum"
diversity_ploidy_pga$project_id[diversity_ploidy_pga$project_id=="TCGA-STAD"] <- "Stomach"

#diversity_ploidy_pga$project_id <- total$project_id

##################### PLOT

#my_comparisons <- list( c("Colon", "Oesophagus"), c("Colon", "Pancreas"), c("Colon", "Rectum"), c("Colon", "Stomach"), c("Colon", "Duodenum"), c("Pancreas", "Oesophagus"), c("Rectum", "Oesophagus"), c("Stomach", "Oesophagus"), c("Duodenum", "Oesophagus"), c("Pancreas", "Rectum"), c("Stomach", "Pancreas"), c("Duodenum", "Pancreas"), c("Stomach", "Rectum"), c("Duodenum", "Rectum"), c("Duodenum", "Stomach") )

my_comparisons <- list( c("Colon", "Duodenum"),c("Colon", "Oesophagus"), c("Colon", "Pancreas"), c("Colon", "Rectum"), c("Colon", "Stomach"), c("Duodenum", "Oesophagus"),c("Duodenum", "Pancreas"), c("Duodenum", "Rectum"),c("Duodenum", "Stomach"),c("Pancreas", "Oesophagus"), c("Rectum", "Oesophagus"), c("Stomach", "Oesophagus"),c("Pancreas", "Rectum"), c("Stomach", "Pancreas"),c("Stomach", "Rectum"))

pdf("Figure1_stats.pdf", width=7, height=5)
ggplot(diversity_ploidy_pga, aes(x=project_id, y=PIC, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+ylim(c(0,2.4))+
  ylab("Diversity")+ xlab(NULL)+
  theme_half_open(12)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 10))+         
  stat_compare_means(comparisons = my_comparisons)   
dev.off()

pdf("Figure2_stats.pdf", width=7, height=5)
ggplot(diversity_ploidy_pga, aes(x=project_id, y=GD, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+ylim(c(0,14.5))+
  ylab("Ploidy")+ xlab(NULL)+
  theme_half_open(12)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 10))+
  stat_compare_means(comparisons = my_comparisons)

dev.off()

pdf("Figure1.pdf", width=7, height=5)
ggplot(diversity_ploidy_pga, aes(x=project_id, y=PIC, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+ylim(c(0,1))+
  ylab("Diversity")+ xlab(NULL)+
  theme_half_open(12)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 20))
dev.off()

pdf("Figure2.pdf", width=7, height=5)
ggplot(diversity_ploidy_pga, aes(x=project_id, y=GD, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+ylim(c(1,7))+
  ylab("Ploidy")+ xlab(NULL)+
  theme_half_open(12)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 20))
dev.off()

pdf("Figure3.pdf", width=7, height=4)
ggplot(diversity_ploidy_pga, aes(x=project_id, y=PGA, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+ylim(c(0,1))+
ylab("PGA")+xlab(NULL)+
  theme_half_open(12)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 20))
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

# PIC score along genom
PIC_SCORE_bin = as.data.frame(apply(copy, 1, get_diversity))

names = total[,c(1,4)]
panc = names$sample[names$project_id == "TCGA-PAAD"]
panc_pic = as.data.frame(apply(copy[,names(copy) %in% panc], 1, get_diversity))
names(panc_pic)[1] <- "PIC_bin"
panc_pic$type <- "Pancreas"
panc_pic$num <- rownames(panc_pic)

coad = names$sample[names$project_id == "TCGA-COAD"]
coad_pic = as.data.frame(apply(copy[,names(copy) %in% coad], 1, get_diversity))
names(coad_pic)[1] <- "PIC_bin"
coad_pic$type <- "Colon"
coad_pic$num <- rownames(coad_pic)

oes = names$sample[names$project_id == "TCGA-ESCA"]
oes_pic = as.data.frame(apply(copy[,names(copy) %in% oes], 1, get_diversity))
names(oes_pic)[1] <- "PIC_bin"
oes_pic$type <- "Oesophagus"
oes_pic$num <- rownames(oes_pic)

stom = names$sample[names$project_id == "TCGA-STAD"]
stom_pic = as.data.frame(apply(copy[,names(copy) %in% stom], 1, get_diversity))
names(stom_pic)[1] <- "PIC_bin"
stom_pic$type <- "Stomach"
stom_pic$num <- rownames(stom_pic)

rec = names$sample[names$project_id == "TCGA-READ"]
rec_pic = as.data.frame(apply(copy[,names(copy) %in% rec], 1, get_diversity))
names(rec_pic)[1] <- "PIC_bin"
rec_pic$type <- "Rectum"
rec_pic$num <- rownames(rec_pic)

duo = names$sample[names$project_id == "Duodenum"]
duo_pic = as.data.frame(apply(copy[,names(copy) %in% duo], 1, get_diversity))
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
    scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
    ylab("Diversity along genome")

ggplot(total_pic_bin, aes(fill=type, y=PIC_bin, x=num)) + 
    geom_bar(stat="identity", position="fill") +
    geom_line(data=total_pic_bin, aes(x=num,y=PIC_bin, group=type), color="black") +
    scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D")) +
    ggtitle("bins_PIC")


ggplot(total_pic_bin, aes(y=PIC_bin, x=num)) +
  geom_line(aes(linetype = type))


ggplot(total_pic_bin, aes(colour=type, y=PIC_bin, x=num, group=type)) + geom_line()   
    geom_bar(stat="identity", position="fill") +
    geom_line(data=total_pic_bin, aes(x=num,y=PIC_bin, group=type))

save(total_pic_bin, file="total_pic_bin.RData")


######### see if samples have been shattered (possible chromothripsis)

chrom = as.data.frame(names(Merge5)[4:ncol(Merge5)])
names(chrom)[1] <- "sample"
chrom$CHROMOTHRIPSIS = 0

for (name in names(Merge5)[4:ncol(Merge5)]) {
	max=0
	for (chr in 1:22) {
		sub=Merge5[Merge5$Chromosome==chr,name]
		if (max < findshattering(sub)) {
			max = findshattering(sub)
		}
	chrom$CHROMOTHRIPSIS[chrom$sample==name] <- max
	}
}

new_merged = merge(diversity_ploidy_pga, chrom, by="sample")

## plot

pdf("Figure4.pdf", width=7, height=5)
ggplot(new_merged, aes(x=project_id, y=CHROMOTHRIPSIS, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+
ylab("Segments")+xlab(NULL)+
  theme_half_open(12)+theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 20))
dev.off()


total_pic_bin

pdf("Figure5.pdf", width=7, height=5)
ggplot(total_pic_bin, aes(x=type, y=PIC_bin, fill=type)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+
ylab("Diversity")+xlab(NULL)+
  theme_half_open(12)+theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 20))
dev.off()


# Plot along genome
total_pic_bin %>%
  ggplot( aes(x=num, y=PIC_bin, group=type, color=type)) +
    geom_line() +
    ggtitle("Population-wde diversity") +
    scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
    ylab("Diversity along genome")


# Plot along genome
duo_pic %>%
  ggplot( aes(x=num, y=PIC_bin, group=type, color=type)) +
    geom_line() +
    ggtitle("Population-wide diversity") +
    scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
    ylab("Diversity along genome")







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


### same analysis as above, however, COAD is split into left and right

LSCC$case_submitter_id
RSCC$case_submitter_id

diversity_ploidy_pga$project_id_extra <- diversity_ploidy_pga$project_id
diversity_ploidy_pga$project_id_extra[diversity_ploidy_pga$sample %in% LSCC$case_submitter_id] <- "Left Colon"
diversity_ploidy_pga$project_id_extra[diversity_ploidy_pga$sample %in% RSCC$case_submitter_id] <- "Right Colon"
complete = diversity_ploidy_pga[diversity_ploidy_pga$project_id_extra!="Colon",]
complete = distinct(complete)

pdf("Figure1_complete.pdf", width=7, height=5)
ggplot(complete, aes(x=project_id_extra, y=PIC, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+ylim(c(0,1))+
  ylab("Diversity")+ xlab(NULL)+
  theme_half_open(12)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 20))
dev.off()

pdf("Figure2_complete.pdf", width=7, height=5)
ggplot(complete, aes(x=project_id_extra, y=GD, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+ylim(c(1,7))+
  ylab("Ploidy")+ xlab(NULL)+
  theme_half_open(12)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 20))
dev.off()

pdf("Figure3_complete.pdf", width=7, height=4)
ggplot(complete, aes(x=project_id_extra, y=PGA, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+ylim(c(0,1))+
ylab("PGA")+xlab(NULL)+
  theme_half_open(12)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 20))
dev.off()

pdf("Figure4.pdf", width=7, height=5)
ggplot(new_merged, aes(x=project_id, y=CHROMOTHRIPSIS, fill=project_id)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+
ylab("Segments")+xlab(NULL)+
  theme_half_open(12)+theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 20))
dev.off()


pdf("Figure5_with_gd.pdf", width=7, height=5)
ggplot(total_pic_bin, aes(x=type, y=PIC_bin, fill=type)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+
ylab("Diversity")+xlab(NULL)+
  theme_half_open(12)+theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 20))
dev.off()


total_pic_bin$type_extra <- total_pic_bin$type

merge(total_pic_bin, complete, by="sample")

pdf("Figure5_with_gd.pdf", width=7, height=5)
ggplot(total_pic_bin, aes(x=type, y=PIC_bin, fill=type)) +
  geom_violin(trim=TRUE)+scale_fill_manual(values=c("#6283A9","#BDD7E7","#F0F0F0","#FCAE91", "#B9574E", "#76000D"))+
geom_boxplot(width=0.1, fill="white")+
ylab("Diversity")+xlab(NULL)+
  theme_half_open(12)+theme(axis.text.x = element_text(angle = 45, vjust = 0.5),text = element_text(size = 20))
dev.off()




##############################
##############################

# amplifications

filtered





 
