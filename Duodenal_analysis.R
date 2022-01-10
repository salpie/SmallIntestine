files = Sys.glob("/Users/nowins01/Documents/Duodenal_project/lpWGS/jupyter/*pcf_segments.txt")


total_ploidy = data.frame()
total_PGA = data.frame()
total_gd = data.frame()
total_amplifications = data.frame()
total_highamplifications = data.frame()

for (i in 1:length(files)) {
    
    test = read.table(files[i], header=T)
    hmm = as.data.frame(lapply(test, rep, test$Num_Bins))
    if (sum(hmm$Copies < 1)>0) {
        hmm$Copies <- hmm$Copies+1
    }
    ploidy = mean(hmm$Copies)
    total_ploidy = rbind(total_ploidy, ploidy)
    PGA = length(hmm$Copies[hmm$Copies!=2])/length(hmm$Copies)
    total_PGA = rbind(total_PGA, PGA)
    gd = median(hmm$Copies)
    total_gd = rbind(total_gd, gd)
    if (sum(hmm$Copies < 1)>0) {
        hmm$Copies <- hmm$Copies+1
        print(files[i])
    }
    amplifications = sum(hmm$Copies > 5)
    total_amplifications = rbind(total_amplifications, amplifications)
    highamplifications = sum(hmm$Copies > 7)
    total_highamplifications = rbind(total_highamplifications, highamplifications)
 }


load("/Users/nowins01/Documents/HIPEC/QDNAseq/RData/Segment_Mean_multipcf.RData")

    test=read.table(files[1], sep="\t", header=T)
	new_segments_ACE <- as.data.frame(lapply(test, rep, test$Num_Bins))
    new_segments_ACE[,c(1:3)] <- Segment_Mean[,c(1:3)]
    comp = new_segments_ACE[,c(1:3,8)]

for (i in files[-1]) {
    test=read.table(i, sep="\t", header=T)
	new_segments_ACE <- as.data.frame(lapply(test, rep, test$Num_Bins))
    new_segments_ACE[,c(1:3)] <- Segment_Mean[,c(1:3)]
    comp = cbind(comp, new_segments_ACE[,8])
}


names= gsub("/Users/nowins01/Documents/Duodenal_project/lpWGS/jupyter/", "", files)
names= gsub("_single_pcf_segments.txt", "", names)


names(comp)[4:ncol(comp)] <- names
save(comp, file="Complete_Duodenal_Samples.RData")

