library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr)
library(forcats)

setwd("/home/jhuang/DATA/Data_Emilia_MeDIP")

mydat <- read.csv2("KEGG_4.csv", sep="\t", header=TRUE)
mydat$GeneRatio <- sapply(mydat$enrichmentRatio, function(x) eval(parse(text=x)))
#mydat$FoldEnrichment <- as.numeric(mydat$FoldEnrichment)
mydat$P.adjust <- as.numeric(mydat$FDR)
mydat$Comparison <- factor(mydat$Comparison, levels=c("hypomethylated","hypermethylated"))    # solution for point 2
mydat$description <- factor(mydat$description, levels=rev(c("Axon guidance","Long-term depression","Vascular smooth muscle contraction","Oxytocin signaling pathway","AGE-RAGE signaling pathway in diabetic complications","Proteoglycans in cancer","Tight junction","Circadian entrainment","Glutamatergic synapse","Hypertrophic cardiomyopathy (HCM)","Gastric acid secretion","Phosphatidylinositol signaling system","Fc gamma R-mediated phagocytosis","Platelet activation","Dilated cardiomyopathy (DCM)","Rap1 signaling pathway","Calcium signaling pathway","Phospholipase D signaling pathway","Adrenergic signaling in cardiomyocytes","Arrhythmogenic right ventricular cardiomyopathy (ARVC)")))    # solution for point 2
#mydat$Regulation <- factor(mydat$Regulation, levels=c("up","down"))
#
#width=600, , height=500
svg("KEGG_4.svg")
ggplot(mydat, aes(y = description, x = Comparison, size = GeneRatio)) + geom_point(aes(color = P.adjust), alpha = 1.0) +  labs(x = "", y = "") + theme(axis.text.x = element_text(angle = 15, vjust = 0.5)) + theme(axis.text = element_text(size = 14)) + theme(legend.text = element_text(size = 15)) + theme(legend.title = element_text(size = 15)) + scale_size(range = c(1, 8)) + guides(color = guide_legend(override.aes = list(size = 6)))    # solution for point 1
dev.off()



