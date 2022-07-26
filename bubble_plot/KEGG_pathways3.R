library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr)
library(forcats)

mydat <- read.csv2("KEGG_pathways3.txt", sep="\t", header=TRUE)
mydat$GeneRatio <- sapply(mydat$GeneRatio, function(x) eval(parse(text=x)))
#mydat$FoldEnrichment <- as.numeric(mydat$FoldEnrichment)
mydat$P.adjust <- as.numeric(mydat$p.adjust)
mydat$Comparison <- factor(mydat$Comparison, levels=c("WAC vs mock","WAP vs mock","deltaYopT vs mock","deltaYopO vs mock","deltaYopQ vs mock","deltaYopE vs mock","deltaYopH vs mock","pTTSSplusYopQ vs pTTSS","WACplusInh vs WAC","deltaYopQ vs WAC","WAC vs WAP","deltaYopQ vs WAP","deltaYopH vs WAP","deltaYopE vs WAP","deltaYopT vs WAP","pTTSS vs WAP","pTTSSplusYopQ vs WAP"))    # solution for point 2
#mydat$Term <- factor(mydat$Term, levels=rev(c("Cell proliferation","Negative regulation of cell proliferation","Cell division","Mitotic nuclear division","G1/S transition of mitotic cell cycle","DNA replication initiation","DNA replication","Sister chromatid cohesion","DNA repair","Cellular response to DNA damage stimulus","Transcription, DNA-templated","Regulation of transcription, DNA-templated","Positive regulation of transcription, DNA-templated","Positive regulation of transcription from RNA polymerase II promoter","Positive regulation of gene expression","Negative regulation of transcription from RNA polymerase II promoter","rRNA processing","Protein folding","Inflammatory response","Immune response","Innate immune response","Positive regulation of ERK1 and ERK2 cascade","Chemokine-mediated signaling pathway","Chemotaxis","Cell chemotaxis","Neutrophil chemotaxis","Viral process","Response to virus","Defense response to virus","Cellular response to lipopolysaccharide","Type I interferon signaling pathway","Extracellular matrix organization","Cell adhesion","Nervous system development","Angiogenesis","Apoptotic process")))    # solution for point 2
#mydat$Regulation <- factor(mydat$Regulation, levels=c("up","down"))
#
png("KEGG_pathways3.png", width=1060, height=1400)
ggplot(mydat, aes(y = Description, x = Comparison, size = GeneRatio)) + geom_point(aes(color = P.adjust), alpha = 1.0) +  labs(x = "", y = "") + theme(axis.text.x = element_text(angle = 30, vjust = 0.5)) + theme(axis.text = element_text(size = 14)) + theme(legend.text = element_text(size = 15)) + theme(legend.title = element_text(size = 15)) + scale_size(range = c(1, 8)) + guides(color = guide_legend(override.aes = list(size = 6)))    # solution for point 1
dev.off()
