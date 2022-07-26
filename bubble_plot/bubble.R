library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr)
library(forcats)

mydat <- read.csv2("GO_DAVID_Summary_R.csv", sep=",", header=TRUE)
#mydat$GeneRatio <- sapply(mydat$GeneRatio_frac, function(x) eval(parse(text=x)))
mydat$FoldEnrichment <- as.numeric(mydat$FoldEnrichment)
mydat$Comparison <- factor(mydat$Comparison, levels=c("sT 3 dpi","sT 8 dpi","LT 3 dpi","LT 8 dpi","LTtr 3 dpi","LTtr 8 dpi","sT+LT 3 dpi","sT+LTtr 9/12 dpi"))    # solution for point 2
mydat$Term <- factor(mydat$Term, levels=rev(c("Cell proliferation","Negative regulation of cell proliferation","Cell division","Mitotic nuclear division","G1/S transition of mitotic cell cycle","DNA replication initiation","DNA replication","Sister chromatid cohesion","DNA repair","Cellular response to DNA damage stimulus","Transcription, DNA-templated","Regulation of transcription, DNA-templated","Positive regulation of transcription, DNA-templated","Positive regulation of transcription from RNA polymerase II promoter","Positive regulation of gene expression","Negative regulation of transcription from RNA polymerase II promoter","rRNA processing","Protein folding","Inflammatory response","Immune response","Innate immune response","Positive regulation of ERK1 and ERK2 cascade","Chemokine-mediated signaling pathway","Chemotaxis","Cell chemotaxis","Neutrophil chemotaxis","Viral process","Response to virus","Defense response to virus","Cellular response to lipopolysaccharide","Type I interferon signaling pathway","Extracellular matrix organization","Cell adhesion","Nervous system development","Angiogenesis","Apoptotic process")))    # solution for point 2
mydat$Regulation <- factor(mydat$Regulation, levels=c("up","down"))
png("bubble.png", width=1166, height=1067)
ggplot(mydat, aes(y = Term, x = Comparison, size = FoldEnrichment)) + geom_point(aes(color = Regulation), alpha = 1.0) + labs(x = "", y = "") + theme(axis.text.x = element_text(angle = 30, vjust = 0.5)) + theme(axis.text = element_text(size = 20)) + theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20)) + scale_size(range = c(1, 20)) + guides(color = guide_legend(override.aes = list(size = 10)))    # solution for point 1
dev.off()
