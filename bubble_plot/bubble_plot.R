# #-1- https://www.mzes.uni-mannheim.de/socialsciencedatalab/article/datavis/#bubble
# https://github.com/SocialScienceDataLab/Data_Visualization/blob/master/Code/Multivariate.R
# 
# https://stackoverflow.com/questions/60025301/bubble-chart-with-r
# 
# 
# 
#-2- https://stackoverflow.com/questions/60025301/bubble-chart-with-r
library(tidyverse)
library(ggplot2)

df1 %>%
  rownames_to_column(var = "id") %>%
  gather(key, Abundance, -id) %>%
  ggplot(aes(key, id)) +
  geom_point(aes(size = Abundance), colour = "red", fill = "red", shape = 21) +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


##I have another idea, do you think it might be possible to differentiate between upregulated and downregulated GO terms with blue and red color, instead of showing the padj. Values (because I also filter for only significant terms, I don’t need to show the padj. Values in the figure). Maybe you have an idea how to do that.
#https://www.quora.com/What-color-does-blue-and-red-make
#blue + red = violet (or magenta)



#ID	Description	GeneRatio	BgRatio	pvalue	p.adjust	qvalue	geneID	Count
#GO:0042254	ribosome biogenesis	35/413	297/18670	6.13365149901527E-16	2.46879472835365E-12	2.2132797198552E-12	55759/708/9188/9221/25879/54552/23076/10885/10360/22984/55127/2091/55720/5036/6059/51154/50628/55646/5822/4869/9816/51491/4839/285855/54881/9875/23404/64965/705/23195/10171/6152/23560/23517/81875	35


#-3- https://www.biostars.org/p/463626/
library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2) #for plotting
library(forcats) #for plotting

#Toy data.frame
mydat <- structure(list(GO_term = structure(c(2L, 4L, 3L, 1L, 4L), 
                                            .Label = c("Kinase", "Metabolism", "Nucleus", "Photosynthesis"), 
                                            class = "factor"), 
                        Number = c(5, 10, 15, 16, 20), 
                        Class = structure(c(3L, 2L, 1L, 1L, 2L),
                                          .Label = c("hs", "hzs", "start_duf"), class = "factor"), 
                        Type = structure(c(1L, 1L, 2L, 3L, 1L), 
                                         .Label = c("BP", "CC", "MF"), 
                                         class = "factor")), 
                   class = "data.frame", row.names = c(NA, 5L))


#sapply(mydat, class)
mydat <- read.csv2("bubble_input-1.csv", sep=",", header=TRUE)
mydat$GeneRatio <- sapply(mydat$GeneRatio_frac, function(x) eval(parse(text=x)))

#https://statisticsglobe.com/change-font-size-of-ggplot2-plot-in-r-axis-text-main-title-legend
png("bubble.png", width=600, height=600)
ggplot(mydat, aes(y = Description, x = Cmp, size = GeneRatio)) + geom_point(aes(color = Regulation), alpha = 1.0) + labs(x = "", y = "") + theme(axis.text.x = element_text(angle = 15, vjust = 0.5)) + theme(axis.text = element_text(size = 20)) + theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20)) 
dev.off()        

+ geom_tile(aes(width = Inf, fill = Type), alpha = 0.4)
+ scale_fill_manual(values = c("green", "red", "blue"))
                   
                   
#First we group by Type and GO_term, and assign a "yes" to the first row
#and "no" to every other row of the grouping
mydat %<>% 
  group_by(Type, GO_term) %>%
  mutate(typefill = if_else(row_number() == 1, "yes", "no")) %>%
  ungroup()
#Then in the whole data.frame, typefill = "yes" will be replaced by the Type value
#from that row, and typefill = "no" will be replaced with NA
mydat %<>% mutate(typefill = ifelse(typefill == "yes", as.character(Type), NA))



#Plotting, now pass typefill to geom_tile's fill parameter instead of Type
#ggplot(mydat, aes(y = reorder(GO_term, as.numeric(Type)), x = Number, size = Number)) + geom_point(aes(color = Class), alpha = 1.0) + 
ggplot(mydat, aes(y = reorder(GO_term, as.numeric(Type)), x=Number, size = Number)) + geom_point(aes(color = Class), alpha = 1.0) + 
  geom_tile(aes(width = Inf, fill = typefill), alpha = 0.4) + 
  scale_fill_manual(values = c("green", "red", "blue"))

