### check if # of TEs are significantly different ###

library(ggplot2)
library(stringr)
library(RColorBrewer)
library(dplyr)
library(reshape)
library(forcats)
library(janitor)

# set directory
dir <- "/Users/heidiyang/tetragnatha_project/dnaPipeTE_counts_data/Counts"
setwd(dir)

# read tables in
files <- list.files(path=dir, full.names=FALSE, recursive=FALSE)
ldf <- lapply(files, read.table)

## clean up the data ##
# convert list of dfs into columns of single df
compiled_df <-  do.call("cbind", ldf)
# remove duplicated class columns
counts_df <- compiled_df[!duplicated(as.list(compiled_df))]
# name columns specimen name 
new_names <- lapply(files, function(x) {substr(x, 1, nchar(x)-13)})
colnames(counts_df) <- c('class',new_names)
# remove rows with 0s  
counts_df <- counts_df[-c(5, 10, 14),] 
# reset row index
rownames(counts_df) <- seq(length=nrow(counts_df))

## add single or low copy sequences to df ##
# totals of all the sequences to subtract from calculated Total
sum_of_rows <- colSums(counts_df[-12,-1])
# add Single or low copy row
counts_df[nrow(counts_df) + 1,] <- c('Non-repetitive', counts_df[12,-1] - sum_of_rows)

## calculate genome percentage and add to df ##
# get total row, remove from df, and reset row index
totals <- counts_df[12,]
counts_df <- counts_df[-12,]
rownames(counts_df) <- seq(length=nrow(counts_df))
# for each specimen calculate percentage
for (c in 2:ncol(counts_df)) {
  counts_df[,c] <- round(as.numeric(counts_df[,c] / totals[,c]), 4) * 100
}

# get TE average for radiating and non-radiating
radiating_TEs <- counts_df[-c(6,7,8,9,10,11,12,17),-c(2,9,15,26)]
all_radiating_TEs <- colSums(radiating_TEs[,-1])
avg_TE_rad <- mean(all_radiating_TEs)

non_TEs <- counts_df[-c(6,7,8,9,10,11,12,17),c(2,9,15,26)]
all_non_TEs <- colSums(non_TEs[,-1])
avg_TE_non <- mean(all_non_TEs)

# add metadata
counts_df[nrow(counts_df)+1,] <- lapply(colnames(counts_df), function(x) {strsplit(x,"_")[[1]][3]})
counts_df[nrow(counts_df)+1,] <- lapply(colnames(counts_df), function(x) {strsplit(x,"_")[[1]][4]})
counts_df[nrow(counts_df)+1,] <- lapply(colnames(counts_df), function(x) {strsplit(x,"_")[[1]][5]})
counts_df[nrow(counts_df)+1,] <- colnames(counts_df)
counts_df[13:16,1] <- c('island', 'volcano', 'ecomorph', 'specimen')
t_df <- data.frame(t(counts_df))
colnames(t_df) <- t_df[1,]
t_df <- t_df[-1,]
rownames(t_df) <- NULL

# check for normality
shapiro.test(as.numeric(t_df$DNA))

# melt dataframe
m_df <- melt(t_df, id=c("island", 'volcano', 'ecomorph', 'specimen'))
colnames(m_df) <- c("island", 'volcano', 'ecomorph', 'specimen', 'class', 'perc')

# radiating only df
mr_df <- m_df[!(m_df$specimen =="T_acuta_Maui_EM_lightbrown_071" | m_df$specimen == 'T_filiciphila_Maui_EM_outgroup_074' | m_df$specimen == 'T_maka_Kauai_Ka_lightbrown_052' | m_df$specimen == 'T_stellarobusta_Maui_EM_maroon_001'),]
mr_df$rad_or_not <- 'rad'
mn_df <- m_df[(m_df$specimen =="T_acuta_Maui_EM_lightbrown_071" | m_df$specimen == 'T_filiciphila_Maui_EM_outgroup_074' | m_df$specimen == 'T_maka_Kauai_Ka_lightbrown_052' | m_df$specimen == 'T_stellarobusta_Maui_EM_maroon_001'),]
mn_df$rad_or_not <- 'not'
all_df <- rbind(mr_df, mn_df)
all_df$perc <- as.numeric(all_df$perc)

## time for da stats! ##
volc_ltr <- kruskal.test(perc ~ volcano, data = mr_df[mr_df$class == 'LTR',])

volc_line <- kruskal.test(perc ~ volcano, data = mr_df[mr_df$class == 'LINE',])

volc_sine <- kruskal.test(perc ~ volcano, data = mr_df[mr_df$class == 'SINE',])

volc_dna <- kruskal.test(perc ~ volcano, data = mr_df[mr_df$class == 'DNA',])

volc_hel <- kruskal.test(perc ~ volcano, data = mr_df[mr_df$class == 'Helitron',])
s_volc_hel <- summary(volc_hel)

volc_rRNA <- kruskal.test(perc ~ volcano, data = mr_df[mr_df$class == 'rRNA',])
s_volc_rRNA <- summary(volc_rRNA)

volc_lc <- kruskal.test(perc ~ volcano, data = mr_df[mr_df$class == 'Low_Complexity',])
s_volc_lc <- summary(volc_lc)

volc_sat <- kruskal.test(perc ~ volcano, data = mr_df[mr_df$class == 'Satellite',])
s_volc_sat <- summary(volc_sat)

volc_sr <- kruskal.test(perc ~ volcano, data = mr_df[mr_df$class == 'Simple_repeat',])
s_volc_sr <- summary(volc_sr)

# for the web-builder vs. spiny leg
mw <-wilcox.test(perc ~ rad_or_not, data=all_df[all_df$class == 'Simple_repeat',], na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
print(mw)

# for the young TEs average testing

