# load libraries
library(ggplot2)
library(phytools)
library(stringr)
library(janitor)
library(RColorBrewer)
library(tibble)
library(dplyr)
library(tidyverse)

# read in tree
skmer_tree <- read.tree(text="((((((T_restricta_037:0.018174980259,(T_kamakou_012:0.010249524294,T_kamakou_010:0.011278475706):0.003537957530):0.001856675030,T_perreirai_024:0.018716626566):0.000243759786,((((T_stellarobusta_001:0.014232560772,T_acuta_071:0.018628439228):0.003602556638,(T_maka_052:0.020423480000,T_filiciphila_074:0.018694520000):0.004434656856):0.015690222578,(T_pilosa_028:0.022364093045,(T_mohihi_022:0.020807753827,T_kauaiensis_020:0.017653246173):0.009150107866):0.001129780860):0.005997369111,((((T_quasimodo_064:0.009295299696,T_quasimodo_002:0.010551700304):0.000636522888,(T_quasimodo_005:0.004731476823,T_quasimodo_004:0.004373523177):0.007030647009):0.000789553668,T_quasimodo_033:0.012912822404):0.004610080835,(((T_obscura_006:0.013496084426,T_kukuiki_058:0.013483915574):0.002538304823,T_kikokiko_072:0.015804305944):0.001282938358,(T_anuenue_066:0.001756558883,T_anuenue_048:0.001982441117):0.013962470795):0.003945370729):0.004418508119):0.001258957099):0.003403982199,T_tantalus_039:0.016658297207):0.001566976360,T_brevignatha_059:0.014570826887):0.000578645048,((T_brevignatha_067:0.003191322367,T_brevignatha_007:0.003111677633):0.001666783851,T_brevignatha_050:0.004198479965):0.011108235210,T_waikamoi_062:0.011411133039);")

# root tree to outgroup
skmer_tree <- root(skmer_tree, 7) # the seventh lil nested friends

### from the stackedbar.R script ###

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

## combine others and na ##
counts_df[13,] <- c("Unclassified", (counts_df[10,-1] + counts_df[11,-1]))
total <- counts_df[12,]
counts_df <- counts_df[-c(10,11),]
rownames(counts_df) <- 1:nrow(counts_df)

## add single or low copy sequences to df ##
# totals of all the sequences to subtract from calculated Total
sum_of_rows <- colSums(counts_df[-12,-1])
# add Single or low copy row
counts_df[nrow(counts_df) + 1,] <- c('Non-repetitive', (sum_of_rows - counts_df[10,-1]))

# calculate percentage of genome that are TEs
total <- counts_df[10,]
counts_df <- counts_df[-10,]
counts_df[nrow(counts_df) + 1,] <- total
rownames(counts_df) <- NULL

# make genomic repeat df
genomic_repeat_content <- data.frame(apply(counts_df[,-1], 2, function(x) (sum(x[1:9]) / x[12]) * 100))
genomic_repeat_content["ind"] <- rownames(genomic_repeat_content)
colnames(genomic_repeat_content) <- c("perc", "ind")
grc_df <- genomic_repeat_content#[-c(1,8,14,25),]
list_of_names <- data.frame(sapply(rownames(grc_df), function (x) {strsplit(x, "_")})) 
grc_df$volcano <- sapply(1:ncol(list_of_names), function(i) {paste(list_of_names[c(3,4),i],collapse="_")}) 
grc_df$individual <- sapply(1:ncol(list_of_names), function(i) {paste(list_of_names[c(1,2,6),i],collapse="_")}) 
grc_df$volcano <- replace(grc_df$volcano, c(1,8,14,25), "web-builders")

# relative proportion of TEs 
rel_df <- counts_df[-c(10,11,12),]
rel_df <- rel_df %>%
  adorn_totals("row")
rownames(rel_df) <- seq(length=nrow(rel_df))
# for each specimen calculate percentage
for (c in 2:ncol(rel_df)) {
  rel_df[,c] <- round(as.numeric(rel_df[,c] / rel_df[10,c]), 4) * 100
}
rel_df <- rel_df[-11,]

# may need to transpose data frame
t_rel_df <- data.frame(t(rel_df))
colnames(t_rel_df) <- t_rel_df[1,]
t_rel_df <- t_rel_df[-1, -10]

## stats for volcano age and TE content ##
volcano <- c("Kauai_Ka", "Oahu_Wai", "Oahu_Koo", "Molokai_Mo", "Maui_EM", "Hawaii_Koh", "Hawaii_ML", "Hawaii_Kil", "web-builders")
volcano_ages <- c(4e06, ((3.9e06 + 2.5e06)/2), 1.7e06, 1.3e06, ((1e06 + 7.5e05)/2), 1e06, ((7e05 + 1e06)/2), ((2.1e05 + 2.8e05)/2), 0)
volcano_df <- data.frame(volcano, volcano_ages)
volcano_ages_df <- merge(grc_df, volcano_df, by = "volcano")

# add volcano names and data
list_of_names <- data.frame(sapply(rownames(t_rel_df), function (x) {strsplit(x, "_")})) 
t_rel_df$volcano <- sapply(1:ncol(list_of_names), function(i) {paste(list_of_names[c(3,4),i],collapse="_")}) 
t_rel_df$individual <- sapply(1:ncol(list_of_names), function(i) {paste(list_of_names[c(1,2,6),i],collapse="_")}) 
t_rel_df <- t_rel_df[-c(1,8,14,25),] # remove outgroups
volcano_ages_df <- merge(t_rel_df, volcano_df, by = "volcano")
#specimen_bar_order <- unlist(young_te_df$specimen)

# summarize df
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
plt_df <- data_summary(volcano_ages_df, varname="perc", 
                    groupnames=c("volcano", "volcano_ages"))


# correlation test
vol_corr_df <- volcano_ages_df %>% group_by(volcano) %>% summarise(avg_repeat_content = mean(as.numeric(perc)), volcano_age = mean(volcano_ages))
volcano_ages_df$volcano <- factor(volcano_ages_df$volcano, levels = c("Kauai_Ka", "Oahu_Wai", "Oahu_Koo", "Molokai_Mo", "Maui_EM", "Hawaii_Koh", "Hawaii_ML", "Hawaii_Kil", "web-builders"))
volcano_order <- c("Kauai_Ka", "Oahu_Wai", "Oahu_Koo", "Molokai_Mo", "Maui_EM", "Hawaii_Koh", "Hawaii_ML", "Hawaii_Kil", "web-builders")
volcano_ages_df <- volcano_ages_df[order(match(volcano_ages_df$volcano, volcano_order)),]

corr <- cor.test(vol_corr_df$volcano_age, vol_corr_df$avg_repeat_content, method="pearson")
corr2 <- cor.test(volcano_ages_df$volcano_age, as.numeric(volcano_ages_df$Total), method="pearson")
corr2
ggplot(plt_df, aes(x=volcano_ages, y=perc)) + 
  labs(x = "volcano age", y = "genomic repeat percentage") + 
  scale_x_reverse() + # reverse volcanic age
  geom_point(color="blue") + geom_line(color="black") + theme_bw() +
  geom_errorbar(aes(ymin=perc-sd, ymax=perc+sd), width=.2,
                position=position_dodge(0.05))

# boxplot
ggplot(volcano_ages_df, aes(x=factor(volcano, levels=volcano_order), y=perc, fill=factor(volcano, levels=volcano_order))) + xlab("volcano") + ylab("genomic repeat percentage") +
  geom_boxplot(notch=FALSE, outlier.shape=1) + theme_bw() + theme(axis.text.x = element_text(angle=90)) +
  scale_fill_brewer(palette="BuPu", direction = -1, breaks = c("Kauai_Ka", "Oahu_Wai", "Oahu_Koo", "Molokai_Mo", "Maui_EM", "Hawaii_Koh", "Hawaii_ML", "Hawaii_Kil", "web-builders"), name="volcano")
# for average repeat content per volcano
ggplot(vol_corr_df, aes(x=volcano, y=avg_repeat_content)) + geom_bar(stat="identity", fill="#967bb6") + 
  #scale_color_brewer(palette="BuPu", direction = -1, breaks = c("Kauai_Ka", "Oahu_Koo", "Oahu_Wai", "Molokai_Mo", "Maui_EM", "Hawaii_Koh", "Hawaii_ML", "Hawaii_Kil")) +
  theme_bw() + labs(x = "volcano", y = "average repeat content")


## back to the program here ##
# transpose counts_df
t_counts_df <- data.frame(t(counts_df))
colnames(t_counts_df) <- t_counts_df[1,]
totals <- t_counts_df[, 10]
t_counts_df <- t_counts_df[-1, -c(10,11,12)]

# manipulate strings to match the tree names
list_of_names <- data.frame(sapply(rownames(t_counts_df), function (x) {strsplit(x, "_")})) 
new_names <- sapply(1:ncol(list_of_names), function(i) {paste(list_of_names[c(1,2,6),i],collapse="_")}) # 1, 2, 6 are the desired parts of the string
rownames(t_counts_df) <- new_names

# change the names for the individuals - how to change them for the other ones!!
t_counts_vol_df <- merge(t_counts_df, spec_vol_names, by="specimen")
t_counts_vol_df[t_counts_vol_df == "web-builders"] <- "web-builders_WB"
species <- str_split_i(t_counts_df$specimen, pattern="_", 2)
volc <- data.frame(str_split_fixed(t_counts_vol_df$volcano, pattern="_", 2))
volc["specimen"] <- t_counts_vol_df$specimen
volc["island"] <- sapply(volc$X1, function(x) substr(x, start = 1, stop = 2))
volc["volc_isl"] <- paste("(", volc$island, "-", volc$X2, ")", sep="")
volc_spec_df <- volc %>% 
  mutate(volc_isl=ifelse(grepl("we", volc_isl, fixed=TRUE),gsub("we-","",volc_isl),volc_isl)) %>%
  mutate(volc_isl=ifelse(grepl("Mo", volc_isl, fixed=TRUE),gsub("Mo-","",volc_isl),volc_isl)) %>%
  mutate(volc_isl=ifelse(grepl("Ka", volc_isl, fixed=TRUE),gsub("Ka-","",volc_isl),volc_isl))
volc <- paste("T", species, volc_spec_df$volc_isl, sep="_")
 

# colors for bars
#bar_colors <- colorRampPalette(brewer.pal(8, "Paired"))(12)
bar_colors <- c("#88CCEE", "#CC6677", "#DDCC77", "#009E73", "#332288", "#AA4499", 
                                           "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")  
repeat_colors <- rep(bar_colors, 27) # need to have same pattern for each species

# plot it!
plotTree.barplot(skmer_tree,t_counts_df,
                 args.barplot=list(col=repeat_colors,
                                   border="gray27",
                                   args.axis=list(at=seq(0,1,by=0.2))))
mtext("repeat abundance",1,at=50,line=2.5)

## I don't know how to get this part outside of the figure! ##
legend(x="top", legend=colnames(t_counts_df),
       pch=22,pt.cex=2,pt.bg=repeat_colors, title="Repetitive Element")






