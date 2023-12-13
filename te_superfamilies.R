# load libraries and setwd
library(ggplot2)
library(stringr)
library(viridis)
library(data.table)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)

# needs to have reads_landscape and blast_reads.counts in it
dir <- "/Users/heidiyang/tetragnatha_project/landscape"
setwd(dir)

# get landscape and count files
land_files <- list.files(path=dir, pattern="*_landscape", full.names=FALSE, recursive=FALSE)
count_files <- list.files(path=dir, pattern="*.counts", full.names=FALSE, recursive=FALSE)
all_files <- list(land_files, count_files)
file_df <- as.data.frame(do.call(cbind, all_files))
file_df$specimen <- substr(file_df$V2, 1, nchar(file_df$V2)-21)
rownames(file_df) <- NULL

superfam <- function(row) {
  # get files
  land_name <- as.character(row[1])
  land <- read.table(land_name)
  count_name <- as.character(row[2])
  count <- read.table(count_name)
  # parse sample names to get species and volcano
  specimen <- row[3]
  sep <- strsplit(as.character(specimen), "_")
  spec <- paste(sep[[1]][c(1,2,6)],collapse="_")
  volcano <- paste(sep[[1]][3:4],collapse="_")
  # name columns
  colnames(land) <- c("id", "name", "classification", "fam")
  land <- land[-4] # remove fam bc it's unnecessary
  # calculate divergence
  land$div <- 100-land[,1]
  # split column to get order and superfamily
  land[c('order', 'superfamily')] <- str_split_fixed(land$classification, '/', 2)
  land['order_superfam'] <- str_split_i(land$classification, '-', 1)
  land <- land[(str_sub(land$superfamily, -1)) != "?", ] # remove classifications that are not certain
  # add species name and volcano name - have to be repeated for every row
  land["individual"] <- rep(spec, length(land$superfamily))
  land["volcano"] <- rep(volcano, length(land$superfamily))
  return(land)
}

# get TE superfamilies
land_fam_values <- apply(file_df, 1, function(x) {superfam(x)}) # apply function defined above
land_fam_df <- data.frame(do.call(rbind, land_fam_values)) # bind all of the dfs together
land_fam_df <- land_fam_df %>% # change sister group volcano names
  mutate(volcano = ifelse(individual == "T_acuta_071", "web-builders", volcano)) %>%
  mutate(volcano = ifelse(individual == "T_filiciphila_074", "web-builders", volcano)) %>%
  mutate(volcano = ifelse(individual == "T_maka_052", "web-builders", volcano)) %>%
  mutate(volcano = ifelse(individual == "T_stellarobusta_001", "web-builders", volcano))
col_order <- c("individual", "volcano", "id", "div", "name", "classification", 
               "order_superfam", "order", "superfamily")
land_fam_df <- land_fam_df[, col_order] # organize columns

# fix names
land_fam_df[land_fam_df == "web-builders"] <- "web-builders_WB"
species <- str_split_i(land_fam_df$individual, pattern="_", 2)
volc <- data.frame(str_split_fixed(land_fam_df$volcano, pattern="_", 2))
volc["individual"] <- land_fam_df$individual
volc["island"] <- sapply(volc$X1, function(x) substr(x, start = 1, stop = 2))
volc["volc_isl"] <- paste("(", volc$island, "-", volc$X2, ")", sep="")
vol_order <- c("Kauai_Ka", "Oahu_Wai", "Oahu_Koo", "Molokai_Mo", "Maui_EM", 
               "Hawaii_Koh", "Hawaii_ML", "Hawaii_Kil", "web-builders")
volc <- volc %>% 
  mutate(volc_isl=ifelse(grepl("we", volc_isl, fixed=TRUE),
                         gsub("we-","",volc_isl),volc_isl)) %>%
  mutate(volc_isl=ifelse(grepl("Mo", volc_isl, fixed=TRUE),
                         gsub("Mo-","",volc_isl),volc_isl)) %>%
  mutate(volc_isl=ifelse(grepl("Ka", volc_isl, fixed=TRUE),
                         gsub("Ka-","",volc_isl),volc_isl))
land_fam_df$volc_isl <- paste("T.", species, volc$volc_isl, sep=" ")
land_fam_df <- land_fam_df %>%
  mutate(volc_isl=ifelse(grepl("filiciphila", volc_isl, fixed=TRUE),
                         gsub("filiciphila","filiciphilia",volc_isl),volc_isl)) %>%
  mutate(volc_isl=ifelse(grepl("stella", volc_isl, fixed=TRUE),
                         gsub("stella","stela",volc_isl),volc_isl)) %>%
  mutate(volcano=ifelse(grepl("web-builders_WB", volcano, fixed=TRUE),
                        gsub("_WB","",volcano),volcano)) %>%
  mutate(volcano =  factor(volcano, levels = vol_order)) %>%
  arrange(volcano)
volc_isl_order <- unique(land_fam_df$volc_isl)

# count the number of superfamilies per individual
count_ind_superfam_df <- land_fam_df %>% 
  count(volc_isl, order_superfam)
colnames(count_ind_superfam_df) <- c("individual", "order_superfam", "count")

# get the total number of elements that exist in each of them
count_ind_tes_df <- land_fam_df %>%
  count(volc_isl)
colnames(count_ind_tes_df) <- c("individual", "total")

counts_te_df <- merge(count_ind_superfam_df, count_ind_tes_df, by="individual")
counts_te_df$perc <- (counts_te_df$count / counts_te_df$total) * 100

avg_superfam_perc <- counts_te_df %>%
  group_by(order_superfam) %>%
  summarize(mean = mean(perc, na.rm=TRUE)) %>%
  arrange(desc(mean))

# get the top 20 TE superfamilies - can't include them all
count_superfam_df <- count_ind_superfam_df %>%
  group_by(order_superfam) %>% 
  summarise(count = sum(count)) %>%
  arrange(desc(count)) %>%
  filter(!(order_superfam == "DNA")) %>% # did this bc for some reason some only had order
  slice(1:20)
top_20_TEs <- count_superfam_df$order_superfam # get the top 20 TE superfams and filter the df for them
top_20_TE_df <- count_ind_superfam_df %>%
  filter(order_superfam %in% top_20_TEs) %>%
  arrange(desc(order_superfam))
te_order <- unique(top_20_TE_df$order_superfam)

# get volcano and species orders for graphing
spec_order_df <- land_fam_df %>%
  mutate(volcano =  factor(volcano, levels = vol_order)) %>%
  arrange(volcano)
spec_order <- unique(spec_order_df$individual)


brewer.pal(9, "BuPu") # get the hex codes

# dotplot - I think this makes my eyes hurt lmao
ggplot(top_20_TE_df, aes(x=factor(individual, levels=volc_isl_order), y=factor(order_superfam, levels=te_order), size=count, color=count)) + 
  geom_point(alpha = 0.8) + scale_color_gradient2(low="#F7FCFD", mid="#BFD3E6", high="#88419D") + xlab("individual") + ylab("superfamily") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# heatmap - this one is better
ggplot(top_20_TE_df, aes(x=factor(individual, levels=volc_isl_order), y=factor(order_superfam, levels=te_order), fill=count)) + 
  geom_tile(color = "white", lwd = 0.5, linetype = 1) + scale_fill_gradient2(low="#F7FCFD", mid="#E0ECF4", high="#88419D") + xlab("individual") + ylab("superfamily") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### PERCENTAGE OF TE SUPERFAMS ###
dir <- "/Users/heidiyang/tetragnatha_project/dnaPipeTE_counts_data/Counts" # this is the counts data!
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
counts_df[,-1] <- apply(counts_df[,-1],2,as.numeric)
genomic_repeat_content <- apply(counts_df[,-1], 2, function(x) {sum(x[1:9])})

# get the proportion of DNA hAT to Total
just_hats <- land_fam_df %>%
  filter(order_superfam == "DNA/hAT")
prop_hats <- as.numeric(just_hats$n) / genomic_repeat_content
just_hats["prop_hats"] <- prop_hats
just_hats["perc_hats"] <- prop_hats * 100
just_hats_sort <- just_hats %>%
  arrange(desc(perc_hats))

# test for significant differences in hAT transposons
just_hats_vol_df <- land_fam_df %>%
  filter(order_superfam == "DNA/hAT") %>%
  group_by(volcano) %>%
  summarise(count = n())
kruskal.test(count ~ volcano, data = just_hats_vol_df)

ggplot(aes(x=factor(volcano, levels=vol_order), y=count, fill=factor(volcano, levels=vol_order)),data=just_hats_vol_df) + geom_bar(stat="identity") +
  scale_fill_brewer(palette="BuPu", direction=-1, breaks=vol_order) + labs(fill = "volcano") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("volcano")




