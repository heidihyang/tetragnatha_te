# Repeat Landscape Plot
# from Clement's Landscapes.R script

# load libraries and setwd
library(ggplot2)
library(stringr)
library(viridis)
library(data.table)
library(reshape2)
library(RColorBrewer)
library(dplyr)
# needs to have reads_landscape and blast_reads.counts in it
dir <- "/Users/heidiyang/tetragnatha_project/landscape"
setwd(dir)

# get landscape and count files
land_files <- list.files(path=dir, pattern="*_landscape", 
                         full.names=FALSE, recursive=FALSE)
count_files <- list.files(path=dir, pattern="*.counts", 
                          full.names=FALSE, recursive=FALSE)
all_files <- list(land_files, count_files)
file_df <- as.data.frame(do.call(cbind, all_files))
file_df$specimen <- substr(file_df$V2, 1, nchar(file_df$V2)-21)
rownames(file_df) <- NULL

## different graphs ##
# histogram
ggplot(land, aes(div)) +
  geom_histogram(binwidth=1.1) +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  geom_density(aes(y=1.1*..count..)) +
  #ggtitle(specimen_name) +
  xlab("blastn divergence") +
  ylab("proportion of the genome") +
  scale_x_continuous(limits = c(0, 35)) +
  scale_y_continuous(labels=function(x)signif(x/as.numeric(count)*100, 3), limits=c(0,25000)) +
  theme(legend.position="none")

# repeat landscape histogram plot for each individual
i <- 2
for (i in 1:nrow(file_df)) {
  row <- file_df[i,]
  land_name <- as.character(row[1])
  land <- read.table(land_name)
  count_name <- as.character(row[2])
  count <- read.table(count_name)
  specimen_name <- substr(count_name, 1, nchar(count_name)-21)
  # name columns
  names(land) <- c("id", "name", "fam1", "fam")
  # calculate divergence
  land$div <- 100-land$id
  # split fam1 column to just get superfamily
  land[c('order', 'superfamily')] <- str_split_fixed(land$fam1, '/', 2)
  # remove redundant rows and rows with ?
  land <- land[-c(4)]
  land <- land[(str_sub(land$order, -1)) != "?", ]
  land <- land %>%
    filter(order != "Unknown")
  #land <- land[land$div <= 5 & land$family == 'hAT']
  
  #plot the landscape graph
  lp <- ggplot(land, aes(div, fill=order)) +
    geom_histogram(binwidth=1.1) +
    scale_fill_viridis(discrete = T, name="TE order") +
    theme_bw() +
    ggtitle(specimen_name) +
    xlab("blastn divergence") +
    ylab("% of genome") +
    scale_x_continuous(limits = c(0, 35)) +
    scale_y_continuous(labels=function(x)signif(x/as.numeric(count)*100, 3), limits=c(0,25000))
  
  # save graph
  ggsave(filename=paste(specimen_name, 'landscape_plot.png', sep="_"), plot=lp, height=6, width=8)
}

### repeat landscape density plot ###
# get divergence for all species
div_species <- function(row) {
  land_name <- as.character(row[1])
  land <- read.table(land_name)
  count_name <- as.character(row[2])
  count <- read.table(count_name)
  specimen <- row[3]
  sep <- strsplit(as.character(specimen), "_")
  spec <- paste(sep[[1]][c(1,2,6)],collapse="_")
  volcano <- paste(sep[[1]][3:4],collapse="_")
  # name columns
  colnames(land) <- c("id", "name", "fam1", "fam")
  # calculate divergence
  land$div <- 100-land[,1]
  # split fam1 column to just get superfamily
  land[c('superfamily', 'family')] <- str_split_fixed(land$fam1, '/', 2)
  # remove redundant rows and rows with ?
  land <- land[-c(4)]
  land <- land[(str_sub(land$superfamily, -1)) != "?", ]
  return(c(spec, volcano, land$div))
}

# apply to file_df
land_div_values <- apply(file_df, 1, function(x) {div_species(x)})
land_div_df <- data.frame(do.call(rbind, land_div_values))
colnames(land_div_df)[1] <- "specimen"
colnames(land_div_df)[2] <- "volcano"
# replace name of outgroups
land_div_df$volcano <- replace(land_div_df$volcano, c(1,8,14,25), "web-builders")

# make melted dataframe
land_div_melt <- melt(land_div_df, id=c("specimen", "volcano"))[c(1,2,4)]
div_df <- land_div_melt[order(land_div_melt$specimen),]
div_df$volcano <- factor(div_df$volcano, levels = c("Kauai_Ka", "Oahu_Koo", "Oahu_Wai", "Molokai_Mo", "Maui_EM", "Hawaii_Koh", "Hawaii_ML", "Hawaii_Kil", "web-builders"))
div_df$value <- as.numeric(div_df$value)
div_df <- na.omit(div_df)

# graph!
ggplot(aes(x=value, group=specimen, color=volcano), data=div_df)+ geom_density() + 
  labs(x = "blastn divergence", y = "Proportion of genome") + 
  scale_color_brewer(palette="BuPu", direction = -1, breaks = c("Kauai_Ka", "Oahu_Koo", "Oahu_Wai", "Molokai_Mo", "Maui_EM", "Hawaii_Koh", "Hawaii_ML", "Hawaii_Kil", "web-builders")) +
  theme_classic()
