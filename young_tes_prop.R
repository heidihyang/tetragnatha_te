# load libraries and setwd
library(ggplot2)
library(stringr)
library(viridis)
library(data.table)
library(dplyr)
library(RColorBrewer)
# needs to have reads_landscape and blast_reads.counts in it
dir <- "/Users/heidiyang/tetragnatha_project/landscape/"
setwd(dir)

# get landscape and count files
land_files <- list.files(path=dir, pattern="*_landscape", full.names=FALSE, recursive=FALSE)
count_files <- list.files(path=dir, pattern="*.counts", full.names=FALSE, recursive=FALSE)
all_files <- list(land_files, count_files)
file_df <- as.data.frame(do.call(cbind, all_files))

# for loop for young te proportions
young_dudes <- list()
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
  land[c('superfamily', 'family')] <- str_split_fixed(land$fam1, '/', 2)
  # remove redundant rows and rows with ?
  land <- land[-c(4)]
  land <- land[(str_sub(land$superfamily, -1)) != "?", ]
  # get the young TEs less than 5% divergent from TE contig
  land_young <- land[land$div < 5.00,]
  # proportion of young TEs
  young_prop <- nrow(land_young)/nrow(land) * 100
  # add to dataframe
  young_dudes[[i]] <- c(specimen_name, young_prop)
}

# fix dataframe
young_props <- data.frame(do.call(rbind, young_dudes))
young_props$X2 <- as.numeric(young_props$X2)
colnames(young_props) <- c("specimen", "young_prop")
spec_vol_names <- young_props[,-2]

# get columns in right format and rows in right order
young_props$volcano <- sapply(young_props$specimen, function(x) {paste(strsplit(x,"_")[[1]][3:4], collapse="_")})
young_props$specimen <- sapply(young_props$specimen, function(x) {paste(strsplit(x,"_")[[1]][c(1,2,6)], collapse="_")})
young_props$volcano <- replace(young_props$volcano, c(1,8,14,25), "web-builders")

# fix specimen order
vol_order <- c("Kauai_Ka", "Oahu_Wai", "Oahu_Koo", "Molokai_Mo", "Maui_EM", "Hawaii_Koh", "Hawaii_ML", "Hawaii_Kil", "web-builders")
young_props <- young_props%>%
  mutate(volcano =  factor(volcano, levels = vol_order)) %>%
  arrange(volcano)
spec_order <- young_props$specimen
young_props$volcano <- as.character(young_props$volcano)

young_props$volcano[young_props$volcano == "web-builders"] <- "web-builders_WB"
species <- str_split_i(young_props$specimen, pattern="_", 2)
volc <- data.frame(str_split_fixed(young_props$volcano, pattern="_", 2))
volc["specimen"] <- young_props$specimen
volc["island"] <- sapply(volc$X1, function(x) substr(x, start = 1, stop = 2))
volc["volc_isl"] <- paste("(", volc$island, "-", volc$X2, ")", sep="")
volc <- volc %>% 
  mutate(volc_isl=ifelse(grepl("we", volc_isl, fixed=TRUE),
                         gsub("we-","",volc_isl),volc_isl)) %>%
  mutate(volc_isl=ifelse(grepl("Mo", volc_isl, fixed=TRUE),
                         gsub("Mo-","",volc_isl),volc_isl)) %>%
  mutate(volc_isl=ifelse(grepl("Ka", volc_isl, fixed=TRUE),
                         gsub("Ka-","",volc_isl),volc_isl))
young_props$volc_isl <- paste("T.", species, volc$volc_isl, sep=" ")
young_props <- young_props %>%
  mutate(volc_isl=ifelse(grepl("filiciphila", volc_isl, fixed=TRUE),
                         gsub("filiciphila","filiciphilia",volc_isl),volc_isl)) %>%
  mutate(volc_isl=ifelse(grepl("stella", volc_isl, fixed=TRUE),
                         gsub("stella","stela",volc_isl),volc_isl)) %>%
  mutate(volcano=ifelse(grepl("web-builders_WB", volcano, fixed=TRUE),
                        gsub("_WB","",volcano),volcano)) %>%
  mutate(volcano =  factor(volcano, levels = vol_order)) %>%
  arrange(volcano)
volc_isl_order <- young_props$volc_isl


# plot!
ggplot(young_props, aes(x=factor(volc_isl, levels=volc_isl_order), 
                        y=young_prop, fill=volcano)) + geom_bar(stat="identity") + 
  labs(x = "individual", y = "% young TEs") + 
  scale_fill_brewer(palette="BuPu", direction = -1, breaks = 
                      c("Kauai_Ka", "Oahu_Wai", "Oahu_Koo", "Molokai_Mo", 
                        "Maui_EM", "Hawaii_Koh", "Hawaii_ML", "Hawaii_Kil", 
                        "web-builders")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(young_props, aes(x=volcano, y=young_prop, fill=volcano)) + geom_boxplot() + 
  labs(x = "volcano", y = "% young TEs") + 
  scale_fill_brewer(palette="BuPu", direction = -1, breaks = 
                      c("Kauai_Ka", "Oahu_Wai", "Oahu_Koo", "Molokai_Mo", 
                        "Maui_EM", "Hawaii_Koh", "Hawaii_ML", "Hawaii_Kil", 
                        "web-builders")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### UNNECESSARY STUFF ###
# get averages of young TEs for radiating and outgroup individuals - don't run again!
young_rad <- young_props[!(young_props$X1 =="T_acuta_Maui_EM_lightbrown_071" | young_props$X1 == 'T_filiciphila_Maui_EM_outgroup_074' | young_props$X1 == 'T_maka_Kauai_Ka_lightbrown_052' | young_props$X1 == 'T_stellarobusta_Maui_EM_maroon_001'),]
young_rad_mean <- mean(young_rad$X2)

young_sis <- young_props[(young_props$X1 =="T_acuta_Maui_EM_lightbrown_071" | young_props$X1 == 'T_filiciphila_Maui_EM_outgroup_074' | young_props$X1 == 'T_maka_Kauai_Ka_lightbrown_052' | young_props$X1 == 'T_stellarobusta_Maui_EM_maroon_001'),]
young_sis_mean <- mean(young_out$X2)

# idk what I was doing here but this is a thing
mr_df <- young_props[!(young_props$specimen =="T_acuta_Maui_EM_lightbrown_071" | young_props$specimen == 'T_filiciphila_Maui_EM_outgroup_074' | young_props$specimen == 'T_maka_Kauai_Ka_lightbrown_052' | young_props$specimen == 'T_stellarobusta_Maui_EM_maroon_001'),]
mr_df$spiny_or_web <- 'spiny_leg'
mn_df <- young_props[(young_props$specimen =="T_acuta_Maui_EM_lightbrown_071" | young_props$specimen == 'T_filiciphila_Maui_EM_outgroup_074' | young_props$specimen == 'T_maka_Kauai_Ka_lightbrown_052' | young_props$specimen == 'T_stellarobusta_Maui_EM_maroon_001'),]
mn_df$spiny_or_web <- 'web_builder'
all_df <- rbind(mr_df, mn_df)
all_df$young_prop <- as.numeric(all_df$young_prop)

my_df <- melt(setDT(mr_df), id=c("island", 'volcano', 'ecomorph', 'specimen', 'spiny_or_web'))

# da stats time
yw <- kruskal.test(value ~ unlist(volcano), data = my_df)
print(yw)

yp <- ggplot(young_props, aes(specimen, young_prop)) + geom_bar(stat="identity", fill="lightblue") + xlab('specimen') + ylab('% young TEs') + 
  theme_bw() + theme(axis.text.x = element_text(angle=90))
