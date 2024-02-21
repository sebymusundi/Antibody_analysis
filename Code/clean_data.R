rm(list=ls())

# load pacakges 
library(tidyverse)
library(readxl)
library(FactoMineR)
library(factoextra)
library(EnvStats)
library(ggpubr)
library(ggrepel)
library(heatmaply)
library(heatmap3)
library(data.table)
library(viridis)
library(RColorBrewer)
library(ggpval)  
library(corrplot)
library(psych)
library(heatmap3)
library(gplots)
library(heatmaply)
library(viridis)

# Load data

## antibody reactivity data
az_proteins <-  read.csv("Data/az_Proteinswithmetadata.2023.csv")

## protein families data 
protein_families <- read_tsv("Data/no.antigens.participants.2023.tsv", 
                             show_col_types = FALSE)

# Convert ESPs to become BSPs

protein_families$Family[protein_families$Family=="ESP"] = "BSP"


# make colum names for protein families unique 
colnames(protein_families) <- make.names(colnames(protein_families), 
                                         unique = T)


## Make protein names unique

protein_families$GeneID.and.Domain <- make.names(protein_families$GeneID.and.Domain,
                                                 unique = T)

## convert to data frame 

protein_families <- as.data.frame(protein_families)

