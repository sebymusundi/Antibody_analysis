##############################################################################################################
# Plot volcano plot to compare antibody responses between primigravida and 
# multigravida 
################################################################################
# load packages 
library(tidyverse)
library(ggrepel)

# Clean up data to only use multigravida 2 and 3 

az_proteins_volcano <- az_proteins%>%
  filter(Gravida %in% c(1,2,3))



# Classify gravidity into multigravid and primigravid 

az_proteins_volcano$Gravida[az_proteins_volcano$Gravida==1]="Primigravida"
az_proteins_volcano$Gravida[az_proteins_volcano$Gravida==2]="Multigravida"
az_proteins_volcano$Gravida[az_proteins_volcano$Gravida==3]="Multigravida"


# make Gravida a factor

az_proteins_volcano <- az_proteins_volcano%>%
  mutate(Gravida=factor(Gravida, 
                        levels = c("Primigravida", "Multigravida")))

# Identify proteins whose seropositivity was less than 10%

seroposite_proteins <- data_new %>%
  filter(seroprevalence>10)

# Maake the dataframe longer

full_data <- az_proteins_volcano%>%
  pivot_longer(cols = PF3D7_0100100_DBL.0.11:PF3D7_1116100)

# select only the seroreactive proteins for downstream analysis

full_data <- full_data%>%
  filter(name %in% c(seroposite_proteins$Protein_IDs))


# select primagravida data

primigravida <- full_data%>%
  filter(Gravida=="Primigravida")%>%
  dplyr::select(Study.Number, name, value)%>%
  pivot_wider(names_from = Study.Number, 
              values_from = value) %>%
  as.data.frame()

# Ensure data is in numeric form

primigravida[,2:13] <- as.data.frame(sapply(primigravida[,2:13], as.numeric))


# Assign row names to primigravida 

rownames(primigravida) <- primigravida$name


# Remove the first name

primigravida <- primigravida[-1]


# Assign all rows to become primigravida 

colnames(primigravida)[1:12] <- "primigravida"


# Check if all values are greater than background 

primagravida.check <-matrix(data = NA, ncol = 1, nrow = 678)

colnames(primagravida.check)="sum"

for (i in 1:678){
  primagravida.check[i,] <- sum(primigravida[i,]>1.562500)
}



# Select multigravidas

multigravida <- full_data%>%
  filter(Gravida=="Multigravida")%>%
  dplyr::select(Study.Number, name, value)%>%
  pivot_wider(names_from = Study.Number, 
              values_from = value) %>%
  as.data.frame()



multigravida[,2:31] <- as.data.frame(sapply(multigravida[,2:31], as.numeric))

# Assign rownames to multigravida set

rownames(multigravida) <- multigravida$name


# Remove first row 

multigravida <- multigravida[-1]

# Assign all rows to become multigravida 

colnames(multigravida)[1:30] <- "multigravida"

# Check reactivity 


multigravida.check <-matrix(data = NA, ncol = 1, nrow = 678)

colnames(multigravida.check)="sum"

for (i in 1:678){
  multigravida.check[i,] <- sum(multigravida[i,]>1.562500)
}



# Combine two dataset

combined.dataset <- cbind(primigravida,multigravida)


# Create a matrix

raw_values=matrix(data=NA, nrow = 678, ncol = 1)

colnames(raw_values) <- "P-values"

# Determine raw p-values 

for (i in 1:678){
  raw_values[i] <- t.test(combined.dataset[i,1:12], 
                          combined.dataset[i,13:42], 
                          var.equal=T)$p.value
}



# Find log of all values 

combined.dataset <- log(combined.dataset)

# Find the log mean of each variable 
mean.primagravida <- apply(combined.dataset[1:12],1, geoMean)
mean.multigravida <- apply(combined.dataset[13:42],1, geoMean)

# Calculate the log fold change 

fold.change <- mean.primagravida-mean.multigravida 

# make a dataframe containing log fold change values and p-values 
final.dataset <- as.data.frame(cbind(fold.change,raw_values))

# Indicate protein name as a column 
final.dataset$Protein_ID <- rownames(final.dataset)

# Remove rownames
rownames(final.dataset) <- NULL

# Check column names
colnames(final.dataset)[2] <- "raw_pvalues"

# classify p-values 
final.dataset$Significance[final.dataset$raw_pvalues < 0.05 & final.dataset$fold.change > 0 ] <- "p-values < 0.05 & fold change > 0"
final.dataset$Significance[final.dataset$raw_pvalues < 0.05 & final.dataset$fold.change < 0 ] <- "p-values < 0.05 & fold change < 0"
final.dataset$Significance[final.dataset$raw_pvalues > 0.05] <- "p-values > 0.05"

# volcano plot
png("volcanoplot.png", height = 3000, width =4800, res = 300 )
x_limits <- c(0.25, 0.5)
y_limits <- c(1.3, NA)

x_lim_2 <- c(-0.25, -0.5)
y_lim_2 <- c(1.3, NA)

ggplot(final.dataset, aes(x=fold.change, 
                          y=-log10(raw_pvalues), col=Significance))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red")+
  xlab("log 2 fold change")+
  ylab("- log 10(p-values)")+
  theme_classic()+
  scale_color_manual(values=c('black', "red", "green"))+
  theme(legend.position= "right" )+
  geom_label_repel(aes(label=Protein_ID), 
                   subset(final.dataset, raw_pvalues < 0.05 & fold.change> 0.0), size=3.0, label.size = 0.4, 
                   box.padding = 0, direction = "y", hjust=1,  xlim = x_limits, 
                   ylim = y_limits)+
  geom_label_repel(aes(label=Protein_ID), 
                   subset(final.dataset, raw_pvalues < 0.05 & fold.change < 0.0), size=3.0, label.size = 0.4, 
                   box.padding = 0, direction = "y", hjust=1,  xlim = x_lim_2, 
                   ylim = y_lim_2)+
  ggtitle("Volcano Plot: Primigravida vs Multigravida Antibody Responses Differential Analysis")+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))+
  theme(axis.title.y = element_text(face = "bold", vjust = +2.5, size = 12),
        axis.title.x = element_text(face = "bold", vjust = -0.75, size = 12))+
  theme(axis.text = element_text(colour = "black", size = 10))

dev.off()
