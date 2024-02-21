###############################################################################
#
###############################################################################

rm(list = ls())



########################################################################
# Determining seroprevalance based on protein families
################################################################################

# make new matrix

data_new <- matrix(data = NA, nrow = 698, ncol = 5)

# assign matrix Column names

colnames(data_new) <- c(" antigen_number", "participants_no","seroprevalence", 
                        "Family", "Proteins")


# Loop for determining number of proteins with value greater than threshold across participants 

for (i in 1:698){
  data_new[i] <- sum(protein_families[i,5:57]>1.5625)
}

# Fill second column of new dataframe  

data_new[,2]=53

# Determine seroprevalence

data_new[,3]=(data_new[,1]/data_new[,2])*100

# include protein families in the new dataframe 

data_new[,4] <- protein_families$Family

data_new[,5] <- protein_families$GeneID.and.Domain

# Convert to dataframe 

data_new <- as.data.frame(data_new)

# Convert variables to become factors or numerials

data_new$Family <- as.factor(data_new$Family)

data_new$seroprevalence <- as.numeric(data_new$seroprevalence)

data_new$Protein_IDs <- protein_families$GeneID.and.Domain


# Make a plot for the seroprevalence 

seroplot <- ggplot(data_new, aes(Family, seroprevalence))+
  geom_boxplot(fill="grey90")+
  theme_classic()+
  geom_dotplot(binaxis = "y", stackdir = "center", stackratio = 1.5, 
               dotsize = 1.4,binwidth = 0.5)+
  geom_hline(yintercept = 10, linetype="dashed", colour="red")+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.title=element_text(colour = "black", face="bold", size = 11))+
  xlab(" Protein groups")+
  ylab("Seroprevalence (%)")+
  scale_y_continuous(breaks = seq(0,100,10))


#seroplot


seroprevalence_plot <- ggboxplot(data_new,
                                 x="Family", 
                                 y="seroprevalence", 
                                 add = "jitter", 
                                 color = "Family", 
                                 palette = "Dark2")+
  geom_hline(yintercept = 10, linetype="dashed", colour="red")+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.title=element_text(colour = "black", face="bold", size = 11))+
  xlab(" Protein groups")+
  ylab("Seroprevalence (%)")+
  scale_y_continuous(breaks = seq(0,100,10))+
  theme(legend.position = "none")
seroprevalence_plot

# Find the individual number of each protein family used 
protein_family_numbers <- data_new %>%
  group_by(Family) %>%
  summarise(count=n())

#figEdits
#sero1 = seroplot + 
# theme(axis.text = element_text(size = 12))
#sero1                  

#sero2=sero1 + 
# theme(axis.title = element_text(size = 15))
#sero2

#sero3 = sero2 + 
# theme(axis.title = element_text(face="bold"))
#sero3

#sero4=sero3 + theme(
# axis.title.y = element_text(vjust = +2.0),
#axis.title.x = element_text(vjust = -0.75)) 
##sero4

#sero5= sero4 + theme(panel.grid.major = element_blank(), 
#               panel.grid.minor = element_blank(),
#              panel.background = element_blank(), 
#            axis.line = element_line(colour = "black"))

#sero5

# Assign custom color
#use this if you indicate at aes color = family
#sero6=sero5+scale_color_manual(values=c("black",
# "black",
# "blue",
#"brown",
# "grey",
#"purple",
# "#00AFBB"))
#sero6

#ggsave("seroplot.tiff", height = 6, width = 8, dpi = 300)
ggsave("Results/seroprevalence_plot.tiff", height = 6, width = 8, dpi = 300)

################################################################################
# Relationship between gravidity and Antibody Breadth
################################################################################
# Select specific rows and columns 

az_proteins <- az_proteins[5:57,]

rownames(az_proteins) <- NULL
# convert columns to numeric form 

az_proteins[,33:730] <- as.numeric(unlist(az_proteins[,33:730]))

# assign protein names to columns 

colnames(az_proteins)[33:730] <- protein_families$GeneID.and.Domain

# Matrix of participants 

gravidity_matrix=matrix(data=NA, nrow = 53, ncol =2 )

# Assign column names 

colnames(gravidity_matrix) <- c("Gravidity", "Antigen_number")

# Loop 

for(i in 1:53){
  gravidity_matrix[,1] <- az_proteins[1:53,4]
  gravidity_matrix[i,2] <- sum(az_proteins[i,33:730]>1.5625)
}


# Convert to become dataframe 

gravidity_matrix <- as.data.frame(gravidity_matrix)

# Convert gravidity to become a factor

gravidity_matrix$Gravidity <- as.factor(gravidity_matrix$Gravidity)

# Plot a curve for gravidity versus number of antigens 

gravidity.vs.antibody_breadth <- ggplot(gravidity_matrix, aes(Gravidity, Antigen_number))+
  geom_boxplot()+
  geom_jitter(width = 0.05)+
  theme_classic()+
  ylab("Antibody breadth")+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.title = element_text(colour = "black", face = "bold", size = 12))+
  scale_y_continuous(labels=seq(0,700,100), breaks = seq(0,700,100))

#gravidity.vs.antibody_breadth

# Calculating antibody breadth and plotting p-values

my_comparisons <- list(c("1", "2"), c("1","3"), c("2","3"))
gravidity_vs_antibody_breadth <- ggboxplot(gravidity_matrix, 
                                           x="Gravidity", 
                                           y="Antigen_number", 
                                           xlab="Gravidity", 
                                           ylab = "Antibody breadth", 
                                           add = "jitter", 
                                           color="Gravidity", 
                                           palette = "Dark2")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(exact=FALSE))+
  stat_compare_means(label.y = 900) +
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.title = element_text(colour = "black", face = "bold", size = 11))+
  scale_y_continuous(labels= seq(0,900,100), 
                     breaks = seq(0,900,100), 
                     limits = c(0,900))


gravidity_vs_antibody_breadth
ggsave("Results/gravidity_vs_antibody_breadth.tiff", height = 6, width = 8, dpi = 300)
#figEdits
#version1 = gravidity.vs.antibody_breadth + 
#  theme(axis.text = element_text(size = 12))
#version1                  

#version2=version1 + 
#  theme(axis.title = element_text(size = 15))
#version2

#version3 = version2 + 
# theme(axis.title = element_text(face="bold"))
#version3

#Version4=version3 + theme(
#  axis.title.y = element_text(vjust = +2.0),
#  axis.title.x = element_text(vjust = -0.75)) 
#Version4


#version5= Version4 + theme(panel.grid.major = element_blank(), 
#                           panel.grid.minor = element_blank(),
#                           panel.background = element_blank(), 
#                           axis.line = element_line(colour = "black"))
#version5

# Assign custom color
#to use if you add color =gravida on aes
#version6=version5+scale_color_manual(values=c("blue",
#"blue",
#"blue",
# "red",
# "red"))

#version6

# Assign brewer color
#ch+scale_color_brewer(palette="Dark2")

# Assign gray scale
#ch+scale_color_grey() + theme_classic()


#ggsave("gravidity.vs.antibody_breadth.tiff", height = 6, width = 8, dpi = 300)

###############################################################################
# Scatterplot for antibody breadth versus age of participants 
###############################################################################

# Make matrix 
age_matrix <- matrix(data = NA, nrow = 53, ncol = 2)

# Assign column names to matrix 

colnames(age_matrix) <- c("Age", "Antibody_breadth")

# enter data

age_matrix[,1] <- az_proteins[,3]
age_matrix[,2] <- gravidity_matrix[,2]


# Convert to data frame
age_matrix <- as.data.frame(age_matrix)


# Plot scatter plot for age versus antibody breadth

antibody_breadth.vs.age <- ggplot(age_matrix, aes(Age, Antibody_breadth))+
  geom_point()+
  geom_smooth(method = "lm", color="blue")+
  theme_classic()+
  ylab("Antibody breadth")+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.title = element_text(colour = "black", face = "bold", size = 12))+
  stat_regline_equation(label.x =34)+
  stat_cor(label.x = 34, label.y = 580)


antibody_breadth.vs.age

#figEdits

#vers1 = antibody_breadth.vs.age + 
#  theme(axis.text = element_text(size = 12, face = "bold", color = "black"))
#vers1                  

#vers2=vers1 + 
#  theme(axis.title = element_text(size = 15))
#vers2

#Vers3 = vers2 + 
#  theme(axis.title = element_text(face="bold"))
#Vers3

#Vers4=Vers3 + theme(
#  axis.title.y = element_text(vjust = +2.0),
#  axis.title.x = element_text(vjust = -0.75)
#) 

#Vers4

ggsave("Results/antibody_breadth.vs.age.tiff", height = 6, width = 8, dpi = 300)


###############################################################################
# Scatterplot for antibody breadth versus Hemoglobin of participants 
########################################################################
# Make matrix 
hb_matrix <- matrix(data = NA, nrow = 53, ncol = 2)

# Assign column names to matrix 

colnames(hb_matrix) <- c("HB", "Antibody_breadth")

# enter data

hb_matrix[,1] <- az_proteins[,18]
hb_matrix[,2] <- gravidity_matrix[,2]
hb_matrix


# Convert to data frame
hb_matrix <- as.data.frame(hb_matrix)


antibody_breadth.vs.hb <- ggplot(hb_matrix, aes(HB, 
                                                Antibody_breadth))+geom_point(colour = "black", size = 2)+ 
  geom_smooth(method="lm", col="blue")+
  stat_regline_equation(label.x = 13, label.y = 640)+
  stat_cor(method = "spearman", label.x = 13, 
           label.y = 610)+
  xlab("HB [g/dl]") +
  ylab('Antibody breadth')+
  theme_classic()+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.title = element_text(colour = "black", face = "bold", size = 12))

antibody_breadth.vs.hb

#versi1 = antibody_breadth.vs.hb + 
#  theme(axis.text = element_text(size = 12, face = "bold", color = "black"))
#versi1                  
#versi2=versi1 + 
#  theme(axis.title = element_text(size = 15))
#versi2
#versi3 = versi2 + 
# theme(axis.title = element_text(face="bold"))
#versi3
#versi4=versi3 + theme(
#  axis.title.y = element_text(vjust = +2.0),
# axis.title.x = element_text(vjust = -0.75)) 
#versi4

ggsave("Results/antibody_breadth.vs.hb.tiff", height = 6, width = 8, dpi = 300)




##############COMBINED FIGURE 2

#library(cowplot)
#combined.figure2 <- plot_grid (Vers4,version5,versi4, 
#                               nrow = 2, labels = c("A", "B", "C", 
#                                                   label_size = 12))

#combined.figure2


#ggsave("combined Fig 2.tiff", height = 10, width = 18)

###################################################################################################################
# Heatmap for all proteins 
##################################################################################################################

library(ggbeeswarm)

# Extract data for primigravid 1,2  and 3
az_proteins_heatmap <- az_proteins%>%
  filter(Gravida %in% c(1,2,3))

# Classify multigravid 2 and 3 together
az_proteins_heatmap$Gravida[az_proteins_heatmap$Gravida==1] ="Primigravid"
az_proteins_heatmap$Gravida[az_proteins_heatmap$Gravida==2] ="Multigravida"
az_proteins_heatmap$Gravida[az_proteins_heatmap$Gravida==3] ="Multigravida"



# summary of statistics 
summary_stats <- az_proteins_heatmap %>%
  select(Age....yrs., Birth.Weight.of.Newborn, HB..g.dl., Gravida) %>%
  pivot_longer(cols = 1:3)%>%
  mutate(Gravida=as.factor(Gravida))



# extract birthweight information
summary_stats_bw <- summary_stats %>%
  filter(name== "Birth.Weight.of.Newborn") %>%
  mutate(Gravida=as.factor(Gravida)) %>%
ggplot()+
  geom_beeswarm(aes(Gravida, value, color= Gravida), cex = 4, size=2)+
  theme_classic()+
  ylab("Birth weight(Kg)")+
  ggtitle("Birth weight vs Gravida")+
  theme(axis.text.x = element_text(size = 10, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title =element_text(colour="black", face = "bold"), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5))+
  scale_color_brewer(palette = "Set1")+
  scale_y_continuous(labels =seq(2.4,4.0,0.4),
                     breaks = seq(2.4,4.0,0.4),
                     limits = c(2.4,4.0))
  
summary_stats_bw

# Compare hemoglobin levels versus gravida 

summary_stats_hb <- summary_stats %>%
  filter(name== "HB..g.dl.") %>%
  mutate(Gravida=as.factor(Gravida)) %>%
  ggplot()+
  geom_beeswarm(aes(Gravida, value, color= Gravida), cex = 4, size=2,
                priority = "density")+
  theme_classic()+
  ylab("Haemoglobin level (g/dl)")+
  ggtitle("Haemoglobin vs Gravida")+
  theme(axis.text.x = element_text(size = 10, colour = "black"), 
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title =element_text(colour="black", face = "bold"), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5))+
  scale_color_brewer(palette = "Set1")+
  scale_y_continuous(labels =seq(10,18,2),
                     breaks = seq(10,18,2),
                     limits = c(10,18))
summary_stats_hb

library(cowplot)

summary_plots <- plot_grid(summary_stats_bw, summary_stats_hb)

ggsave("Results/summary_plots.tiff", height = 6, width = 10, dpi = 300)

# Make the table longer 
az_proteins_heatmap <-  az_proteins_heatmap %>% 
  pivot_longer(cols = 33:730, 
               names_to = "Protein_ID", 
               values_to = "Ab_responses") %>% 
  arrange(Protein_ID)

# Extract protein family groups for all gravida selected 
# Ensure all proteins are listed as factors 
#az_proteins_heatmap$Protein_ID <- as.factor(az_proteins_heatmap$Protein_ID)

# Check number of proteins
#levels(az_proteins_heatmap$Protein_ID)

#check structure of data
#str(az_proteins_heatmap)

# sort data based on protein ID
#az_proteins_heatmap$Protein_ID <- sort(az_proteins_heatmap$Protein_ID, decreasing = FALSE)
#protein_families$GeneID.and.Domain <- sort(protein_families$GeneID.and.Domain, decreasing = FALSE)



# Repeat the protein IDs 23916 (42 samples * 698 proteins) 

protein_class <- as.data.frame(rep(protein_families[,4], 42))
families <- as.data.frame(rep(protein_families[,2], 42))
# Assign column names 

colnames(protein_class) <- "protein_class"
colnames(families) <- "protein_family"

# sort protein class to have a specific order
#protein_class$protein_class <- sort(protein_class$protein_class)


# write loop for extracting values 
ordered_protein_class<- order(protein_class$protein_class)

# extract protein families 
extracted_protein_families <- as.data.frame(families$protein_family[ordered_protein_class])

# assign required column names 
colnames(extracted_protein_families) <- "protein_family"

# join columns 
az_proteins_heatmap <- cbind(az_proteins_heatmap, extracted_protein_families)

# sort data based on gravidity
az_proteins_heatmap <-  az_proteins_heatmap %>%
  arrange(desc(Gravida))


# Extract heatmap data for primigravid 
data_heatmap <- az_proteins_heatmap %>%
  select(Study.Number, Protein_ID, Ab_responses) %>%
  pivot_wider(values_from = Ab_responses, names_from = Study.Number) %>%
  as.data.frame()

# 

# Assign the first column as a rownames

rownames(data_heatmap) <-data_heatmap[,1]

# Extract data for heatmap plot

data_heatmap <- data_heatmap[, 2:43]

# Check the structure of the data 
str(data_heatmap)

# transformed heatmap
t_heatmap <-  t(data_heatmap)

t_heatmap <- as.data.frame(t_heatmap)

# extract the protein families protein names
protein_fam <- protein_families[,c(2,4)]

# sort order of protein families
protein_fam <- protein_fam %>%
  arrange(Family)

# create order 
protein_order <- as.character(protein_fam$GeneID.and.Domain)

# order transformed heatmap based on the arranged protein families 
t_heatmap <- t_heatmap %>%
  select(all_of(protein_order))

# Create colour codes for protein families 
color_vector_1 <- rep("red" , 158)
color_vector_2 <- rep("blue", 108)
color_vector_3 <- rep("green", 163)
color_vector_4 <- rep("orange", 182)
color_vector_5 <- rep("purple", 54)
color_vector_6 <- rep("pink", 33)

# Concatenate or join the color vectors into a single vector
row_colors_vector <- c(color_vector_1, color_vector_2, color_vector_3, 
                       color_vector_4, color_vector_5, color_vector_6)

# Apply colour codes for  primigravid and multigravid women 

color_vector_7 <- rep("yellow", 12)
color_vector_8 <- rep("grey", 30)

column_color_vector <- c(color_vector_7, color_vector_8)


# plot heatmap
heatmap_col <-  colorRampPalette(c("black", "red", "green"))(100)

png("heatmap.png", height = 3000, width =6000, res = 300)

heatmap.2(as.matrix(t_heatmap), 
          col= heatmap_col,
          dendrogram = "none", 
          scale = "none",
          trace="none", 
          sepcolor="black", 
          key.title="Antibody responses", 
          key.xlab = "value", 
          key.ylab = "",
          keysize = 1.2, 
          Rowv = FALSE, 
          Colv = FALSE, 
          key.par = list(mar = c(4, 1, 1, 5)), 
          sep=c(0,0), 
          cexRow = 0.1, 
          cexCol = 0.1, 
          ColSideColors = row_colors_vector, 
          RowSideColors = column_color_vector, 
          labRow = NA, 
          labCol = NA, 
          margins = c(5,4), 
          xlab = "Antigens", 
          ylab = "Participants")

legend("topright", legend =c("BSP", "CIDR", "DBL", "RIFIN", "STEVOR", "SURFIN", 
                             "Primigravid", "Multigravid" ), 
       fill = c("red", "blue", "green", "orange", "purple", "pink", 'yellow', "grey"), 
       border = FALSE, y.intersp = 0.9, cex = 0.8, ncol = 4)

dev.off()


#################################################################################################################
# Plot volcano plot to compare antibody responses between primigravida and 
# multigravida 
################################################################################

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

logisitic_regression <- az_proteins_volcano[c(4,33:730), ]



#vol.plot
#vol.plot2 = vol.plot+ theme(axis.text=element_text(size=12),
#                           axis.title=element_text(size=14,face="bold"))
#vol.plot2
#ggsave("volcano.plot.tiff", height = 10, width = 18, dpi = 300)

###############################################################################
# lINEAR REGRESSION 
#################################################################################

data_regression <- az_proteins_heatmap %>%
  select(Study.Number, Protein_ID, Ab_responses) %>%
  pivot_wider(values_from = Ab_responses, names_from = Study.Number) %>%
  as.data.frame()






#########################################################################################
#SAVE_ALL_PROTEINS

write.csv(final.dataset, file= "Volcano.all.proteins.csv")

p.all.proteins<- final.dataset[,2]
p.all.proteins
p.all.proteins <- as.numeric(unlist(p.all.proteins))
p.all.proteins

adj_p.all.proteins = p.adjust(p.all.proteins, method = "fdr") # qvalue of mydata

adj_p.all.proteins
write.csv(adj_p.all.proteins, file= "Volcano.all.proteins.adjusted_pvalues.csv")


#SAVE_Significant proteins ONLY
sig.proteins <- final.dataset%>%
  filter(raw_pvalues<0.05)
write.csv(sig.proteins, file= "Volcano.sig.proteins.csv")


#adjust p values for sig proteins
sig.proteins 

p.sig.proteins<-sig.proteins[,2]
p.sig.proteins
p.sig.proteins <- as.numeric(unlist(p.sig.proteins))
p.sig.proteins

adj_p.sig.proteins = p.adjust(p.sig.proteins, method = "fdr") # qvalue of mydata

adj_p.sig.proteins
write.csv(adj_p.sig.proteins, file= "Volcano.sig.proteins.adjusted_pvalues.csv")

######################################
# Select protein that were significant for mutigravida 

data_multigravida <- final.dataset%>%
  filter(fold.change<0.0 & raw_pvalues <0.05)

# filter the proteins from the original dataset

multigravida_data <- full_data%>%
  filter(name %in% c(data_multigravida$Protein_ID))

# Convert to required class

multigravida_data$Gravida <- as.factor(multigravida_data$Gravida)
multigravida_data$value <- as.numeric(multigravida_data$value)

# Plot significant proteins in multi-gravida

multigravida_significant <- ggboxplot(multigravida_data, 
                                      x="Gravida", 
                                      y="value", 
                                      color= "Gravida", 
                                      add="jitter", 
                                      palette = "jco")+
  facet_wrap(~name)+
  ylab("Antibody response")+
  labs(fill="Protein")+
  theme(axis.text = element_text(colour = "black"))+
  stat_compare_means(aes(group=Gravida), method = "t.test")+
  theme(axis.title = element_text(colour = "black", face = "bold"))

multigravida_significant

ggsave("multigravida_significant.tiff", height = 8, width = 10, dpi = 300)

###############################################################################
# Assess overall gravidity of 3 significant multigravida proteins 

# Make the data-frame longer and  filter out the three significant multigravida 
# proteins with their gravidity

multigravida_analysis <- az_proteins %>%
  pivot_longer(cols = 33:730, 
               names_to = "Protein_ID", 
               values_to = "Ab_response") %>%
  filter(Protein_ID %in% c(data_multigravida$Protein_ID))%>%
  select(Study.Number, Gravida, Ab_response, Protein_ID)%>%
  mutate(Gravida=as.factor(Gravida))


# Plot individual graphs assessing gravidity of the three proteins. 

plot_1_gravida <- ggplot(multigravida_analysis, aes(Gravida, Ab_response,
                                                    color=Protein_ID))+
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  theme_classic()+
  ylab("Antibody response")+
  xlab("Gravida")+
  facet_wrap(~Protein_ID)+
  theme(axis.text = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face = "bold"))+
  theme(legend.position = "none")+
  scale_color_brewer(palette="Dark2")

plot_1_gravida

ggsave("plot_1_gravida.tiff", height = 6, width = 8, dpi = 300)

# or 
plot_2_gravida <- ggboxplot(multigravida_analysis, 
                            x="Gravida", 
                            y="Ab_response",
                            palette = "Dark2",
                            color = "Protein_ID", 
                            add = "jitter") +
  facet_wrap(~Protein_ID)+
  theme(axis.text = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black"))+
  theme(legend.position = "none")+
  xlab("Gravida")+
  ylab("Antibody response")

plot_2_gravida  

###############################################################################


# Extract data for significant proteins for primigravida

data_primigrivida <- final.dataset%>%
  filter(fold.change>0.0 & raw_pvalues <0.05)

# Extract heatmap data for primigravida

heatmap_data <- full_data%>%
  filter(name %in% c(data_primigrivida$Protein_ID))%>%
  arrange(Gravida)%>%
  select(Study.Number, name, value, Gravida)%>%
  filter(Gravida=='Primigravida')%>%
  as.data.frame()


# Develop heatmap using geom tiles
Heatmap_primigravida <- ggplot(heatmap_data, aes(Study.Number, name, 
                                                 fill=value))+
  geom_tile(color="grey")+
  ylab("Proteins")+
  xlab("")+
  scale_fill_gradient(low = "white", high = "red")+
  theme(axis.text.y = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face = "bold"))+
  theme(axis.text.x = element_blank())+
  theme(legend.position = "none")+
  ggtitle("Primigravida")+
  theme(plot.title= element_text(hjust = 0.5, size = 12))

Heatmap_primigravida


# Extract heatmap data for primigravida

heatmap_data_multigravida <- full_data%>%
  filter(name %in% c(data_primigrivida$Protein_ID))%>%
  arrange(Gravida)%>%
  select(Study.Number, name, value, Gravida)%>%
  filter(Gravida=="Multigravida")%>%
  as.data.frame()


# Develop heatmap using geom tiles
Heatmap_multigravida <- ggplot(heatmap_data_multigravida, aes(Study.Number, name, 
                                                              fill=value))+
  geom_tile(color="grey")+
  ylab("")+
  xlab("")+
  scale_fill_gradient(low = "white", high = "red")+
  theme(axis.text.y = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face = "bold"))+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank())+
  ggtitle("Multigravida")+
  theme(plot.title= element_text(hjust = 0.5, size = 12))

Heatmap_multigravida


library(cowplot)

plot_grid(Heatmap_primigravida, Heatmap_multigravida)

##############################################################################
# PCA Analysis
###################################################################################

# select data for PCA analysis 

pca_data <- az_proteins_volcano[,33:730]

# Assign row names for PCA data 

rownames(pca_data) <- az_proteins_volcano$Study.Number


# Change entire dataframe to become numeric

pca_data <- as.data.frame(sapply(pca_data, as.numeric))

# Run PCA analysis

res.pca <- PCA(pca_data, scale.unit = TRUE, ncp = 6)

# plot the scree plot for PCA analysis 

pca.screeplot <- fviz_eig(res.pca, addlabels =T, ylim=c(0,25))+
  theme_classic()+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(colour = "black"))

pca.screeplot

# categorize pca plot based on gravidity 

dim1.biplot <- fviz_pca_biplot(res.pca, col.ind= az_proteins_volcano$Gravida, 
                               geom.ind = "point", 
                               label = "ind", 
                               palette=c("red", "green"), 
                               addEllipses = TRUE, 
                               legend.title="Gravidity", 
                               select.var = list(contrib = 5), 
                               pointshape=19, 
                               col.var = "black")+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(colour = "black"))

dim1.biplot

dim2.3.biplot <- fviz_pca_biplot(res.pca, col.ind= az_proteins_volcano$Gravida, 
                                 geom.ind = "point", 
                                 label = "ind", 
                                 palette=c("red", "green"), 
                                 addEllipses = TRUE, 
                                 legend.title="Gravidity", 
                                 select.var = list(contrib = 5), 
                                 pointshape=19, 
                                 col.var = "black", axes = c(2,3))+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(colour = "black"))

dim2.3.biplot


dim3.4biplot <- fviz_pca_biplot(res.pca, col.ind= az_proteins_volcano$Gravida, 
                                geom.ind = "point", 
                                label = "ind", 
                                palette=c("red", "green"), 
                                addEllipses = TRUE, 
                                legend.title="Gravidity", 
                                select.var = list(contrib = 5), 
                                pointshape=19, 
                                col.var = "black",
                                axes = c(3,4))+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(colour = "black"))

dim3.4biplot

dim4.5.biplot <- fviz_pca_biplot(res.pca, col.ind= az_proteins_volcano$Gravida, 
                                 geom.ind = "point", 
                                 label = "ind", 
                                 palette=c("red", "green"), 
                                 addEllipses = TRUE, 
                                 legend.title="Gravidity", 
                                 select.var = list(contrib = 5), 
                                 pointshape=19, 
                                 col.var = "black",
                                 axes = c(4,5))+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(colour = "black"))

dim4.5.biplot

dim5.6.biplot <- fviz_pca_biplot(res.pca, col.ind= az_proteins_volcano$Gravida, 
                                 geom.ind = "point", 
                                 label = "ind", 
                                 palette=c("red", "green"), 
                                 addEllipses = TRUE, 
                                 legend.title="Gravidity", 
                                 select.var = list(contrib = 5), 
                                 pointshape=19, 
                                 col.var = "black",
                                 axes = c(5,6))+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(colour = "black"))

dim5.6.biplot


library(cowplot)
combined.pca <- plot_grid(pca.screeplot, dim1.biplot, 
                          dim2.3.biplot, dim3.4biplot, dim4.5.biplot, 
                          dim5.6.biplot, nrow = 2, labels = c("A", "B", "C", "D", 
                                                              "E", "F"))


combined.pca
ggsave("combined.pca.tiff", height = 12, width = 20)


###################################################################

#VAR2CSA 
VAR2CSA <- full_data%>%
  filter(name %in% c("PF3D7_1200600_DBLpam1",
                     "PF3D7_1200600_DBLpam2", 
                     "PF3D7_1200600_CIDRpam",
                     "PF3D7_1200600_DBLpam3",
                     "PF3D7_1200600_DBLepam4",
                     "PF3D7_1200600_DBLepam5", 
                     "PF3D7_1200600_DBLe10"))
PF3D7_1200600_DBLe10 <- VAR2CSA%>%
  filter(name=="PF3D7_1200600_DBLe10")

var2csa_plot <- ggboxplot(VAR2CSA, 
                          x="Gravida", 
                          y="value", 
                          add = "jitter", 
                          palette = "jco", 
                          color = "Gravida", 
                          xlab="Gravida", 
                          ylab = "Reactivity")+
  facet_wrap(~name, scales = "free", ncol = 4)+
  stat_compare_means(aes(group=Gravida), method = "t.test", label.y = 30)+
  scale_y_continuous(breaks = seq(0,30,5))+
  theme(axis.text = element_text(size = 10))

ggsave("var2csa_plot.tiff", height = 8, width = 12, dpi = 300)

###############################################################################
#THE END
