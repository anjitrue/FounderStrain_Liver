library(plyr)
library(tidyr)

library(ggplot2)
library(factoextra)

#Definne color palette
colors <- c("#206F94", # teal
            "#F47F72", #coral
            "#75C69D", #baby green
            "#5C66AF", #light purpleblue
            "#2CA8E0", #coon blue
            "#1F6F94", #darker blue
            "#ED237A", #hot pink)
            "#5C66AF", #light purpleblue
            "#2A4093", #dark purpleblue
            "#2C9379", #dark baby green
            "#83D5F7", #light coon blue
            "#93211E", #dark red
            "#E73C25", #red
            "#81143A" #dark pink
)

#Read csv
FS_Metabolomics<- read.csv("P:/EAT_20190812_DO_Search05/FounderStrain/QuantResults_202012111815107115_editKAO_forR.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(FS_Metabolomics) <- make.unique(FS_Metabolomics$Feature.ID)

FS_Metabolomics_df <- FS_Metabolomics[,-c(1)]

#Create dataframes for additional information retrieved in Nick's Software
FS_LfQ <- FS_Metabolomics_df[,grep("_LFQ", colnames(FS_Metabolomics_df))]
FS_tier <- FS_Metabolomics_df[,grep("_tier", colnames(FS_Metabolomics_df))]
FS_quantMZ <- FS_Metabolomics_df[,grep("_quantMZ", colnames(FS_Metabolomics_df))]
FS_RT <- FS_Metabolomics_df[,grep("_RT", colnames(FS_Metabolomics_df))]

#R
sample.names <- colnames(FS_LfQ)

#Create a meta document for the samples
Sample_meta <- as.data.frame(sample.names) %>% separate(sample.names, into = paste("V", 1:7, sep = "_"))

#Rename column names
colnames(Sample_meta) <- c("Date", "Scientist_Initials", "Sample_Type", "Strain", "Sample_ID", 
                           "Sample_Prep_Number", "Sample_Queue_Number")
#Use if else statement to check if the Date contains an X as the first character, then remove
  Sample_meta$Date <- ifelse(substr(Sample_meta$Date,1,1) == "X", sub("^.", "", Sample_meta$Date), Sample_meta$Date)
  
  Sample_meta$Strain <- sub("Cast", "CAST", Sample_meta$Strain)
#Concatonate the Sample_ID to contain strain and number information
  Sample_meta$Sample_ID <- paste0(Sample_meta$Strain,"-", Sample_meta$Sample_ID)
  
#Determine the number of mice in each strain
  table(Sample_meta$Strain)

#Write the meta table out into a csv and save in project folder
write.csv(Sample_meta, "P:/EAT_20190812_DO_Search05/FounderStrain/Meta_Founder_Strain_Metabolomics_Liver.csv")

#Transform metabolomics dataframe for PCA analysis
t_FS_Metabolomics_df <- t(FS_LfQ)


rownames(t_FS_Metabolomics_df) <- make.unique(Sample_meta$Sample_ID)


save(colors,FS_LfQ, FS_tier, FS_quantMZ, FS_RT, Sample_meta, FS_Metabolomics_df, t_FS_Metabolomics_df, 
     file = "P:/EAT_20190812_DO_Search05/FounderStrain/DataAnalysis/Founder_Strain_Metabolomics_data.Rdata")


