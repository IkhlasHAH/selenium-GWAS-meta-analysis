# selenium-GWAS-meta-analysis
R script for selenium mesurments GWAS meta-analysis


#excel spreadsheets combination using R for the selenium data 

dtEvans <- excel_sheets("evans.xlsx") %>% map_df(~read_xlsx("Evans.xlsx"))

dtmarta <- excel_sheets("marta.xlsx") %>% map_df(~read_xlsx("marta.xlsx"))

dtYangPlasma <- excel_sheets("YangPlasma.xlsx") %>% map_df(~read_xlsx("YangPlasma.xlsx"))

dtYangSerum <- excel_sheets("YangSerum.xlsx") %>% map_df(~read_xlsx("YangSerum.xlsx"))


#############How to change columns' names############### 

colnames(dtevansF)
colnames(dtEvansF) <- c("SNP", "A1", "A2", "BETA", "SE", "p", "direction of effect", "CHR", "BP", " ", " ")


colnames(dtmarta)
colnames(dtmarta) <- c("SNP", "A1", "A2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "BETA", "SE", "p", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal", "CHR", "BP")

colnames(dtYangPlasma)
colnames(dtYangPlasma) <- c("CHR", "BP", "SNP", "A1", "A2", "EAF", "sample size", "BETA", "SE", "p")

colnames(dtYangSerum)
colnames(dtYangSerum) <- c("CHR", "BP", "SNP", "A1", "A2", "EAF", "sample size", "BETA", "SE", "p")


###################move these files to decuments as txt ####################### 

write.table (x=dtevansF, file='dtevansF.txt', sep=' ', row.names=FALSE, quote=FALSE)

write.table (x=dtmarta, file='dtMarta.txt', sep=' ', row.names=FALSE, quote=FALSE)

write.table (x=dtYangPlasma, file='dtYangPlasma.txt', sep=' ', row.names=FALSE, quote=FALSE)

write.table (x=dtYangSerum, file='dtYangSerum.txt', sep=' ', row.names=FALSE, quote=FALSE)


################download the file to documents########## 

write.table (x=bedformat, file='evans.bed', sep=' ', row.names=FALSE, quote=FALSE)





#########remove the header from the dataset table#################

#just delet it from your txt file 

##########remove the first column from a dataset#########################


# Assuming your matrix is called 'mat'
mat <- mat[, -1]  # Remove the first column

# Assuming your data frame is called 'df'
df <- df[-1, ]  # Remove the first row


########### to remove rows with NA results#########

bedformat <- na.omit(bedformat)



#######to retrieve SNPs posotions from NCBI for NA SNPs#####

BiocManager::install("BSgenome")
library(BSgenome)
BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
chromosomevec = as.character(seq(1, 22, by = 1))
all_snps_info = snpsBySeqname(snps, chromosomevec)

head(all_snps_info)

RSIDEvans <- makeGRangesFromDataFrame(new_df, seqnames.field = "CHR", start.field = "posGRCh37", end.field = "posGRCh37", keep.extra.columns = TRUE)




#############to upload data from your computer into R, as a txt file##############

gwasdata = read.table("C:/Users/hp/Documents/dtmarta.txt")

ucsc-SNPS-data = read.delim("C:/Users/hp/Documents/ucsc-SNPs.txt") 



#adding columns to the data 

dtEvans %>% mutate(id = paste0(CHR, ':', BP, '_', A1, '/', A2))

dtevansF <- dtEvans %>% mutate(id = paste0(CHR, ':', BP, '_', A1, '/', A2))


dtYangPlasma %>% mutate(id = paste0(CHR, ':', BP, '_', A1, '/', A2))

dtYangPlasmaF <- dtYangPlasma %>% mutate(id = paste0(CHR, ':', BP, '_', A1, '/', A2))

dtYangSerum %>% mutate(id = paste0(CHR, ':', BP, '_', A1, '/', A2))

dtYangSerumF <- dtYangPlasma %>% mutate(id = paste0(CHR, ':', BP, '_', A1, '/', A2))



##############adding columns to the data without creating new file, rather, use data$newcolumnname by doing so, you are going only to add a column 


dtmarta$BP <-str_split(dtmarta$SNP, ":")%>% sapply('[',1)

dtmarta$BP <-str_split(dtmarta$SNP, "_")%>% sapply('[',1)

dtmarta$bp <-str_split(dtmarta$BP, "_")%>% sapply('[',1)

dtmarta$bp <- substr(dtmarta$SNP, 
                         start = (regexpr(":", dtmarta$SNP) + 1),
                         stop = (regexpr("_", dtmarta$SNP) - 1))



#######to delete a column in R########### 

dtmarta <- subset(dtmarta, select = -BP)



###########take evans's NA SNPs from decuments as txt file and upload it on R###########


install.packages("openxlsx")
library(openxlsx)

setwd("C:/Users/hp/Documents/")

dataNA <- read.table(file = "C:/Users/hp/Documents/NA.txt", sep = "\t", header = TRUE)    

write.xlsx(data, "path/to/your/output/file.xlsx")      


###########create new dataframe from with a column from an existing ne################

new_df <- data.frame(rsID = dtevansF$SNP)


write.table (x=new_df, file='new_df.txt', sep=' ', row.names=FALSE, quote=FALSE)

############data filtering############## 


View(dtAragam)
> dtlevin_1 <- filter(dtLevin, effect_allele_frequency >= 0.05)

> QC_histogram(dataset, data_col = 1,
               +              save_name = "dataset", save_dir = getwd(),
               +              export_outliers = FALSE,
               +              filter_FRQ = NULL, filter_cal = NULL,
               +              filter_HWE = NULL, filter_imp = NULL,
               +              filter_NA = TRUE,
               +              filter_NA_FRQ = filter_NA, filter_NA_cal = filter_NA,
               +              filter_NA_HWE = filter_NA, filter_NA_imp = filter_NA,
               +              breaks = "Sturges",
               +              graph_name = colnames(dataset)[data_col],
               +              header_translations, check_impstatus = FALSE,
               +              ignore_impstatus = FALSE,
               +              T_strings = c("1", "TRUE", "yes", "YES", "y", "Y"),
               +              F_strings = c("0", "FALSE", "no", "NO", "n", "N"),
               +              NA_strings = c(NA, "NA", ".", "-"), ...)

#how to change column names 

colnames(dtlevin_1)
colnames(dtlevin_1) <- c("chr", "pb", "Efallele", "Othallele", "beta", "SE", "EAF", "pvale", "rsID")


#manhatin plot command 
manhattan_plot(dataset = dtlevin_1, chr = "chromosome", pvalue = "pvale", position = 'pb', fileName = "levinplot.png", plot.title = "manhatin Plot", plot.subtitle = "Levin manhatin Plot", p.threshold = 0.01, sig.threshold.log = -log10(5 * 10^-8), beta = NULL, std.error = NULL, check.columns = TRUE)


#how to start plink on command promt 
C:\Users\hp>cd C:\Users\hp\Desktop\plink_win64_20230116



> colnames(dtlevin_1) <- c("chr", "pb", "Efallele", "Othallele", "beta", "SE", "EAF", "pvale", "rsID")
> View(dtlevin_1)
> View(dtlevin_1)
> colnames(dtlevin_1) <- c("CHR", "BP", "A1", "A2", "BETA", "SE", "F", "P", "SNP")
> write.table(x=dtlevin_1, file = "dtlevin_1.txt", sep = " ")
> write.table(x=dt, file='dt.txt', sep=' ', row.names=FALSE, quote=FALSE)
> write.table(x=dtlevin_1, file='delevin_1.txt', sep=' ', row.names=FALSE, quote=FALSE)


> colnames(dtsonia_1) <- c("SNP", "CHR", "BP", "A1", "A2", "F", "BETA", "SE", "P")
> write.table(x= dtsonia_1, file = 'dtsonia_1.txt', sep = ' ', row.names = FALSE, quote = FALSE)
> write.table(x=dtlevin_1, file='dtlevin_1.txt', sep=' ', row.names=FALSE, quote=FALSE)


##########creat new columns###########

dtEvans %>% mutate(id=paste0(CHR,':',BP,'_',A1,'/',A2)

dtEvans %>% mutate(id = paste0(CHR, ':', BP, '_', A1, '/', A2))
                   
                   
dtEvans %>% mutate(id=CHR,':',BP,'_',A1,'/',A2)


########bayan's script#####
#Get the different unique combinations of signs (+, -, 0, ?) in the Direction column (e.g. --, ++, +-, -+, +0, -0, -? etc.). Here you can also see how many variants has which combination of signs.
summary(as.factor(meta$Direction))

#Will return the row numbers for those variants that have these values in their direction column. Please make sure that all options from the summary above is included in at least and only one of these lines:
variants_in_all <- which(meta$Direction %in% c("---", "--+", "--0", "-+-", "-++", "-+0", "+--", 
                                               "+-+", "+-0", "++-", "+++", "++0"))
variants_hunt_only <- which(meta$Direction %in% c("-??", "+??", "0??"))
variants_moba_only <- which(meta$Direction %in% c("?+?","?-?","?0?"))
variants_pivus_only <- which(meta$Direction %in% c("??+","??-","??0"))
variants_hunt_moba <- which(meta$Direction %in% c("-+?", "+-?", "++?", "--?"))
variants_hunt_pivus <- which(meta$Direction %in% c("-?-", "-?+", "+?-", "+?+"))
variants_moba_pivus <- which(meta$Direction %in% c("?-+", "?+-", "?++", "?--", "?+0"))

meta$N <- NA
meta$N[variants_in_all] <- n_hunt+n_moba+n_pivus
meta$N[variants_moba_only] <- n_moba
meta$N[variants_hunt_only] <- n_hunt
meta$N[variants_pivus_only] <- n_pivus
meta$N[variants_hunt_moba] <- n_hunt+n_moba
meta$N[variants_hunt_pivus] <- n_hunt+n_pivus
meta$N[variants_moba_pivus] <- n_moba+n_pivus

summary(as.factor(meta$N))
head (meta)




######change the column names of  na rows evans data#####

colnames(na_rows)
colnames(na_rows) <- c("name", "A1", "A2", "BETA", "SE", "p", "direction of effect", "chrom", "BP", "id", " ", " ")

colnames(dtevansF) <- c("SNP", "A1", "A2", "BETA", "SE", "p", "direction of effect", "CHR", "BP", "id", " ", " ")
colnames(ucscSNPs) <- c("bin", "chrom", "BP", "chromEnd", "name", "class")

#######delete the last two columns of a dataset#########

na_rows <- na_rows[, -((ncol(na_rows)-1):ncol(na_rows))]

dtevansF <- dtevansF[, -((ncol(dtevansF)-1):ncol(dtevansF))]


##############Seperate NA data to another data frame###########

na_rows <- FEvans[apply(FEvans, 1, function(row) any(is.na(row))), ]

marta-na_rows <- dtMarta[apply(dtMarta, 1, function(row) any(is.na(row))), ]

#########remove certain columns from dataset###########

ucscSNPs <- subset(ucscSNPs, select = -c(score, strand, refNCBI, refUCSC, observed, molType, valid, avHet, avHetSE, func, locType, weight))

ucscSNPs <- subset(ucscSNPs, select = -c(bin, chromEnd, class))


######to remove the duplicated rows #########

results1 <- result %>% unique()
results2 <- result %>% select(c(name, chrom, BP)) %>% unique()

tmp = ucscSNPs %>% select(c(name, chrom, BP)) %>% unique()

ucsc1 <- unique(ucscSNPs)

ucscSNPs <- unique(ucscSNPs$name)



##########################---new----###################################
#########################---start--####################################
################----for the thesis script-----#########################

#########uploading data files from your computer to R as csv##########

library(tidyverse)
library(readxl)
library(data.table)

FYangPlasma = fread("C:/Users/DELL/Documents/yangplasma.tsv")

FYangSerum = fread("C:/Users/DELL/Documents/yangserum.tsv")

FEvans = fread("C:/Users/DELL/Documents/dtevansF.txt")                                      

FMarta = fread("C:/Users/DELL/Documents/marta.tbl")

ucscSNPs = fread("C:/Users/DELL/Documents/thesis data/ucscSNPs.txt")  #this file is also uploaded no fill out the NA values in Evans data

lift_over_results = fread("C:/Users/DELL/Documents/lift_over_results.bed")

####after I lost the txt file of evans data, i wrote a script to convert 3 excel sheets to a txt file


######combine diffrenet exel spreadsheets into one file######

file_path <- "C:/Users/DELL/Documents/Evans.xlsx"

sheet1 <- read_excel("C:/Users/DELL/Documents/Evans.xlsx", sheet = 1)
sheet2 <- read_excel("C:/Users/DELL/Documents/Evans.xlsx", sheet = 2)
sheet3 <- read_excel("C:/Users/DELL/Documents/Evans.xlsx", sheet = 3)


sheet1$chromosome <- as.character(sheet1$chromosome)
sheet2$chromosome <- as.character(sheet2$chromosome)
sheet3$chromosome <- as.character(sheet3$chromosome)


convert_columns_to_character <- function(df) {
  df[] <- lapply(df, as.character)
  return(df)
}

sheet1 <- convert_columns_to_character(sheet1)
sheet2 <- convert_columns_to_character(sheet2)
sheet3 <- convert_columns_to_character(sheet3)


FEvans <- bind_rows(sheet1, sheet2, sheet3)


write.table(FEvans, "FEvans.txt", sep = "\t", row.names = FALSE)


###########to find NA values from ucscSNPs###### this is the script that i have wrote

##############Seperate NA data to another data frame###########

na_rows <- FEvans[apply(FEvans, 1, function(row) any(is.na(row))), ]

######designing a function##### 

fill_na_from_ucsc_SNPS <- function(na_rows, ucscSNPs) {
  merged_data <- merge(na_rows, ucscSNPs, by = "name", all.x = TRUE)
  filled_data <- merged_data %>%
    mutate(chrom = ifelse(is.na(chrom.x), chrom.y, chrom.x)) %>%
    mutate(BP = ifelse(is.na(BP.x), BP.y, BP.x)) %>%
    select(-chrom.x, -chrom.y, -BP.x, -BP.y)  # Remove redundant columns
  
  return(filled_data)
}
####applying the function####
result <- fill_na_from_ucsc_SNPS(na_rows, ucscSNPs) #the chromosome values were filled as (chrx) rather than 23



#######script from dr fouad aso for filling the navalues####

##############load libraries###########
library(data.table)
library(dplyr)
setwd("C:/Users/fouad/thesis")

## load input
evans <- fread('dtevansF.txt');
dbsnp <- fread('ucscSNPS.txt');

## fix evans colnames
colnames(FEvans) <-
  c('SNP','A1','A2','BETA','SE','p','doe','CHR','BP','id');


## I don't understand why every row in evans is repeated 3 times?!
dim(evans)
evans <- unique(evans);
dim(evans)

## split evans into NA at BP and nonNA
naEvans <- FEvans[is.na(FEvans$CHR),];
Evans <- FEvans[!is.na(FEvans$CHR),];

## preprocess dbsnp file (select required columns, remove indels, then

dbsnp <- ucscSNPs %>% select(c(name, chrom, chromStart, chromEnd, class)) %>% filter(class == 'single') %>% unique;
colnames(dbsnp) <- c('SNP','CHR', 'start', 'end', 'class')


#dbsnp <- dbsnp %>% select(-c(end)) %>% unique;
dbsnp <- dbsnp %>% filter(!grepl('_', CHR));
dbsnp$CHR <- gsub('chr','', dbsnp$CHR);
dbsnp$CHR <- gsub('X', '23', dbsnp$CHR);
dbsnp$CHR <- gsub('Y', '24', dbsnp$CHR);
dbsnp$CHR <- gsub('M', '25', dbsnp$CHR);
dbsnp$CHR <- as.integer(dbsnp$CHR);

## merge
df <- left_join(naEvans, dbsnp, by=c('SNP'))


###to combine df(filled na values) and evans rows ####

####remove CHR.x and BP columns
df <- df %>% select(-CHR.x, -BP, -class)

#### change column's names of df
colnames(df) <- c('SNP', 'A1', 'A2', 'BETA', 'SE', 'p', 'doe', 'CHR', 'BP')

###combine df and evans
convert_columns_to_character <- function(df) {
  df[] <- lapply(df, as.character)
  return(df)
}

df <- convert_columns_to_character(df)
Evans <- convert_columns_to_character(Evans)

combined_df <- bind_rows(df, Evans)



###check if there is still any NA 

na_rows <- combined_df[apply(combined_df, 1, function(row) any(is.na(row))), ] 

### there were still 19 snps with NA values i am going to remove them

#####remove the 19 rows that contains NA values from combined_df

combined_clean <- na.omit(combined_df)

write_delim(combined_clean, "C:/Users/DELL/Documents/evans_filled_NA.txt", delim = "\t")


################# bed format #################

combined_clean$CHR <- as.character(combined_clean$CHR)
combined_clean$BP <- as.numeric(combined_clean$BP)

bed_dfF <- combined_clean %>%
mutate(start = BP, end = BP + 1, name = paste("SNP_", row_number(), sep = "")) %>%
select(CHR, start, end, name)


write_delim(bed_dfF, "C:/Users/DELL/Documents/FEvans_bed.txt", delim = "\t")


#####take only the head of the bed data to see if t is working in the liftover
####check if it is working with the chromosome name 23 or x 

new_bed <- head(bed_df)
write_delim(bed_dfF, "C:/Users/DELL/Documents/bed_dfF.bed", separate())

write.table (x=new_bed, file='bed.txt', sep=, row.names=FALSE, quote=FALSE)

#####add chr before each chromosome number
bed_dfF$CHR <- paste("chr", bed_dfF$CHR, sep = "")
bed_dfF$CHR[bed_dfF$CHR == "chrX"] <- "chrx"
write.table (x=bed_dfF, file='bed.txt', sep=' ', row.names=FALSE, quote=FALSE)


bed_3columns <- subset(bed_dfF, select = -name)
write.table (x=bed_3columns, file='bed_3_columns.bed', sep=' ', row.names=FALSE, quote=FALSE, header=false)



#######to fix the bed file for the liftover

# Identify rows with non-integer values in the third column
problematic_rows <- which(!grepl("^\\d+$", bed_3columns$end))

# Print the problematic rows
print(bed_3columns[problematic_rows, ])

# Convert scientific notation to integer if needed
bed_3columns$end <- as.integer(bed_3columns$end)
bed_3columns$start <- as.integer(bed_3columns$start)
bed_3columns$CHR <- as.integer(bed_3columns$CHR)

write.table (x=bed_3columns, file='bed_3columns.txt', sep=' ', row.names=FALSE, quote=FALSE)


####بدي أجرب كمان مرة لليفت أوفر

bed_final <- bed_dfF[, c("CHR", "start", "end")] #remove the snp_name column

bed_final_clean <- na.omit(bed_final)  #remove NA values 

custom_format <- apply(bed_final_clean, 1, function(row) {
  paste(row[1], ":", row[2], "-", row[3], sep = "")
})     #creat a custom format for saving TXT

file_path <- "C:/Users/DELL/Documents/bed_final_clean.txt"
writeLines(custom_format, con = file_path)



########################

Se_ALSPAC_autosomes = fread("C:/Users/DELL/Documents/Se_ALSPAC_autosomes.txt")

Se_ALSPAC_ChrX = fread("C:/Users/DELL/Documents/Se_ALSPAC_ChrX.txt")

Se_QIMR_autosomes = fread("C:/Users/DELL/Documents/Se_QIMR_autosomes.txt")

Se_QIMR_ChrX = fread("C:/Users/DELL/Documents/Se_QIMR_ChrX.txt")



#######remove the not-converted snps from the original bed file
###I have two files, called bed_clean and file1
###i want to delete all of file1 rows from bed_clean on Rstudio

notconverted3Columns = fread("C:/Users/DELL/Desktop/file (1).txt")
bed_final_clean = fread("C:/Users/DELL/Documents/bed_final_clean.txt")


bed_final_clean$start <- as.numeric(bed_final_clean$start)
notconverted3Columns$start <- as.numeric(notconverted3Columns$start)

bed_final_clean$end <- as.numeric(bed_final_clean$end)
notconverted3Columns$end <- as.numeric(notconverted3Columns$end)

bed_without_converted <- anti_join(bed_final_clean, notconverted3Columns, by = names(bed_final_clean))

bed_without_converted$CHR <- as.character(bed_without_converted$CHR)
bed_without_converted$start <- as.numeric(bed_without_converted$start)
bed_without_converted$end <- as.numeric(bed_without_converted$end)

bed_without_converted <- apply(bed_without_converted, 2, trimws)

custom_format <- apply(bed_without_converted, 1, function(row) {
  paste(row[1], ":", row[2], "-", row[3], sep = "")
})     #creat a custom format for saving TXT

file_path <- "C:/Users/DELL/Documents/bed_out_n_converted.txt"
writeLines(custom_format, con = file_path)

###i have a data set that on my desk top each row looks like this: 
###chrx:115823324-115823325, 
###but i want to upload it on R as a three columns data frame that 
###the first contains chrx, 
###the second  this: 115823324 and 
###the third this: 115823325

notconverted <- read.table("C:/Users/DELL/Desktop/file (1).txt", header = FALSE, stringsAsFactors = FALSE)

notconverted3Columns <- notconverted %>%
  separate(V1, into = c("CHR", "position"), sep = ":") %>%
  separate(position, into = c("start", "end"), sep = "-")




##i have a txt file on desktop, the file data is arranged like this v1:v2-v3.
##and i want to upload it on R but i want it to be in three seperated columns
###and columns names are CHR, start, end 


# Set the file path (adjust the path according to your actual desktop location)
file_path <- "C:/Users/DELL/Downloads/liftover_results_out_notconverted.bed"

# Read the file into a dataframe
LO_results_out_nconverted <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Rename the column to 'V1' for processing
names(LO_results_out_nconverted) <- "V1"

# Use separate to split the data into three columns
LO_results_out_nconverted <- separate(LO_results_out_nconverted, col = "V1", into = c("CHR", "start", "end"), sep = "[:-]")

# View the resulting dataframe
print(data)


############
###i have two data sets, the first called LO the second called BED. 
###I  want to creat a data set that contains all of the 3 columns from LO 
###and one column that called start from BED in Rstudio 
###the common key to merge on is the arrangement 
###they are both have the same arrangement and the same rows number 

bed_without_converted <- as.data.frame(bed_without_converted)
LO_results_out_nconverted <- as.data.frame(LO_results_out_nconverted)

LO_results <- cbind(LO_results_out_nconverted, start_BED = bed_without_converted$start)

##################################################################

#have two data sets on R, the first one is called combined_clean and 
#the other is called LO_results both have BP, CHR and converted_BP columns. 
#However, the converted_BP values in combined_clean dataset are NA, 
#so what i want to do is to fill those NA values in  combined_clean$converted_BP 
#depending on the mutual values of BP and CHR columns in both datasets 

combined_clean$converted_BP <- NA #created new column called converted_BP in combined_clean dataset 

colnames(LO_results) <- c("CHR", "converted_BP", "end", "BP") #made the columns names identical in both datasets


LO_results$BP <- as.numeric(as.character(LO_results$BP)) #idintified BP cariables as numeric

LO_results$CHR <- gsub("chr", "", LO_results$CHR) #delet "chr" before each number in CHR column

combined_clean_filled <- combined_clean %>%
  left_join(LO_results, by = c("BP", "CHR")) %>%
  mutate(converted_BP = ifelse(is.na(converted_BP.x), converted_BP.y, converted_BP.x)) %>%
  select(-converted_BP.x, -converted_BP.y)

#delete all of the NA values 
combined_clean_filled <- combined_clean_filled[complete.cases(combined_clean_filled$converted_BP), ]
###now i have converted Evans data that is called "combined_clean_filled" 
###that contains onlt 2543142 SNPs

FEvans_converted <- combined_clean_filled

###change columns names to make them identical in all datasets

colnames(FEvans_converted)
colnames(FEvans_converted) <- c("rsID", "A1", "A2", "OR", "SE", "p", "doe", "CHR", "BP_36", "BP","SNP")


colnames(FMarta)
colnames(FMarta) <- c("SNPID", "A1", "A2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "OR", "SE", "p", "Doe", "HetISq", "HetChiSq", "HetDf", "HetPVal", "CHR", "BP")

colnames(FYangPlasma)
colnames(FYangPlasma) <- c("CHR", "BP", "rsID", "A1", "A2", "EAF", "sample size", "OR", "SE", "p", "id", "SNP")

colnames(FYangSerum)
colnames(FYangSerum) <- c("CHR", "BP", "rsID", "A1", "A2", "EAF", "sample size", "OR", "SE", "p", "id", "SNP")



#adding ID columns to the data such as marta's data

####Evans_converted <- subset(FEvans_converted, select = -SNP)
FEvans_converted <- FEvans_converted %>%
  mutate(SNP = paste0(CHR, ':', BP, '_', A1, '/', A2)) ####adding id column 
FEvans_converted$id <- toupper(FEvans_converted$id)###change the letters in "id" column to capital letters
FEvans_converted$A1 <- toupper(FEvans_converted$A1)
FEvans_converted$A2 <- toupper(FEvans_converted$A2)

write.table (x=FEvans_converted, file='FEvans_converted.txt', sep=' ', row.names=FALSE, quote=FALSE)

FYangPlasma <- FYangPlasma %>%
  mutate(SNP = paste0(CHR, ':', BP, '_', A1, '/', A2))
FYangPlasma$SNP <- toupper(FYangPlasma$id) ###change to capital letters

FYangSerum <- FYangSerum %>%
  mutate(SNP = paste0(CHR, ':', BP, '_', A1, '/', A2))
FYangSerum$SNP <- toupper(FYangSerum$id) ###change to capital letters


FMarta <- FMarta %>%
  mutate(SNP = paste0(CHR, ':', BP, '_', A1, '/', A2))
FMarta$A1 <- toupper(FMarta$A1)
FMarta$A2 <- toupper(FMarta$A2)

####creat CHR column in Marta's data 

FMarta <- FMarta %>%
  mutate(CHR = sapply(str_split(SNP, ":"), `[`, 1))

####Add BP column to Marta's dataset

FMarta <- FMarta %>%
  mutate(BP = sapply(str_split(SNP, "[:_]"), `[`, 2))

FMarta$A1 <- toupper(FMarta$A1)
FMarta$A2 <- toupper(FMarta$A2)

##download the updated vresions of data to run the meta analysis on plink 

write.table (x=FMarta, file='FinalMarta.txt', sep=' ', row.names=FALSE, quote=FALSE)

write.table (x=FEvans_converted, file='FinalEvans.txt', sep=' ', row.names=FALSE, quote=FALSE)

write.table (x=FYangPlasma, file='FinalYangPlasma.txt', sep=' ', row.names=FALSE, quote=FALSE)

write.table (x=FYangSerum, file='FinalYangSerum.txt', sep=' ', row.names=FALSE, quote=FALSE)


####Preparing the data for Metal meta analysis###
#######solve the direction of effect to find N (sample size)


setwd("C:/Users/DELL/Downloads/plink")
evans <- read.table("FinalEvans.txt", head =T)

evan_output <- "FinalEvans.txt"
n_qimr <- 2603
n_alspac <- 2874

evans <- fread(evan_output)
summary(as.factor(evans$doe))
variants_in_all <- which(evans$doe %in% c("--", "-+", "+-", "++", "0-", "0+", "0"))
variants_alspac_only <- which(evans$doe %in% c("-?", "+?", "0?"))
variants_qimr_only <- which(evans$doe %in% c("?+","?-","?0"))

evans$N <- NA
evans$N[variants_in_all] <- n_alspac+n_qimr
evans$N[variants_alspac_only] <- n_alspac
evans$N[variants_qimr_only] <- n_qimr


summary(as.factor(evans$N))

colnames(evans)
colnames(evans) <- c("rsID", "Allele1", "Allele2", "EFFECT", "STDERR", "PVALUE", "doe", "CHR", "BP_36", "BP","MARKER", "N")
write.table (x=evans, file='MEvans.txt', sep=' ', row.names=FALSE, quote=FALSE)


##Marta 

setwd("C:/Users/DELL/Downloads/plink")
marta <- read.table("FinalMarta.txt", head =T)

marta_output <- "FinalMarta.txt"
n_hunt <- 2819
n_moba <- 2812
marta <- fread(marta_output)


summary(as.factor(marta$Doe))

variants_in_all <- which(marta$Doe %in% c("--", "-+", "-0", "+-", "++", "+0", "00"))
variants_hunt_only <- which(marta$Doe %in% c("-?", "+?", "0?"))
variants_moba_only <- which(marta$Doe %in% c("?+","?-","?0"))

marta$N <- NA
marta$N[variants_in_all] <- n_hunt+n_moba
marta$N[variants_moba_only] <- n_moba
marta$N[variants_hunt_only] <- n_hunt

summary(as.factor(marta$N))

colnames(marta)
colnames(marta) <- c("SNPID", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "EFFECT", "STDERR", "PVALUE", "Doe", "HetISq", "HetChiSq", "HetDf", "HetPVal", "CHR", "BP", "MARKER", "N")



write.table (x=marta, file='MMarta.txt', sep=' ', row.names=FALSE, quote=FALSE)


yangP <- FYangPlasma
colnames(yangP)
colnames(yangP) <- c("CHR", "BP", "rsID", "Allele1", "Allele2", "EAF", "sample size", "EFFECT", "STDERR", "PVALUE", "MARKER")
write.table (x=yangP, file='MyangP.txt', sep=' ', row.names=FALSE, quote=FALSE)

yangS <- FYangSerum
colnames(yangS)
colnames(yangS) <- c("CHR", "BP", "rsID", "Allele1", "Allele2", "EAF", "sample size", "EFFECT", "STDERR", "PVALUE", "MARKER")
write.table (x=yangS, file='MyangS.txt', sep=' ', row.names=FALSE, quote=FALSE)

#######

FYangPlasma %>% filter(BP == 39814369)
FYangSerum %>% filter(BP == 39814369)
FMarta %>% filter(BP == 39814369)
FEvans_converted %>% filter(BP == 39814369)



#####I want to cheack how many base pairs are common among all datasets
##### but first i want to creat a new column that contains only CHR=BP in all datasets


FMarta <- FMarta %>%
  mutate(CHRBP = paste0(CHR, ':', BP))

FEvans_converted <- FEvans_converted %>%
  mutate(CHRBP = paste0(CHR, ':', BP))

FYangPlasma <- FYangPlasma %>%
  mutate(CHRBP = paste0(CHR, ':', BP))

FYangSerum <- FYangSerum %>%
  mutate(CHRBP = paste0(CHR, ':', BP))

# Perform joins

common_data <- FMarta %>%
  inner_join(FEvans_converted, by = "CHRBP") %>%
  inner_join(FYangPlasma, by = "CHRBP") %>%
  inner_join(FYangSerum, by = "CHRBP")

common_data_SNP <- FMarta %>%
  inner_join(FEvans_converted, by = "SNP") %>%
  inner_join(FYangPlasma, by = "SNP") %>%
  inner_join(FYangSerum, by = "SNP")

############
result <- left_join(FEvans_converted, FMarta, by = "CHRBP")
result <- left_join(result, FYangPlasma, by = "CHRBP")
result <- left_join(result, FYangSerum, by = "CHRBP")

result_SNP <- left_join(FEvans_converted, FMarta, by = "SNP")
result_SNP <- left_join(result_SNP, FYangPlasma, by = "SNP")
result_SNP <- left_join(result_SNP, FYangSerum, by = "SNP")

####result file#### 


Meta_data = fread("C:/Users/DELL/Downloads/plink/metaanalysis results/plink 600 000.meta")


Meta_data2 = fread("C:/Users/DELL/Downloads/plink/metaanalysis results/plink.meta")

meta_data3 = fread("C:/Users/DELL/Downloads/plink/plink.meta")

Meta_data %>% filter(P(R) == 0)

FMarta %>% filter(BP == 73510597)
FEvans_converted %>% filter(BP == 73510597)


#################
colnames (meta_data3) [1] = "Chromosome"
colnames (meta_data3) [2] = "Position"
colnames (meta_data3) [3] = "rsID"
colnames (meta_data3) [4] = "Allele1"
colnames (meta_data3) [5] = "Allele2"
colnames (meta_data3) [6] = "Weight"
colnames (meta_data3) [8] = "P-value"
colnames (meta_data3) [10] = "Beta"

write.table (x=meta_data3, file='metafinalresults.txt', sep=' ', row.names=FALSE, quote=FALSE)
