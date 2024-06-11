#########################################################
##Script made by: 
##Local: Universidade do Algarve
##Last Modification: 09/04/2024
##Project: Convert a .geno file into a usable structure file (.str) 
###################################################################################################################################

geno_to_str <- function(geno_input){

#### Load Package
packages = c("admixr","admixtools","dplyr","tidyverse","tidyr")
  
  ## Load or install all necessary packages
  package.check <- lapply(
    packages,
    FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        BiocManager::install(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    }
  )
  

#######
## Load the eigenstra
df <- read_eigenstrat(geno_input)
name_file <-as.character(deparse(substitute(geno_input)))
name_file <- gsub('"',"" ,test)
## Manually convert .geno file into a structure file
## Load Eigenstrat
#df <- read_eigenstrat("HOAH_7.5perc_MD")
#df <- geno_input
#######
## Covenvert into a easier workable table
df1 <- as.data.frame(df$geno)
df2 <- as.data.frame(t(df1))

#######
## Transpose table to match structure table 
geno <- as.data.frame(t(df1))

#######
## Repeat every row 2 times, to create 2 alleles
geno<-geno[rep(seq_len(nrow(geno)), each = 2), ]  

#######
## Replace all NAs by -9 (Means missing data)
geno[is.na(geno)] <- -9

#######
## Place the original object in a temporary object
geno.test<- geno

#######
## Identify each row
geno.test$IDs <- rownames(geno)

#######
## Identify the alleles, so we can separate into 2 tables, so its easier to prepare the structure (.str) file
geno.test$filter <- rep(c(1,2), times=length(df$ind$X1))

#######
## Grab the first allele and place it in one table
geno.Hap1<- filter(geno.test,grepl(paste(c("1"),collapse = '|'),filter))

#######
## Grab the second allele and place it in another table 
geno.Hap2<- filter(geno.test,!grepl(paste(c("1"),collapse = '|'),filter))

#######
## Transform all the "0s" in the first table into "1s" -> This means all Founder Homozigous are "0" in Geno and will be "1/1" in .str
geno.Hap1[geno.Hap1=="0"]<-1

#######
## Transform all the "1s" in the first table into "2s" -> This means all Heterozigous are "1" in Geno and will be "1/2" in .str
## So it means in the first table they will still be 1, but in the second table, all "1s" are transformed into "2s"
geno.Hap2[geno.Hap2=="1"]<-2

#######
## Transform all the "0s" in the second table into "1s" -> This means all Founder Homozigous are "0" in Geno and will be "1/1" in .str
geno.Hap2[geno.Hap2=="0"]<-1

#######
## PLEASE NOTE: No change are made in "2s" that are in both tables
## Because all Mutation Homozigous are "2" in geno and will be "2/2" in .str
#######

#######
## Concatenate both tables together to create an unfinished .str file
final_geno<-rbind(geno.Hap1,geno.Hap2)

#######
## Reordering the table to match the orginal order from the Eigenstrat .geno file
final_geno<-final_geno[order(match(final_geno[,"IDs"],geno.test[,"IDs"])),]

#######
## Clean the table to only have the SNPs present
final_geno <- final_geno[,1:length(df$snp$SNP)]

#######
## Create an annotation object that has a column with all the individuals repeated 2 times
final_anno <- as.data.frame(rep(colnames(df1),each=2))

#######
## Concatenate both the annotation object and the previously created unfinished .str file
final_file <- cbind(final_anno,final_geno)

#######
## Remove the header of the first collumn
colnames(final_file)[1] <-""

#######
## Save the finished .str file 
write.table(final_file ,paste0(name_file, "_tst.str",sep=""), row.names =F,sep="\t",quote=F)
#mapply (write.table, final_file, paste0(test, "_test.str", sep =""),row.names=F,sep="\t",quote=F)
}

## EOF
#######