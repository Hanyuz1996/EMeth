
library("GEOquery")

setwd("~/research/TCGA/_microE/methylation/data")

# ------------------------------------------------------------
# GSE35069: Reinius et al.
# ------------------------------------------------------------

gse3 = getGEO("GSE35069", destdir="~/research/TCGA/_microE/methylation/data")

show(gse3)
gse3 = gse3[[1]]

mDat3 = as.data.frame(gse3)
dim(mDat3)
mDat3[1:2,1:5]

pDat3 = pData(phenoData(gse3))
dim(pDat3)
names(pDat3)
pDat3[1:2,]

table(pDat3$source_name_ch1)

sDat3 = pDat3[,c("geo_accession", "source_name_ch1")]

names(sDat3) = c("id", "label")
sDat3[1:2,]
table(sDat3$label)

i450 = pData(featureData(gse3))
dim(i450)
i450[1:5,1:5]

# ------------------------------------------------------------
# select probes to be used
# ------------------------------------------------------------

names(i450)
i450[1:2,]

table(i450$Genome_Build)
table(i450$CHR)

table(i450$Methyl27_Loci)
table(i450$Probe_SNPs=="")
table(i450$Probe_SNPs_10=="")

w2use = (i450$Genome_Build=="37" & i450$CHR != "")
w2use = w2use & (i450$Probe_SNPs=="" & i450$Probe_SNPs_10=="")

table(w2use)

info = i450[which(w2use),]
dim(info)
info[1:2,]

info$ID = as.character(info$ID)
table(info$ID %in% names(mDat3))

setwd("~/research/Deconvolution/data")

write.table(info, file = "i450_filtered_probes.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# write out data
# ------------------------------------------------------------

mDat3 = mDat3[,match(info$ID,names(mDat3))]

dim(mDat3)
mDat3 = t(mDat3)
dim(mDat3)
mDat3[1:2,1:5]

mDat3 = data.frame(id=rownames(mDat3), mDat3)
dim(mDat3)
mDat3[1:2,1:5]

write.table(mDat3, file = "Reinius_methy_data.txt", append = FALSE,
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(sDat3, file = "Reinius_methy_samples.txt", append = FALSE,
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# ------------------------------------------------------------
# CD4T_CD8T_E-GEOD-71957
# ------------------------------------------------------------
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k")

sams = read.table("CD4T_CD8T_E-GEOD-71957/E-GEOD-71957.hyb.sdrf.txt", 
                  header=TRUE, sep="\t", as.is=TRUE, quote="")
dim(sams)
sams[1:2,]
names(sams)
names(sams)[8] = "diagonsis1"
names(sams)[9] = "diagonsis2"

length(unique(sams$Source.Name))

table(sams$Characteristics..organism.part., useNA="ifany")
table(sams$Characteristics..sex., useNA="ifany")
summary(sams$Characteristics..age., useNA="ifany")

table(sams$Characteristics..cell.type., useNA="ifany")
table(sams$Characteristics..cell.type., sams$FactorValue..cell.type)
table(sams$Material.Type, useNA="ifany")
table(sams$diagonsis1, sams$Material.Type, useNA="ifany")
table(sams$diagonsis2, sams$Material.Type, useNA="ifany")
table(sams$diagonsis1, sams$diagonsis2, useNA="ifany")

unique(sams$diagonsis1)

ww1 = which(sams$diagonsis1 == unique(sams$diagonsis1)[1])
sams$diagonsis1[ww1] = sams$diagonsis2[ww1]
table(sams$diagonsis1, useNA="ifany")

w2kp = which(sams$diagonsis1 == "Healthy")
length(w2kp)
sams = sams[w2kp,]
table(sams$Material.Type, useNA="ifany")
table(sams$Characteristics..cell.type., sams$Material.Type, useNA="ifany")

# ------------------------------------------------------------
# write out sample information
# ------------------------------------------------------------

names(sams)
cols2kp = c(1,3,4,5,7,8,13,16,22,29:30,40)
sams2kp = sams[,cols2kp]
dim(sams2kp)
sams2kp[1:4,]

## some samples are repated because the list both chanel (Cy3 and Cy5)
sams2kp = unique(sams2kp)
dim(sams2kp)

names(sams2kp) = gsub("Characteristics..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("Comment..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("\\.$", "", names(sams2kp), perl=TRUE)
names(sams2kp) = gsub(".", "_", names(sams2kp), fixed=TRUE)

dim(sams2kp)
sams2kp[1:4,]
apply(sams2kp[,c(5:9,11)], 2, table)

table(sams2kp$Material_Type, sams2kp$Array_Design_REF)

samsMethy = sams2kp[which(sams2kp$Material_Type=="genomic DNA"),]
dim(samsMethy)
samsMethy[1:4,]

table(gsub(" 1", "", samsMethy$Source_Name) == samsMethy$Assay_Name)

samsRNA = sams2kp[which(sams2kp$Material_Type=="total RNA"),]
dim(samsRNA)
samsRNA[1:4,]

setwd("~/research/Deconvolution/data")

write.table(samsMethy, file = "CD4T_CD8T_E-GEOD-71957_methy_sample.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

write.table(samsRNA, file = "CD4T_CD8T_E-GEOD-71957_RNA_sample.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# ------------------------------------------------------------
# read in methylation data
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k/CD4T_CD8T_E-GEOD-71957")
fls = list.files(pattern="_sample_table.txt", recursive=TRUE)
length(fls)

mDat = NULL

for(i in 1:nrow(samsMethy)){
  cat(i, date(), "\n")
  fnm   = samsMethy$Derived_Array_Data_File[i]
  fpath = fls[grep(fnm, fls)]
  
  if(length(fpath) != 1){
    stop("unexpected file path\n")
  }
  
  dati = read.table(fpath, header=TRUE, sep="\t", as.is=TRUE, na.string="null")
  dim(dati)
  dati[1:2,]
  
  summary(dati$Sample_31_CD8_Detection_Pval)
  table(dati$Sample_31_CD8_Detection_Pval > 0.01)
 
  if(i==1){
    probes = dati$Reporter.Identifier
  }else{
    if(length(probes) != nrow(dati)){
      stop("probe number are differnt\n")
    }else if(any(probes != dati$Reporter.Identifier)){
      stop("probes are differnt\n")
    }
  }
 
  mDat = cbind(mDat, dati$VALUE) 
}

dim(mDat)
rownames(mDat) = probes
colnames(mDat) = samsMethy$Assay_Name

table(probes %in% info$ID)
table(info$ID %in% probes)

w2kp = which(probes %in% info$ID)
mDat = mDat[w2kp,]

dim(mDat)
mDat[1:2,1:5]

mDat = data.frame(id=rownames(mDat), mDat)
dim(mDat)
mDat[1:2,1:5]

setwd("~/research/Deconvolution/data")

write.table(mDat, file = "CD4T_CD8T_E-GEOD-71957_methy_data.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# ------------------------------------------------------------
# ------------------------------------------------------------
# CD4T_E-GEOD-72364
# ------------------------------------------------------------
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k")

sams = read.table("CD4T_E-GEOD-72364/E-GEOD-72364.sdrf.txt", 
                  header=TRUE, sep="\t", as.is=TRUE, quote="")
dim(sams)
sams[1:2,]
names(sams)

length(unique(sams$Source.Name))

table(sams$Characteristics..sex., useNA="ifany")
summary(sams$Characteristics..age., useNA="ifany")

table(sams$Characteristics..cell.type., useNA="ifany")
table(sams$Material.Type, useNA="ifany")
table(sams$Characteristics..disease.state., useNA="ifany")

w2kp = which(sams$Characteristics..disease.state. == "normal")
length(w2kp)
sams = sams[w2kp,]
table(sams$Characteristics..cell.type., useNA="ifany")

# ------------------------------------------------------------
# write out sample information
# ------------------------------------------------------------

names(sams)
cols2kp = c(1, 3:6, 12, 18, 25:26, 36)
sams2kp = sams[,cols2kp]
dim(sams2kp)
sams2kp[1:4,]

## some samples are repated because they list both chanel (Cy3 and Cy5)
sams2kp = unique(sams2kp)
dim(sams2kp)

names(sams2kp) = gsub("Characteristics..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("Comment..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("\\.$", "", names(sams2kp), perl=TRUE)
names(sams2kp) = gsub(".", "_", names(sams2kp), fixed=TRUE)

dim(sams2kp)
sams2kp[1:4,]
apply(sams2kp[,c(2,4:7,9)], 2, table)

setwd("~/research/Deconvolution/data")

write.table(sams2kp, file = "CD4T_E-GEOD-72364_methy_sample.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# ------------------------------------------------------------
# read in methylation data
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k/CD4T_E-GEOD-72364")
fls = list.files(pattern="_sample_table.txt", recursive=TRUE)
length(fls)

mDat = NULL

for(i in 1:nrow(sams2kp)){
  cat(i, date(), "\n")
  fnm   = sams2kp$Derived_Array_Data_File[i]
  fpath = fls[grep(fnm, fls)]
  
  if(length(fpath) != 1){
    stop("unexpected file path\n")
  }
  
  dati = read.table(fpath, header=TRUE, sep="\t", as.is=TRUE, na.string="null")
  dim(dati)
  dati[1:2,]
  
  summary(dati$Detection.Pval)
  table(dati$Detection.Pval > 0.01)
  
  if(i==1){
    probes = dati$Reporter.Identifier
  }else{
    if(length(probes) != nrow(dati)){
      stop("probe number are differnt\n")
    }else if(any(probes != dati$Reporter.Identifier)){
      stop("probes are differnt\n")
    }
  }
  
  mDat = cbind(mDat, dati$VALUE) 
}

dim(mDat)
rownames(mDat) = probes
colnames(mDat) = sams2kp$Assay_Name

table(probes %in% info$ID)
table(info$ID %in% probes)

w2kp = which(probes %in% info$ID)
mDat = mDat[w2kp,]

dim(mDat)
mDat[1:2,1:5]

mDat = data.frame(id=rownames(mDat), mDat)
dim(mDat)
mDat[1:2,1:5]

setwd("~/research/Deconvolution/data")

write.table(mDat, file = "CD4T_E-GEOD-72364_methy_data.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# ------------------------------------------------------------
# ------------------------------------------------------------
# mono_CD4T_E-GEOD-56047
# ------------------------------------------------------------
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k")

sams = read.table("mono_CD4T_E-GEOD-56047/E-GEOD-56047.sdrf.txt", 
                  header=TRUE, sep="\t", as.is=TRUE, quote="")
dim(sams)
sams[1:2,]
names(sams)

length(unique(sams$Source.Name))

table(sams$Characteristics..racegendersite., useNA="ifany")
table(sams$Characteristics..well.id., useNA="ifany")

summary(sams$Characteristics..age., useNA="ifany")

apply(sams[,c(6,9:11,16)], 2, summary)
summary(rowSums(sams[,c(6,9:11,16)], na.rm=TRUE))

table(sams$Characteristics..cell.type., useNA="ifany")
table(sams$Material.Type, useNA="ifany")


table(sams$Material.Type, sams$Characteristics..well.id., useNA="ifany")
table(sams$Material.Type, sams$Characteristics..racegendersite., useNA="ifany")

# ------------------------------------------------------------
# write out sample information
# ------------------------------------------------------------

names(sams)
cols2kp = c(1, 3:6, 9:11, 16, 15, 17:18, 22, 29:30, 38)
sams2kp = sams[,cols2kp]
dim(sams2kp)
sams2kp[1:4,]

## some samples are repated because they list both chanel (Cy3 and Cy5)
sams2kp = unique(sams2kp)
dim(sams2kp)

names(sams2kp) = gsub("Characteristics..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("Comment..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("\\.$", "", names(sams2kp), perl=TRUE)
names(sams2kp) = gsub(".", "_", names(sams2kp), fixed=TRUE)

dim(sams2kp)
sams2kp[1:4,]
names(sams2kp)

apply(sams2kp[,c(2,10:13)], 2, table)

samsMethy = sams2kp[which(sams2kp$Material_Type=="genomic DNA"),]
dim(samsMethy)
samsMethy[1:4,]

samsMethy$altID = gsub("_peripheral_CD4 [methylation]", "", 
                       samsMethy$Sample_title, fixed=TRUE)

samsMethy$altID = gsub("_peripheral_CD14 [methylation]", "", 
                       samsMethy$altID, fixed=TRUE)

table(gsub(" 1", "", samsMethy$Source_Name) == samsMethy$Assay_Name)

samsRNA = sams2kp[which(sams2kp$Material_Type=="total RNA"),]
dim(samsRNA)
samsRNA[1:4,]


setwd("~/research/Deconvolution/data")

write.table(samsMethy, file = "mono_CD4T_E-GEOD-56047_methy_sample.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

write.table(samsRNA, file = "mono_CD4T_E-GEOD-56047_RNA_sample.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# ------------------------------------------------------------
# read in methylation data
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k/mono_CD4T_E-GEOD-56047")
fls = list.files(pattern="_sample_table.txt", recursive=TRUE)
length(fls)

## ok, this is not funny.. there is no array data for methylation
## hm... let's go to GEO then.
table(samsMethy$Derived_Array_Data_File)

# ------------------------------------------------------------
# read in methylation data from GEO? nope, it does not work
# ------------------------------------------------------------

# d1 = "~/research/Deconvolution/original_data/methy_450k/mono_CD4T_E-GEOD-56047"
# gse1 = getGEO("GSE56581", destdir=d1)

# ----------------------------------------------------------------
# mannually download data from GEO in two files,first CD4+ T cells
# ----------------------------------------------------------------

# first extract usealbe columns of the file
library(data.table)

setwd("~/research/Deconvolution/original_data/methy_450k/mono_CD4T_E-GEOD-56047")
CD4 = read.table("GSE56581_methylome_normalized.txt", sep="\t", 
            as.is=TRUE, header=TRUE, nrows=10)

dim(CD4)
CD4[1:5,1:10]
cols = c(1,grep("Mvalue", names(CD4)))
cols

date()
CD4 = fread("GSE56581_methylome_normalized.txt", sep="\t", 
            header=TRUE, select=cols)
date()

dim(CD4)
CD4[1:2,1:10]

table(CD4$ID_REF %in% info$ID)
table(info$ID %in% CD4$ID_REF)

w2kp = which(CD4$ID_REF %in% info$ID)
CD4 = CD4[w2kp,]

dim(CD4)
CD4[1:2,1:5]

samsMethyCD4 = samsMethy[which(samsMethy$Sample_source_name == "CD4+ cell"),]
dim(samsMethyCD4)
samsMethyCD4[1:2,]

colnames(CD4) = gsub(".Mvalue", "", colnames(CD4), fixed=TRUE)

dim(CD4)
CD4[1:2,1:5]

table(samsMethyCD4$altID %in% colnames(CD4))
table(colnames(CD4) %in% samsMethyCD4$altID)

cnms = samsMethyCD4$Assay_Name[match(samsMethyCD4$altID,colnames(CD4)[-1])]
colnames(CD4)[-1] = cnms

dim(CD4)
CD4[1:2,1:5]

for(k in 2:ncol(CD4)){
  CD4[[names(CD4)[k]]] = round(CD4[[names(CD4)[k]]], 7)
}

dim(CD4)
CD4[1:2,1:5]

setwd("~/research/Deconvolution/data")

write.table(CD4, file = "CD4T_E-GEOD-56047_methy_data.txt", append = FALSE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ----------------------------------------------------------------
# next CD14 monocytes
# ----------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k/mono_CD4T_E-GEOD-56047")
CD14 = read.table("GSE56046_methylome_normalized.txt", sep="\t", 
                 as.is=TRUE, header=TRUE, nrows=10)

dim(CD14)
CD14[1:5,1:10]
cols = c(1,grep("Mvalue", names(CD14)))
summary(cols)
cols[1:20]

date()
CD14 = fread("GSE56046_methylome_normalized.txt", sep="\t", 
            header=TRUE, select=cols)
date()

dim(CD14)
CD14[1:2,1:10]

table(CD14$ID_REF %in% info$ID)
table(info$ID %in% CD14$ID_REF)

w2kp = which(CD14$ID_REF %in% info$ID)
CD14 = CD14[w2kp,]

dim(CD14)
CD14[1:2,1:5]

samsMethyCD14 = samsMethy[which(samsMethy$Sample_source_name == "CD14+ cell"),]
dim(samsMethyCD14)
samsMethyCD14[1:2,]

colnames(CD14) = gsub(".Mvalue", "", colnames(CD14), fixed=TRUE)

dim(CD14)
CD14[1:2,1:5]

table(samsMethyCD14$altID %in% colnames(CD14))
table(colnames(CD14) %in% samsMethyCD14$altID)

cnms = samsMethyCD14$Assay_Name[match(samsMethyCD14$altID,colnames(CD14)[-1])]
colnames(CD14)[-1] = cnms

dim(CD14)
CD14[1:2,1:5]

for(k in 2:ncol(CD14)){
  CD14[[names(CD14)[k]]] = round(CD14[[names(CD14)[k]]], 7)
}

setwd("~/research/Deconvolution/data")

write.table(CD14, file = "CD14T_E-GEOD-56047_methy_data.txt", append = FALSE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# ------------------------------------------------------------
# mono_DC_MAC_E-GEOD-71837
# ------------------------------------------------------------
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k")

sams = read.table("mono_DC_MAC_E-GEOD-71837/E-GEOD-71837.sdrf.txt", 
                  header=TRUE, sep="\t", as.is=TRUE, quote="")
dim(sams)
sams[1:2,]
names(sams)

length(unique(sams$Source.Name))

table(sams$Characteristics..cell.type., useNA="ifany")
table(sams$Characteristics..donor.id., useNA="ifany")
table(sams$Characteristics..treated.with., useNA="ifany")
table(sams$Material.Type, useNA="ifany")

table(sams$Characteristics..donor.id., sams$Characteristics..cell.type.)
table(sams$Characteristics..donor.id., sams$Characteristics..treated.with.)
table(sams$Characteristics..cell.type., sams$Characteristics..treated.with.)

## well, this is a small data set, 
## we will only use monocyte, mature DCs, and mature macrophages (MACs)

cellTypes2use = "peripheral blood isolated CD14+ cells (monocytes)"
cellTypes2use = c(cellTypes2use, "activated/mature DCs", "activated/mature MACs")

sams = sams[which(sams$Characteristics..cell.type. %in% cellTypes2use),]
dim(sams)

names(sams)
cols2kp = c(1, 3:6, 10, 14, 21, 30)
sams2kp = sams[,cols2kp]
dim(sams2kp)
sams2kp[1:4,]

## some samples are repated because they list both chanel (Cy3 and Cy5)
sams2kp = unique(sams2kp)
dim(sams2kp)

names(sams2kp) = gsub("Characteristics..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("Comment..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("\\.$", "", names(sams2kp), perl=TRUE)
names(sams2kp) = gsub(".", "_", names(sams2kp), fixed=TRUE)

dim(sams2kp)
sams2kp[1:4,]
apply(sams2kp[,c(4:6)], 2, table)

setwd("~/research/Deconvolution/data")

write.table(sams2kp, file = "mono_DC_MAC_E-GEOD-71837_methy_sample.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# ------------------------------------------------------------
# read in methylation data
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k/mono_DC_MAC_E-GEOD-71837")
fls = list.files(pattern="_sample_table.txt", recursive=TRUE)
length(fls)

mDat = NULL

for(i in 1:nrow(sams2kp)){
  cat(i, date(), "\n")
  fnm   = sams2kp$Derived_Array_Data_File[i]
  fpath = fls[grep(fnm, fls)]
  
  if(length(fpath) != 1){
    stop("unexpected file path\n")
  }
  
  dati = read.table(fpath, header=TRUE, sep="\t", as.is=TRUE, na.string="null")
  dim(dati)
  dati[1:2,]
  
  summary(dati$Detection.Pval)
  table(dati$Detection.Pval > 0.01)
  
  if(i==1){
    probes = dati$Reporter.Identifier
  }else{
    if(length(probes) != nrow(dati)){
      stop("probe number are differnt\n")
    }else if(any(probes != dati$Reporter.Identifier)){
      stop("probes are differnt\n")
    }
  }
  
  mDat = cbind(mDat, dati$VALUE) 
}

dim(mDat)
rownames(mDat) = probes
colnames(mDat) = sams2kp$Assay_Name

table(probes %in% info$ID)
table(info$ID %in% probes)

w2kp = which(probes %in% info$ID)
mDat = mDat[w2kp,]

dim(mDat)
mDat[1:2,1:5]

mDat = data.frame(id=rownames(mDat), mDat)
dim(mDat)
mDat[1:2,1:5]

setwd("~/research/Deconvolution/data")

write.table(mDat, file = "mono_DC_MAC_E-GEOD-71837_methy_data.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# ------------------------------------------------------------
# ------------------------------------------------------------
# neut_E-GEOD-65097
# ------------------------------------------------------------
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k")

sams = read.table("neut_E-GEOD-65097/E-GEOD-65097.sdrf.txt", 
                  header=TRUE, sep="\t", as.is=TRUE, quote="")
dim(sams)
sams[1:2,]
names(sams)

length(unique(sams$Source.Name))

table(sams$Comment..Sample_source_name., useNA="ifany")
table(sams$Characteristics..cell.type., useNA="ifany")
table(sams$Comment..Sample_source_name., sams$Characteristics..cell.type.)

table(sams$Characteristics..disease.state., useNA="ifany")
table(sams$Characteristics..disease.state., sams$Characteristics..cell.type.)

table(sams$Material.Type, useNA="ifany")

## well, this is a small data set, 
## we will only use monocyte, mature DCs, and mature macrophages (MACs)

cellTypes2use = "normal neutrophils"
sams = sams[which(sams$Characteristics..cell.type. %in% cellTypes2use),]
dim(sams)

names(sams)
cols2kp = c(1, 3:5, 8:10, 23, 30, 39)
sams2kp = sams[,cols2kp]
dim(sams2kp)
sams2kp[1:4,]

## some samples are repated because they list both chanel (Cy3 and Cy5)
sams2kp = unique(sams2kp)
dim(sams2kp)

names(sams2kp) = gsub("Characteristics..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("Comment..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("\\.$", "", names(sams2kp), perl=TRUE)
names(sams2kp) = gsub(".", "_", names(sams2kp), fixed=TRUE)

dim(sams2kp)
sams2kp[1:4,]
apply(sams2kp[,c(2,4:8)], 2, table)

setwd("~/research/Deconvolution/data")

write.table(sams2kp, file = "neut_E-GEOD-65097_methy_sample.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# ------------------------------------------------------------
# read in methylation data
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k/neut_E-GEOD-65097")
fls = list.files(pattern="_sample_table.txt", recursive=TRUE)
length(fls)

mDat = NULL

for(i in 1:nrow(sams2kp)){
  cat(i, date(), "\n")
  fnm   = sams2kp$Derived_Array_Data_File[i]
  fpath = fls[grep(fnm, fls)]
  
  if(length(fpath) != 1){
    stop("unexpected file path\n")
  }
  
  dati = read.table(fpath, header=TRUE, sep="\t", as.is=TRUE, na.string="null")
  dim(dati)
  dati[1:2,]
  
  summary(dati$Detection.Pval)
  table(dati$Detection.Pval > 0.01)
  
  if(i==1){
    probes = dati$Reporter.Identifier
  }else{
    if(length(probes) != nrow(dati)){
      stop("probe number are differnt\n")
    }else if(any(probes != dati$Reporter.Identifier)){
      stop("probes are differnt\n")
    }
  }
  
  mDat = cbind(mDat, dati$VALUE) 
}

dim(mDat)
rownames(mDat) = probes
colnames(mDat) = sams2kp$Assay_Name

table(probes %in% info$ID)
table(info$ID %in% probes)

w2kp = which(probes %in% info$ID)
mDat = mDat[w2kp,]

dim(mDat)
mDat[1:2,1:5]

mDat = data.frame(id=rownames(mDat), mDat)
dim(mDat)
mDat[1:2,1:5]

setwd("~/research/Deconvolution/data")

write.table(mDat, file = "neut_E-GEOD-65097_methy_data.txt", append = FALSE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# ------------------------------------------------------------
# NK_E-GEOD-66562
# ------------------------------------------------------------
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k")

sams = read.table("NK_E-GEOD-66562/E-GEOD-66562.sdrf.txt", 
                  header=TRUE, sep="\t", as.is=TRUE, quote="")
dim(sams)
sams[1:2,]
names(sams)

length(unique(sams$Source.Name))

table(sams$Comment..Sample_source_name., useNA="ifany")
table(sams$Characteristics..cell.type., useNA="ifany")
table(sams$Comment..Sample_source_name. == sams$Characteristics..cell.type.)

table(sams$Characteristics..donor.status., useNA="ifany")
table(sams$Material.Type, useNA="ifany")

## well, this is not ideal since all are CMV donor
## we will all types of NK cells

sams = sams[grep("NK cells", sams$Characteristics..cell.type),]
dim(sams)

names(sams)
cols2kp = c(1, 3:7, 11, 17, 24:25, 35)
sams2kp = sams[,cols2kp]
dim(sams2kp)
sams2kp[1:4,]

## some samples are repated because they list both chanel (Cy3 and Cy5)
sams2kp = unique(sams2kp)
dim(sams2kp)

names(sams2kp) = gsub("Characteristics..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("Comment..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("\\.$", "", names(sams2kp), perl=TRUE)
names(sams2kp) = gsub(".", "_", names(sams2kp), fixed=TRUE)

dim(sams2kp)
sams2kp[1:4,]
apply(sams2kp[,c(4:8,10)], 2, table)

table(sams2kp$cell_type, sams2kp$donor_id)

setwd("~/research/Deconvolution/data")

write.table(sams2kp, file = "NK_E-GEOD-66562_methy_sample.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# ------------------------------------------------------------
# read in methylation data
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k/NK_E-GEOD-66562")
fls = list.files(pattern="_sample_table.txt", recursive=TRUE)
length(fls)

mDat = NULL

for(i in 1:nrow(sams2kp)){
  cat(i, date(), "\n")
  fnm   = sams2kp$Derived_Array_Data_File[i]
  fpath = fls[grep(fnm, fls)]
  
  if(length(fpath) != 1){
    stop("unexpected file path\n")
  }
  
  dati = read.table(fpath, header=TRUE, sep="\t", as.is=TRUE, na.string="null")
  dim(dati)
  dati[1:2,]
  
  summary(dati$Detection.Pval)
  table(dati$Detection.Pval > 0.01)
  
  if(i==1){
    probes = dati$Reporter.Identifier
  }else{
    if(length(probes) != nrow(dati)){
      stop("probe number are differnt\n")
    }else if(any(probes != dati$Reporter.Identifier)){
      stop("probes are differnt\n")
    }
  }
  
  mDat = cbind(mDat, dati$VALUE) 
}

dim(mDat)
rownames(mDat) = probes
colnames(mDat) = sams2kp$Assay_Name

table(probes %in% info$ID)
table(info$ID %in% probes)

w2kp = which(probes %in% info$ID)
mDat = mDat[w2kp,]

dim(mDat)
mDat[1:2,1:5]

mDat = data.frame(id=rownames(mDat), mDat)
dim(mDat)
mDat[1:2,1:5]

setwd("~/research/Deconvolution/data")

write.table(mDat, file = "NK_E-GEOD-66562_methy_data.txt", append = FALSE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# ------------------------------------------------------------
# NK_E-GEOD-66562
# ------------------------------------------------------------
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k")

sams = read.table("Treg_CD4T_E-GEOD-49667/E-GEOD-49667.sdrf.txt", 
                  header=TRUE, sep="\t", as.is=TRUE, quote="")
dim(sams)
sams[1:2,]
names(sams)

length(unique(sams$Source.Name))

table(sams$Comment..Sample_source_name., useNA="ifany")
table(sams$Characteristics..cell.type., useNA="ifany")
table(sams$Comment..Sample_source_name. == sams$Characteristics..cell.type.)

table(sams$Characteristics..disease.state., useNA="ifany")
table(sams$Material.Type, useNA="ifany")
table(sams$Characteristics..sex., useNA="ifany")
table(sams$Characteristics..donor.id., useNA="ifany")

table(sams$Characteristics..donor.id., sams$Characteristics..cell.type.)

sams = sams[grep("rTreg", sams$Characteristics..cell.type),]
dim(sams)

names(sams)
cols2kp = c(1, 3:6, 9, 13, 23, 30,39)
sams2kp = sams[,cols2kp]
dim(sams2kp)
sams2kp

## some samples are repated because they list both chanel (Cy3 and Cy5)
sams2kp = unique(sams2kp)
dim(sams2kp)

names(sams2kp) = gsub("Characteristics..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("Comment..", "", names(sams2kp), fixed=TRUE)
names(sams2kp) = gsub("\\.$", "", names(sams2kp), perl=TRUE)
names(sams2kp) = gsub(".", "_", names(sams2kp), fixed=TRUE)

dim(sams2kp)
sams2kp

table(sams2kp$cell_type, sams2kp$donor_id)

setwd("~/research/Deconvolution/data")

write.table(sams2kp, file = "Treg_E-GEOD-49667_methy_sample.txt", 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# ------------------------------------------------------------
# read in methylation data
# ------------------------------------------------------------

setwd("~/research/Deconvolution/original_data/methy_450k/Treg_CD4T_E-GEOD-49667")
fls = list.files(pattern="_sample_table.txt", recursive=TRUE)
length(fls)

mDat = NULL

for(i in 1:nrow(sams2kp)){
  cat(i, date(), "\n")
  fnm   = sams2kp$Derived_Array_Data_File[i]
  fpath = fls[grep(fnm, fls)]
  
  if(length(fpath) != 1){
    stop("unexpected file path\n")
  }
  
  dati = read.table(fpath, header=TRUE, sep="\t", as.is=TRUE, na.string="null")
  dim(dati)
  dati[1:2,]
  
  summary(dati$Detection.Pval)
  table(dati$Detection.Pval > 0.01)
  
  if(i==1){
    probes = dati$Reporter.Identifier
  }else{
    if(length(probes) != nrow(dati)){
      stop("probe number are differnt\n")
    }else if(any(probes != dati$Reporter.Identifier)){
      stop("probes are differnt\n")
    }
  }
  
  mDat = cbind(mDat, dati$VALUE) 
}

dim(mDat)
rownames(mDat) = probes
colnames(mDat) = sams2kp$Assay_Name

table(probes %in% info$ID)
table(info$ID %in% probes)

w2kp = which(probes %in% info$ID)
mDat = mDat[w2kp,]

dim(mDat)
mDat[1:2,]

mDat = data.frame(id=rownames(mDat), mDat)
dim(mDat)
mDat[1:2,1:5]

setwd("~/research/Deconvolution/data")

write.table(mDat, file = "Treg_E-GEOD-49667_methy_data.txt", append = FALSE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


quit(save="no")


