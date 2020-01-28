library(nnls)
library(data.table)

source('../../source/_lib.R')
setwd("~/Hutch-Research/Data/Processed")

#---------------------------------------------------------------
# Read cell type specific methylation data
#---------------------------------------------------------------


info = read.table("methylation_pure_ct_info.txt", as.is=TRUE, sep="\t",
                  header=TRUE, quote="")
dim(info)
info[1:2,]

dat = fread("methylation_pure_ct_rmPC2_data.txt")
dim(dat)
dat[1:2,1:5]

sam = read.table("methylation_pure_ct_sample.txt", as.is=TRUE,
                 sep="\t", header=TRUE)
dim(sam)
sam[1:2,]

table(names(dat) == sam$id)

dat = data.matrix(dat)
table(sam$label)
rownames(dat) = info$ID

# ------------------------------------------------------------
# read DNA methylation data
# ------------------------------------------------------------

#setwd("~/Hutch-Research/TCGA/COAD/")

datM = fread(file = "methylation_betaValue.txt", header=TRUE)
dim(datM)
datM[1:2, 1:5]

infoM = read.table(file = "methylation_info.txt", sep = "\t",
                   header = TRUE, as.is = TRUE, quote="")
dim(infoM)
infoM[1:2, 1:5]
table(infoM$CHR)

table(datM$id == infoM$Name)

datM = data.matrix(datM[,-1])
rownames(datM) = infoM$Name
dim(datM)
datM[1:2, 1:5]

# ------------------------------------------------------------
# take intersection
# ------------------------------------------------------------

table(info$ID == info$Name)
cpgs = intersect(info$ID, infoM$Name)
length(cpgs)

mat1  = match(cpgs, infoM$Name)
datM  = datM[mat1,]
infoM = infoM[mat1,]

dim(datM)
datM[1:2, 1:5]

dim(infoM)
infoM[1:2, 1:5]

mat2 = match(cpgs, info$ID)
dat  = dat[mat2,]
info = info[mat2,]

dim(dat)
dat[1:2, 1:5]

dim(info)
info[1:2, 1:5]

cellTypes = unique(sam[,2])
# ------------------------------------------------------------
# Splitting data into two parts
# ------------------------------------------------------------
dat_est = dat_gen = sam_est = sam_gen = NULL

uct = unique(sam$label)
uct

for(i in 1:length(uct)){
  dati = which(sam[,2]==cellTypes[i])
  halfi= dati[1:(length(dati)/2)]
  
  dat_est = cbind(dat_est,dat[,halfi])
  dat_gen = cbind(dat_gen,dat[,setdiff(dati,halfi)])
  
  sam_est = rbind(sam_est,sam[halfi,])
  sam_gen = rbind(sam_gen,sam[setdiff(dati,halfi),])
}

dim(dat_est)
dim(dat_gen)
table(sam_est[,2])
table(sam_gen[,2])


# ------------------------------------------------------------
# read DNA methylation sample information
# ------------------------------------------------------------

setwd("./COAD")

coad = read.table("patient_coad_short_table_nMut.txt", sep="\t",
                  as.is=TRUE, header=TRUE)
dim(coad)
coad[1:2,]

ff0    = "patient_coad_M_info_hyperMeth.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is = TRUE)
dim(emInfo)
emInfo[1, ]

table(emInfo$bcr_patient_barcode %in% coad$bcr_patient_barcode)
mat1 = match(emInfo$bcr_patient_barcode, coad$bcr_patient_barcode)
emInfo$nMut = coad$nMut[mat1]

table(colnames(datM) == emInfo$patient_id)


# ------------------------------------------------------------
# read in probes to be used
# ------------------------------------------------------------
load('~/Hutch-Research/Data/Processed/p2use_half_14.RData')
if(1<0){panT = c("CD4T", "CD8T", "Treg")
p2use = levels = list()

level1 = list()
level1[["lymphoid"]] = c("B", "CD4T", "CD8T", "Treg", "NK")
level1[["myeloid"]]  = c("Monocyte", "Neutrophil")

levels[[1]] = level1

level1 = list()
level1[["B_NK"]] = c("B", "NK")
level1[["T"]]    = panT

levels[[2]] = level1

level1 = list()
level1[["B"]]  = "B"
level1[["NK"]] = "NK"

levels[[3]] = level1

level1 = list()
level1[["CD4T_Tregs"]] = c("CD4T", "Tregs")
level1[["CD8T"]] = "CD8T"

levels[[4]] = level1

level1 = list()
level1[["CD4T"]] = "CD4T"
level1[["Tregs"]]   = "Tregs"

levels[[5]] = level1

level1 = list()

level1[["Neutrophil"]] = "Neutrophil"
level1[["Monocyte"]]   = "Monocyte"

levels[[6]] = level1



uct = unique(sam$label)
uct

for(k in 1:length(uct)){
  ct1 = uct[k]
  ct2 = paste("Non", ct1)
  
  level1 = list()
  
  level1[[ct1]] = ct1
  level1[[ct2]] = setdiff(uct, ct1)
  
  levels[[k+6]] = level1
}

level1 = list()

level1[["CD4T"]] = "CD4T"
level1[["CD8T"]]   = "CD8T"

levels[[length(uct)+7]] = level1

length(levels)

path  = "../../figures"
dataType = "methylation"

p2use = list()

for(kk in 1:length(levels)){
  cat(kk, date(), "\n")
  level1 = levels[[kk]]
  nms    = names(level1)
  nm1    = sprintf("%s_vs_%s", nms[1], nms[2])
  
  sam2use = which(sam_est$label %in% unlist(level1))
  nms     = names(level1)
  
  ctype1  = rep("A", nrow(sam_est))
  ctype1[which(sam_est$label %in% level1[[nms[1]]])] = "B"
  ctype1[which(sam_est$label %in% level1[[nms[2]]])] = "C"
  
  w2kpB = which(ctype1 == "B")
  w2kpC = which(ctype1 == "C")
  
  len = min(length(w2kpB),length(w2kpC),10)
  cnt = rep(0,dim(dat_est)[1])
  names(cnt) = rownames(dat_est)

  
  for(rep in 1:2){
    cat(rep,"\n")
    B_sample = sample(w2kpB,len)
    C_sample = sample(w2kpC,len)
    sam_sub = sam_est[c(B_sample,C_sample),]
    dat_sub = dat_est[,c(B_sample,C_sample)]
    temp    = probeSelect(level1, sam_sub, dat_sub, path=path, dataType=dataType)
    cnt[temp$id] = cnt[temp$id] + 1
  }
  p2use_nm1    = names(which(cnt>=1))
  p2use[[nm1]] = sample(p2use_nm1,min(100,length(p2use_nm1)))
}

nms = names(p2use)
for(nm1 in nms){
  p1 = p2use[[nm1]]
  p2use[[nm1]] = p1[which(p1 %in% cpgs)]
}

lapply(p2use, dim)
}

genes = NULL
for(lb1 in names(p2use)){
  genes = c(genes, p2use[[lb1]])
}

length(genes)
length(unique(genes))
genes = unique(genes)

table(genes %in% rownames(datM))
table(genes %in% rownames(dat))

X  = dat_est[match(genes, rownames(dat_est)),]
dim(X)
colnames(X)

X_gen  = dat_gen[match(genes, rownames(dat_gen)),]
dim(X_gen)
colnames(X_gen)

dim(sam_est)
table(sam_est$id == colnames(X))
table(sam_est$label)

setwd('~/Desktop/EMeth/pipelines/CombinedStudy')
# ------------------------------------------------------------
# extract methylation data from tumor samples
# ------------------------------------------------------------

#ys = datM[match(rownames(X), rownames(datM)),]
#dim(ys)
#ys[1:2,1:5]
#cellTypes = unique(sam[,2])
