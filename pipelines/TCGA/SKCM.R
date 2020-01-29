library('data.table')
library('stringr')
source("~/test/source/_lib.R")

setwd("~/Hutch-Research/Data/Processed")

# ------------------------------------------------------------
# read in pure cell type data
# ------------------------------------------------------------

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
setwd("~/Hutch-Research/Data/Real")

datM = fread(file = "skcm_methylation.txt", header=TRUE)
sampinfo = fread('skcm_sample_info.txt',header = FALSE)
cpg <- unlist(datM[,1])
colnames(datM) <- gsub('^X',"",colnames(datM))
datM <- datM[,-1]
rownames(datM) = cpg

tcga_purity <- fread('tcga_purity.txt')
barcode = list()
base_string = '%s-%s-%s'
for(i in 1:length(tcga_purity$array)){
  v1 <- strsplit(tcga_purity$array[i],'-')[[1]][1:3]
  temp <- do.call(sprintf,c(fmt = base_string, as.list(v1)))
  barcode = c(barcode,temp)
}

mat2 <- unlist(intersect(unlist(sampinfo[,2]),barcode))

datM <- subset(datM, select = unlist(sampinfo[which(unlist(sampinfo[,2]) %in% barcode ),3]))

purity_from_file <- tcga_purity[match(mat2,barcode),purity]
filena <- which(is.na(purity_from_file))
if(length(filena)>0){
    purity_from_file <- purity_from_file[-filena]
    datM <- subset(datM, select = -filena)
}
rownames(datM) <- cpg
purity <- purity_from_file
dim(datM)

# ------------------------------------------------------------
# read in probes to be used
# ------------------------------------------------------------

setwd("../../R_batch2")

load("ref_966probes.RData")
genes = probe2use
length(genes)
length(unique(genes))
genes = unique(genes)

table(genes %in% rownames(datM))
table(genes %in% rownames(dat))

X  = dat[match(genes, rownames(dat)),]
dim(X)
colnames(X)

dim(sam)
table(sam$id == colnames(X))
table(sam$label)
cellTypes = unique(sam$label)
# ------------------------------------------------------------
# extract methylation data from tumor samples
# ------------------------------------------------------------

ys = as.matrix(datM[match(rownames(X), rownames(datM)),])
rownames(ys) <- rownames(datM)[match(rownames(X), rownames(datM))]
dim(ys)
ys[1:2,1:5]

ys_na = which(apply(is.na(ys),2,any))
#eta_abs  = emInfo$abs_purity
eta_abs = purity
eta_abs_na = which(is.na(eta_abs))
nona = c(ys_na,eta_abs_na)
nona = as.vector(nona)
ys = ys[,-nona]
eta_abs = eta_abs[-nona]

#-------------------------------------------------------------
# Estimate Mean Matrix mu
#-------------------------------------------------------------
mu = matrix(NA, nrow = dim(X)[1],ncol = length(cellTypes))
s2 = matrix(NA, nrow = dim(X)[1],ncol = length(cellTypes))
row.names(mu) = row.names(s2) = rownames(X)
colnames(mu) = colnames(s2) =cellTypes

for(ct in cellTypes){
  dat.ct = X[,sam[which(sam[,2]==ct),1]]
  mu[,ct] = rowMeans(dat.ct,na.rm=TRUE)
  s2[,ct] = apply(dat.ct,1,sd, na.rm= TRUE)^2
}

#-------------------------------------------------------------
# Read Estimation from Expression Data
#-------------------------------------------------------------
setwd('~/Hutch-Research/Data/Real')
est_expr <- fread('SKCM_composition_cibersortx.txt')
samname <- as.list(est_expr[,1])
samname <- samname$Mixture
for(i in 1:length(samname)){
  samname[i] <- str_replace(samname[i],"X","")
}
length(samname)

com_sample <- intersect(samname,colnames(ys))

est_expr <- est_expr[which(samname %in% com_sample),-c(1,24,25,26)]
rownames(est_expr) <- samname[which(samname %in% com_sample)]

ys <- ys[,com_sample]
eta_abs <- eta_abs[which(colnames(ys) %in% com_sample)]

deconv_expr <- matrix(NA, nrow = nrow(est_expr),ncol = length(cellTypes))
colnames(deconv_expr) <- cellTypes 
rownames(deconv_expr) <- rownames(est_expr)
other <- rowSums(est_expr[,c(10,19,20,21)])
deconv_expr[,"B"] <- rowSums(est_expr[,1:3])/0.4
deconv_expr[,"CD4T"] <- rowSums(est_expr[,5:8])/0.4
deconv_expr[,"CD8T"] <- as.matrix(est_expr[,4])/0.4
deconv_expr[,"Treg"] <- as.matrix(est_expr[,9])/0.4
deconv_expr[,"NK"] <- rowSums(est_expr[,11:12])/0.42
deconv_expr[,"Monocyte"] <- rowSums(est_expr[,13:18])/1.4
deconv_expr[,"Neutrophil"] <- as.matrix(est_expr[,"Neutrophils"])/0.15
#est_expr$`T cells gamma delta` <- est_expr$`T cells gamma delta`
#est_expr$Eosinophils <- est_expr$Eosinophils
#est_expr$`Mast cells resting` <- est_expr$`Mast cells resting`
#est_expr$`Mast cells activated` <- est_expr$`Mast cells activated`
deconv_expr <- deconv_expr / rowSums(deconv_expr)
print(rownames(deconv_expr)[1])
cat('dim of deconv_expr', dim(deconv_expr))
