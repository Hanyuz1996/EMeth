
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

#setwd("~/research/TCGA/COAD/")

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

setwd("../../../R_batch2")

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

ys = datM[match(rownames(X), rownames(datM)),]
dim(ys)
ys[1:2,1:5]

ys_na = which(apply(is.na(ys),2,any))
eta_abs  = emInfo$abs_purity
eta_abs_na = which(is.na(eta_abs))
nona = c(ys_na,eta_abs_na)
ys = ys[,-nona]
eta_abs = eta_abs[-nona]

#-------------------------------------------------------------
# Estimate Mean Matrix mu
#-------------------------------------------------------------
mu = matrix(NA, nrow = dim(X)[1],ncol = length(cellTypes))
s2 = matrix(NA, nrow = dim(X)[1],ncol = length(cellTypes))

row.names(mu) = row.names(s2) =rownames(X)
colnames(mu) = colnames(s2) = cellTypes

for(ct in cellTypes){
  dat.ct = X[,sam[which(sam[,2]==ct),1]]
  mu[,ct] = rowMeans(dat.ct,na.rm=TRUE)
  s2[,ct] = apply(dat.ct,1,sd,na.rm = TRUE)^2
}

#-------------------------------------------------------------
# Read Estimation from Expression Data
#-------------------------------------------------------------
setwd('~/Hutch-Research/Data/Real')
est_expr <- fread('COAD_composition_cibersortx.txt')
samname <- as.list(est_expr[,1])
samname <- samname$Mixture
for(i in 1:length(samname)){
  samname[i] <- str_replace(samname[i],"X","")
}
length(samname)

com_sample <- intersect(samname,colnames(ys))

est_expr <- est_expr[which(samname %in% com_sample),-1]
rownames(est_expr) <- samname[which(samname %in% com_sample)]

deconv_expr <- matrix(NA, nrow = nrow(est_expr),ncol = length(cellTypes))
colnames(deconv_expr) <- cellTypes 
rownames(deconv_expr) <- rownames(est_expr)
other <- rowSums(est_expr[,c(10,19,20,21)])
deconv_expr[,"B"] <- rowSums(est_expr[,1:3])/0.4
deconv_expr[,"CD4T"] <- rowSums(est_expr[,5:8])/0.4
deconv_expr[,"CD8T"] <- as.matrix(est_expr[,4])/0.4
deconv_expr[,"Treg"] <- as.matrix(est_expr[,9])/0.4
deconv_expr[,"NK"] <- rowSums(est_expr[,11:12])/0.42
deconv_expr[,"Monocyte"] <- rowSums(est_expr[,13:18])/1.40
deconv_expr[,"Neutrophil"] <- as.matrix(est_expr[,"Neutrophils"])/0.15
#est_expr$`T cells gamma delta` <- est_expr$`T cells gamma delta`
#est_expr$Eosinophils <- est_expr$Eosinophils
#est_expr$`Mast cells resting` <- est_expr$`Mast cells resting`
#est_expr$`Mast cells activated` <- est_expr$`Mast cells activated`
deconv_expr <- deconv_expr / rowSums(deconv_expr)

print(dim(deconv_expr))

#-----------------------------------------------------------
# Get Barcode for each sample and check purity data
#-----------------------------------------------------------
#tcga_purity <- fread('tcga_purity.txt')
#barcode = list()
#base_string = '%s-%s-%s'
#for(i in 1:length(tcga_purity$array)){
#  v1 <- strsplit(tcga_purity$array[i],'-')[[1]][1:3]
#  temp <- do.call(sprintf,c(fmt = base_string, as.list(v1)))
#  barcode = c(barcode,temp)
#}

#mat2 <- unlist(intersect(emInfo$bcr_patient_barcode,barcode))
#purity_from_file <- tcga_purity[match(mat2,barcode),purity]
#purity_abs <- emInfo[match(mat2,emInfo$bcr_patient_barcode), 'abs_purity']
#filena <- which(is.na(purity_from_file))
#absna <- which(is.na(purity_abs))
#nona <- c(filena,absna)
#purity_from_file <- purity_from_file[-nona]
#purity_abs <- purity_abs[-nona]
# correlation between purity_abs and purity_file = 0.92

#pdf('../../figures/AbsVSFile.pdf')
#tempdata = cbind(purity_from_file, purity_abs)
#colnames(tempdata) <- c("Pan_cancerPurity","AbsolutePurity")
#newplot <- ggplot(data = as.data.frame(tempdata), aes(x=Pan_cancerPurity,y=AbsolutePurity))+ xlim(0,1) + ylim(0,1) +
#  geom_point() + geom_abline(intercept = 0,slope = 1)
#print(newplot)
#dev.off()

#setwd("~/Hutch-Research/Data/Real")
#pc = readRDS('pan_cancer_somatic_mutation_summary.rds')
#purity = fread('./methylation_data/patient_coad_EMC_info_w_plate_nMut.txt')

#mat3 <- unlist(intersect(pc$barcode,barcode))
#purity_est_all <- pc[match(mat3, pc$barcode),'AscatPurity']
#purity_file_all <- as.matrix(tcga_purity[match(mat3,barcode),'purity'])
#estna <- which(is.na(purity_est_all))
#fileallna <- which(is.na(purity_file_all))
#allna <- c(estna,fileallna)
#purity_est_all <- purity_est_all[-allna]
#purity_file_all <- purity_file_all[-allna]
# correlation is 0.86

#pdf('../../figures/EstVSFile.pdf')
#tempdata = cbind(purity_file_all, purity_est_all)
#colnames(tempdata) <- c("Pan_cancerPurity","AscatPurity")
#newplot <- ggplot(data = as.data.frame(tempdata), aes(x=Pan_cancerPurity,y=AscatPurity))+ xlim(0,1) + ylim(0,1) +
 # geom_point() + geom_abline(intercept = 0,slope = 1)
#print(newplot)
#dev.off()










