
source('../../source/_lib.R')
source("GetMethylationMatrix.R")


#----------------------------------------------------------
# Estimate distribution for each CpG and cell type
# Normal distribution for beta-value
#----------------------------------------------------------
use_beta = TRUE

if(! use_beta){
  X = log(X/(1-X))
}

mu = matrix(NA, nrow = dim(X)[1],ncol = length(cellTypes))
row.names(mu) = rownames(X)
colnames(mu) = cellTypes

mu_gen = matrix(NA, nrow = dim(X_gen)[1],ncol = length(cellTypes))
row.names(mu_gen) = rownames(X_gen)
colnames(mu_gen)  = cellTypes

s2 = matrix(NA,nrow=dim(X)[1], ncol = length(cellTypes))
row.names(s2) = rownames(X)
colnames(s2) = cellTypes

s2_gen = matrix(NA,nrow=dim(X_gen)[1], ncol = length(cellTypes))
row.names(s2_gen) = rownames(X_gen)
colnames(s2_gen) = cellTypes
id = 1:nrow(X_gen)

for(ct in cellTypes){
  dat.ct = X[,sam_est[which(sam_est[,2]==ct),1]]
  mu[,ct] = rowMeans(dat.ct,na.rm=TRUE)
  s2[,ct] = apply(dat.ct,1,sd,na.rm=TRUE)^2
  
  dat_gen_ct = X_gen[,sam_gen[which(sam_gen[,2]==ct),1]]
  mu_gen[,ct] = rowMeans(dat_gen_ct,na.rm=TRUE)
  s2_gen[,ct] = apply(dat_gen_ct,1,sd,na.rm=TRUE)^2
  
  mu_not_na  = which(!is.na(mu[,ct]))
  mugen_notna= which(!is.na(mu_gen[,ct]))
  s2_not_na  = which(!is.na(s2[,ct]))
  s2gen_notna= which(!is.na(s2_gen[,ct])) 
 
  id1 = intersect(id,mugen_notna)
  id  = intersect(id1,s2gen_notna)
}

mu     = mu[id,]
mu_gen = mu_gen[id,]
s2     = s2[id,]
s2_gen = s2_gen[id,]

table(is.na(s2))
table(is.na(s2_gen))
#----------------------------------------------------------
# Select Probes by looking at residuals
#----------------------------------------------------------
Rsq  = rep(0,ncol(mu))
resi = matrix(NA, nrow = nrow(mu), ncol = ncol(mu))
for(i in 1:length(cellTypes)){
  temp <- mu[,-i]
  temp <- temp - rowMeans(temp)
  Z   <- as.data.frame(cbind(mu[,i],temp))
  
  fit <- lm(Z[,1] ~ .-1 ,data = Z[,-1])
  Rsq[i] <- summary(fit)$adj.r.squared
  resi[,i] <- (fit$residuals)
}
print(Rsq)

if(1>0){
  for(i in 1:length(cellTypes)){
    pdf(sprintf('../figures/resi_plot_%s.pdf',cellTypes[i]), width = 13.5, height = 6)
    newplot <- hist(resi[,i],breaks = 100)
    print(newplot)
    dev.off()
  } 
}

remove_probe = list()
for(i in 1:length(cellTypes)){
  remove_probe[[cellTypes[i]]] = rownames(mu)[which(abs(resi[,i])<0.05)]
}
remove_probe = unlist(remove_probe)
remove_probe = unique(remove_probe)
remove_probe = intersect(rownames(mu), remove_probe)
length(remove_probe)
new_probe = setdiff(rownames(mu),remove_probe)
table(new_probe %in% rownames(mu))
table(new_probe %in% rownames(mu_gen))

mu = mu[new_probe,]
mu_gen = mu_gen[new_probe,]

s2 = s2[new_probe,]
s2_gen = s2_gen[new_probe,]

mu = as.matrix(mu)
mu_gen = as.matrix(mu_gen)
s2 = as.matrix(s2)
s2_gen = as.matrix(s2_gen)

print(dim(mu))
print(dim(mu_gen))

print(cor(mu))
save(s2,file = 's2.RData')
#-----------------------------------------------------------
# Check multicolinearity
#-----------------------------------------------------------

if(1<0){
pdf("../figures/corplot-mu.pdf")
#par(mar=c(2,2,2,2))
cormu = melt(cor(mu))
colnames(cormu)[3] = 'correlation'
het <- ggplot(data = cormu, aes(x=Var1, y=Var2, fill = correlation)) + 
      geom_tile() + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
print(het)
dev.off()
}
