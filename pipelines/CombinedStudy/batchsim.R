set.seed(1750905)
library('quantreg')
library('gridExtra')
library('cowplot')

dir  = getwd()
source("Simulation.R")
  

simsize  = 10 
simnoise = 2   #(1,2,3,4,5,10)
simpi = 0.1       #(0.01,0.1,0.5,0.9)
alpha = rep(1/7,7)
#penalty = 1000
penalty= (dim(mu)[1])*(10^seq(-2,1,0.5)) 
cellnum = 250
#true.sigma = c(simnoise*0.05,0.05)
reptime = 2      #1:100
methods = c("LaplaceEM","OriEM","svr","ls","rls","qp")
res.err = res.cor = res.rss = array(data = NA, dim = c(length(cellTypes),length(methods),reptime),
                          dimnames = list(cellTypes,methods,1:reptime))
#res.sigma = matrix(NA, nrow = reptime, ncol = 2)
res.pi    = matrix(-1,reptime,4)
res.rho   = array(data = NA, dim = c(simsize,length(cellTypes),length(methods)+1),
                    dimnames = list(1:simsize,cellTypes,c(methods,"true")))
res.nu0   = matrix(-2,reptime,4)
res.iter  = matrix(0,reptime,4)
res.confu = array(data = NA, dim = c(reptime, 4,4),
                  dimnames = list(1:reptime,c("LaplaceEM","MaxVarEM","BinomEM","OriEM"),c("tt","tf","ft","ff")))

res.sigma_c = matrix(-1,reptime,4)
res.sigma_a = matrix(-1,reptime,4)


for(i in 1:reptime){{
  cat("-------------------\n")
  cat(i, date(), "\n") 
result        = runsim(simsize,simnoise,simpi,aber = TRUE,reptime=i,penalty = penalty,cellnum)
true          = result$true
res.rho       = result$rho
res.err[,,i]  = result$err
res.cor[,,i]  = result$cor
res.rss[,,i]  = result$rss
temp   = true$idic
temp_t = which(temp == TRUE)
temp_f = which(temp == FALSE)
total  = length(temp)
 for(j in 1:2){{
     res.nu0[i,j] = cor(result[[j+3]]$nu0, true$nu0)
     res.pi[i,j]  = mean(result[[j+3]]$pi_a)
     res.sigma_c[i,j]  = result[[j+3]]$sigma_c
     res.sigma_a[i,j]  = result[[j+3]]$sigma_c*result[[j+3]]$lambda
     res.iter[i,j]= result[[j+3]]$iter
  
     res.idic = (result[[j+3]]$gamma)>0.5
     res.confu[i,1,j] = mean(res.idic[temp_t] == TRUE)
     res.confu[i,2,j] = mean(res.idic[temp_t] == FALSE)
     res.confu[i,3,j] = mean(res.idic[temp_f] == TRUE)
     res.confu[i,4,j] = mean(res.idic[temp_f] == FALSE)
  }}
  cat(i,'-th repetition ends on', date())
}}

setting = sprintf('%d_%d_%d',cellnum,simsize,simnoise)
save(res.rho,file = sprintf('rho_%s.RData', setting) )
save(res.err,file = sprintf('err_%s.RData', setting))
save(res.cor,file = sprintf('cor_%s.RData', setting))
save(res.nu0,file = sprintf('nu0_%s.RData', setting))

dir.create('../figures')
dir.create(sprintf('../figures/%s',Sys.Date()))
  
pdf(sprintf("../figures/%s/nu0_%d_%d_%d.pdf",Sys.Date(),cellnum,simsize,simnoise), width = 20, height = 8)
nu0_plot = data.frame(cbind(result[[7]]$nu0,true$nu0))
nu0_plot_het = data.frame(cbind(result[[6]]$nu0,true$nu0))
nu0_plot_weight = data.frame(cbind(result[[5]]$nu0,true$nu0))
nu0_plot_laplace = data.frame(cbind(result[[4]]$nu0,true$nu0))

colnames(nu0_plot) = c('nu0_est','nu0_true')
colnames(nu0_plot_het) = c('nu0_est_het','nu0_true')
colnames(nu0_plot_weight) = c('nu0_est_weight','nu0_true')
colnames(nu0_plot_laplace) = c('nu0_est_laplace','nu0_true')

nu0scatter <- ggplot(data = nu0_plot, aes(x=nu0_true, y=nu0_est)) + xlim(0,1) + ylim(0,1) +
  geom_point() + geom_abline(intercept = 0,slope = 1)

nu0scatter_het <- ggplot(data = nu0_plot_het, aes(x=nu0_true, y=nu0_est_het)) + xlim(0,1) + ylim(0,1) +
  geom_point() + geom_abline(intercept = 0,slope = 1)
nu0weight <- ggplot(data = nu0_plot_weight, aes(x=nu0_true, y=nu0_est_weight)) + xlim(0,1) + ylim(0,1) +
  geom_point() + geom_abline(intercept = 0,slope = 1)
nu0laplace <- ggplot(data = nu0_plot_laplace, aes(x=nu0_true, y=nu0_est_laplace)) + xlim(0,1) + ylim(0,1) +
  geom_point() + geom_abline(intercept = 0,slope = 1)
grid.arrange(nu0scatter,nu0scatter_het,nu0weight,nu0laplace,ncol =2)

dev.off()


pdf(sprintf("../figures/%s/rho_%d_%d_%d.pdf",Sys.Date(),cellnum,simsize,simnoise), width = 20, height = 8)
res_rho_plot = data.frame(res.rho[1,,])

EMWeight <-ggplot(data = res_rho_plot, aes(x=true, y=MaxVarEM)) + xlim(0,0.5) + ylim(0,0.5) +
                geom_point() + geom_abline(intercept = 0,slope = 1)

EMLaplace <-ggplot(data = res_rho_plot, aes(x=true, y=LaplaceEM)) + xlim(0,0.5) + ylim(0,0.5) +
                geom_point() + geom_abline(intercept = 0,slope = 1)

EMHet <-ggplot(data = res_rho_plot, aes(x=true, y=BinomEM)) + xlim(0,0.5) + ylim(0,0.5) +
                geom_point() + geom_abline(intercept = 0,slope = 1)
EMOri <-ggplot(data = res_rho_plot, aes(x=true, y=OriEM)) + xlim(0,0.5) + ylim(0,0.5) +
                geom_point() + geom_abline(intercept = 0,slope = 1)
svr <-ggplot(data = res_rho_plot, aes(x=true, y=svr)) + xlim(0,0.5) + ylim(0,0.5) +
                geom_point() + geom_abline(intercept = 0,slope = 1)
ls <-ggplot(data = res_rho_plot, aes(x=true, y=ls)) + xlim(0,0.5) + ylim(0,0.5) +
                geom_point() + geom_abline(intercept = 0,slope = 1)
rls <-ggplot(data = res_rho_plot, aes(x=true, y=rls)) + xlim(0,0.5) + ylim(0,0.5) +
                geom_point() + geom_abline(intercept = 0,slope = 1)
qp  <-ggplot(data = res_rho_plot, aes(x=true, y=qp)) + xlim(0,0.5) + ylim(0,0.5) +
                geom_point() + geom_abline(intercept = 0,slope = 1)
grid.arrange(EMWeight,EMLaplace, EMHet,EMOri,svr,ls,rls,qp,ncol=4)
dev.off()

pdf(sprintf("../figures/%s/error_%d_%d_%d.pdf",Sys.Date(),cellnum,simsize,simnoise), width = 12, height = 4)
err.mean = melt(res.err)
colnames(err.mean) = c("CellTypes","Methods","Repitition","RMSE")  
err_meanplot <- ggplot(data = err.mean, aes(x=CellTypes, 
                    y=RMSE,fill = Methods))+geom_boxplot()+theme_set(theme_cowplot())
print(err_meanplot)
dev.off()
cellTypes <- cellTypes
cormat <- matrix(NA,nrow = length(cellTypes),ncol = length(methods))
colnames(cormat) <- methods
rownames(cormat) <- cellTypes
for(i in 1:length(cellTypes)){{
  cormat[i,] <- sapply(1:length(methods), FUN = function(j){{
    cor(res.rho[,i,methods[j]],true$mix[,i])
  }})
}}

print("Total Error")
print(colMeans(apply(res.err,c(1,2),mean)))
print(colMeans(apply(res.err,c(1,2),sd)))
print("Correlation")
print(colMeans(apply(res.cor,c(1,2),mean)))
print(colMeans(apply(res.cor,c(1,2),sd)))
print("RSS")
print(colMeans(apply(res.rss,c(1,2),mean)))
print(colMeans(apply(res.rss,c(1,2),sd)))



fmt = "$ \\lambda = %d $   & %1.2e(%1.2e) & %1.2e(%1.2e) & %1.2e(%1.2e) & %1.2e(%1.2e) & %1.2e(%1.2e) & %1.2e(%1.2e) \\\\ "
temp = rbind(colMeans(apply(res.err^2,c(1,2),mean)),colMeans(apply(res.err^2,c(1,2),sd)))
num = as.list(matrix(temp,nrow = 1,byrow = FALSE))
rep = c(simnoise,num)
do.call(sprintf,append(fmt,rep))

#bias of pi & nu0 & bias of sigma_a(=simnoise /cellNum) % sigma_c(=1/cellNum) & FCR & FAR & Iter
fmt2 = "\\multirow2}}{{*}}{{$\\lambda = %d$}}  
& OriEM    & %.3f(%.3f) & %.3f(%.3f)  & %.3f(%.3f)  & %.3f(%.3f)  &%.3f(%.3f) &%.3f(%.3f) & %3.2f(%1.2f)  \\\\
& HeteroEM & %.3f(%.3f) & %.3f(%.3f)  & %.3f(%.3f)  & %.3f(%.3f)  &%.3f(%.3f) &%.3f(%.3f) & %3.2f(%1.2f)  \\\\
& MaxVarEM & %.3f(%.3f) & %.3f(%.3f)  & %.3f(%.3f)  & %.3f(%.3f)  &%.3f(%.3f) &%.3f(%.3f) & %3.2f(%1.2f)  \\\\
& LaplaceEM & %.3f(%.3f) & %.3f(%.3f)  & %.3f(%.3f)  & %.3f(%.3f)  &%.3f(%.3f) &%.3f(%.3f) & %3.2f(%1.2f)  \\\\"
temp2 = rbind(colMeans(abs(res.pi-simpi)),apply(abs(res.pi-simpi),2,sd),
              colMeans(res.nu0),apply(abs(res.nu0),2,sd),
              colMeans(abs(res.sigma_a-simnoise/cellnum)),
              apply(abs(res.sigma_a-simnoise/cellnum),2,sd),
              colMeans(abs(res.sigma_c-1/cellnum)),
              apply(abs(res.sigma_c-1/cellnum),2,sd),
              colMeans(res.confu[,,2]),apply(res.confu[,,2],2,sd),
              colMeans(res.confu[,,3]),apply(res.confu[,,3],2,sd),
              colMeans(res.iter),apply(res.iter,2,sd))
num = as.list(matrix(temp2,nrow = 1,byrow = FALSE))
rep = c(simnoise,num)
do.call(sprintf,append(fmt2,rep))

quit(save = 'no')
