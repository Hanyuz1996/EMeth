filterGene <- function(dat,pv,level,label){
  num_of_cell <- length(unique(label))
  num_of_level <- length(level)
  num_of_gene <- dim(dat)[1]
  genes <- list()
  for(i in 1:num_of_level){
    nowcomp <- level[[i]]
    genei <- list()
    datt1 <- subset(dat,select = which(label %in% nowcomp[[1]]))
    #print('datt1')
    #print(datt1[1:5,])
    datt2 <- subset(dat,select = which(label %in% nowcomp[[2]]))
    #print('datt2')
    #print(datt2[1:5,])
    fc <- rep(0,dim(datt1)[1])
    pvv <- rep(-1,dim(datt1)[1])
    for(j in 1:num_of_gene){
      samp1 <- as.matrix(datt1[j,])
      #print(length(samp1))
      #print(samp1)
      samp2 <- as.matrix(datt2[j,])
      #print(length(samp2))
      #print(samp2)
      #cat('mean1',mean(samp1,na.rm=TRUE),'mean2',mean(samp2,na.rm=TRUE),'\n')
      if(sd(samp1)<0.00001 & sd(samp2) < 0.0001){next}
      tt <- t.test(samp1,samp2)
      
      if(tt$p.value < pv){
        genei <- c(genei, j)
        pvv[j] <- tt$p.value
        fc[j] <- log10(mean(samp1,na.rm = TRUE)/mean(samp2,na.rm = TRUE))
      }
    }
    pvvpos <- -log10(1/pvv[unlist(genei)])
    fcpos  <- fc[unlist(genei)]
    print(length(genei))
    
    pdf(sprintf('~/Hutch-Research/figures/Chen/gene_%s_vs_%s.pdf',nowcomp[[1]],do.call(paste,c(as.list(nowcomp[[2]]),sep="+"))))
    tempdata <- as.data.frame(cbind(pvvpos,fcpos))
    colnames(tempdata) <- c("Log10PValue","LogFoldChange")
    newplot <- ggplot(tempdata,aes(y=Log10PValue,x=LogFoldChange))+geom_point()+
               ggtitle(sprintf('%s_vs_%s',nowcomp[[1]],do.call(paste,c(as.list(nowcomp[[2]],sep='+')))))
    print(newplot)
    dev.off()
    
    len <- min(100,length(genei))
    fc <- abs(fc)
    fcind <- order(fc,decreasing = TRUE)
    genei <- intersect(unlist(genei),fcind[1:len])
    genes <- c(genes,genei)
  }
  return(genes)
}
