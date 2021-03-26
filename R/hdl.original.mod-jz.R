#### Experimental modifications and additions to Genomic SEM by Johan Zvrskovec 2021

#### GenomicSEM multivariable HDL function, based on the amazing work by Ning, Pawitan and Shen, Nature Genetics (2020)
### This is a wrapper around the original HDL method and needs the original HDL package loaded.
hdl.original<-function(traits,sample.prev=NA,population.prev=NA,trait.names,LD.path,liabilityScale=FALSE){
  
  #for testing
  # traits = project$sumstats.sel$mungedpath
  # sample.prev = project$sumstats.sel$samplePrevalence
  # population.prev = project$sumstats.sel$populationPrevalence
  # trait.names=project$sumstats.sel$code
  # LD.path=project$folderpath.data.HLD.ld
  
  
  n.traits<-length(traits)
  result.S<-matrix(data=NA,nrow = n.traits, ncol = n.traits)
  result.S.se<-matrix(data=NA,nrow = n.traits, ncol = n.traits)
  result.S_std<-matrix(data=NA,nrow = n.traits, ncol = n.traits)
  result.S_std.se<-matrix(data=NA,nrow = n.traits, ncol = n.traits)
  result.P<-matrix(data=NA,nrow = n.traits, ncol = n.traits)
  
  
  #traiti<-2
  #traitj<-3
  for(traiti in 1:n.traits){
    gwas1.df <- as.data.frame(suppressMessages(read_delim(traits[traiti], "\t", escape_double = FALSE, trim_ws = TRUE,progress = F)))
    for(traitj in 1:n.traits){
      
      if(traiti==traitj){
        hdlres.h2<-NA
        hdlres.h2<-HDL.h2(gwas.df = gwas1.df, LD.path = LD.path)
        
        if(is.na(hdlres.h2)) {
          result.S[traiti,traitj]<-NA_real_
          result.S.se[traiti,traitj]<-NA_real_
          result.S_std[traiti,traitj]<-NA_real_
          result.S_std.se[traiti,traitj]<-NA_real_
          result.P[traiti,traitj]<-NA_real_
        } else {
          result.S[traiti,traitj]<-hdlres.h2$h2
          result.S.se[traiti,traitj]<-hdlres.h2$h2.se
          result.S_std.se[traiti,traitj]<-1
          result.P[traiti,traitj]<-hdlres.h2$P
        }
        
      } else {
        gwas2.df <- as.data.frame(suppressMessages(read_delim(traits[traitj], "\t", escape_double = FALSE, trim_ws = TRUE,progress = F)))
        hdlres.rg<-NA
        hdlres.rg<-HDL.rg(gwas1.df = gwas1.df, gwas2.df = gwas2.df, LD.path = LD.path)
        if(is.na(hdlres.rg)) {
          result.S[traiti,traitj]<-NA_real_
          result.S.se[traiti,traitj]<-NA_real_
          result.S_std[traiti,traitj]<-NA_real_
          result.S_std.se[traiti,traitj]<-NA_real_
          result.P[traiti,traitj]<-NA_real_
        } else {
          result.S[traiti,traitj]<-hdlres.rg$estimates.df[c('Genetic_Covariance'),1]
          result.S.se[traiti,traitj]<-hdlres.rg$estimates.df[c('Genetic_Covariance'),2]
          result.S_std[traiti,traitj]<-hdlres.rg$estimates.df[c('Genetic_Correlation'),1]
          result.S_std.se[traiti,traitj]<-hdlres.rg$estimates.df[c('Genetic_Correlation'),2]
          result.P[traiti,traitj]<-hdlres.rg$P
        }
      }
      
      
    }
  }
  
  ##mod-jz experimental liability scale conversion
  if(liabilityScale){
    liab.S <- matrix(1,nrow=1,ncol=n.traits)
    
    for(z in 1:n.traits){
      pop.prev <- population.prev[z]
      samp.prev <- sample.prev[z]
      
      if(is.na(pop.prev)==F & is.na(samp.prev)==F){
        conversion.factor <- (pop.prev^2*(1-pop.prev)^2)/(samp.prev*(1-samp.prev)* dnorm(qnorm(1-pop.prev))^2)
        liab.S[,z] <- conversion.factor
      }}
    
    result.S.observed <- result.S
    result.S_std.observed <- result.S_std
    
    
    ### Scale S (in result.S) to liability:
    result.S <- diag(as.vector(sqrt(liab.S))) %*% result.S %*% diag(as.vector(sqrt(liab.S)))
    
    #calculate the ratio of the rescaled and original S matrices
    scaleO=as.vector(lowerTriangle((result.S/result.S.observed),diag=T))
    
    #rescale the SEs by the same multiples that the S matrix was rescaled by
    #result.S.se<-as.vector(result.S.se*t(scaleO)) #not working
    
    
    ### Scale standardised S (in result.S_std) to liability:
    result.S_std <- diag(as.vector(sqrt(liab.S))) %*% result.S_std %*% diag(as.vector(sqrt(liab.S)))
    
    #calculate the ratio of the rescaled and original S matrices
    scaleO=as.vector(lowerTriangle((result.S_std/result.S_std.observed),diag=T))
    
    #rescale the SEs by the same multiples that the S matrix was rescaled by
    #result.S_std.se<-as.vector(result.S_std.se*t(scaleO)) #not working
    
    
  }
  
  colnames(result.S) <- trait.names
  colnames(result.S.se) <- trait.names
  colnames(result.S_std) <- trait.names
  colnames(result.S_std.se) <- trait.names
  
  return(list(S=result.S, S.se=result.S.se, S_std=result.S_std, S_std.se=result.S_std.se, P=result.P))
}












