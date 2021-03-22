munge.mod <- function(files,hm3,trait.names=NULL,N,info.filter = .9,maf.filter=0.01,path.dir.output="",doChrSplit=FALSE, doStatistics=TRUE){
  
  
  # for testing
  #library(data.table)
  #library(R.utils)
  #files<-c("/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/cleaned/ALCD03.gz")
  #files<-c("/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/reference.1000G.maf.0.005.txt")
  #hm3<-"/Users/jakz/Documents/local_db/JZ_GED_PHD_C1/data/w_hm3.snplist.flaskapp2018"
  #trait.names=c("REF1KG")
  #info.filter=0.6
  #maf.filter=0.05
  #path.dir.output="/Users/jakz/Documents/local_db/JZ_GED_PHD_ADMIN_GENERAL/data/gwas_sumstats/munged"
  
  length <- length(files)
  filenames <- as.vector(files)
  
  log2<-paste(trait.names,collapse="_")
  
  log.file <- file(file.path(path.dir.output,paste0(log2, "_munge.mod.log")),open="wt")
  
  begin.time <- Sys.time()
  
  cat(print(paste0("The munging of ", length(trait.names), " summary statistics started at ",begin.time), sep = ""),file=log.file,sep="\n",append=TRUE)
  
  cat(print(paste("Reading summary statistics for", paste(files,collapse=" "), ". Please note that this step usually takes a few minutes due to the size of summary statistic files.")),file=log.file,sep="\n",append=TRUE)
  
  ##note that fread is not used here due to formatting differences across summary statistic files
  files = lapply(files, read.table,header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
  cat(print("Reading in reference file"),file=log.file,sep="\n",append=TRUE)
  ref <- fread(hm3,header=T,data.table=F)
  cat(print("All files loaded into R!"),file=log.file,sep="\n",append=TRUE)
  
  for(i in 1:length){
    #for testing!
    #i=1
    
    cat(paste("     "),file=log.file,sep="\n",append=TRUE)
    cat(paste("     "),file=log.file,sep="\n",append=TRUE)
    
    cat(print(paste("Munging file:", filenames[i])),file=log.file,sep="\n",append=TRUE)
    hold_names <- toupper(names(files[[i]]))
    
    names1<-hold_names
    if("SNP" %in% hold_names) cat(print(paste("Interpreting the SNP column as the SNP column.")),file=log.file,sep="\n",append=TRUE)
    hold_names[hold_names %in% c("SNP","SNPID","RSID","RS_NUMBER","RS_NUMBERS", "MARKERNAME", "ID","PREDICTOR")] <- "SNP"
    if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the SNP column.")),file=log.file,sep="\n",append=TRUE)
    
    names1<-hold_names
    if("A1" %in% hold_names) cat(print(paste("Interpreting the A1 column as the A1 column.")),file=log.file,sep="\n",append=TRUE)
    hold_names[hold_names %in%c("A1", "ALLELE1","EFFECT_ALLELE","INC_ALLELE","REFERENCE_ALLELE","EA","REF")] <- "A1"
    if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the A1 column.")),file=log.file,sep="\n",append=TRUE)
    
    names1<-hold_names
    if("A2" %in% hold_names) cat(print(paste("Interpreting the A2 column as the A2 column.")),file=log.file,sep="\n",append=TRUE)
    hold_names[hold_names %in%c("A2","ALLELE2","ALLELE0","OTHER_ALLELE","REF","NON_EFFECT_ALLELE","DEC_ALLELE","OA","NEA", "ALT")]  <- "A2"
    if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the A2 column.")),file=log.file,sep="\n",append=TRUE)
    
    names1<-hold_names
    if("EFFECT" %in% hold_names) cat(print(paste("Interpreting the effect column as the effect column.")),file=log.file,sep="\n",append=TRUE)
    hold_names[hold_names %in%c("OR","B","BETA","LOG_ODDS","EFFECTS","EFFECT","SIGNED_SUMSTAT", "Z","ZSCORE","EST","ZSTAT","ZSTATISTIC")] <- "effect"
    if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the effect column.")),file=log.file,sep="\n",append=TRUE)
    
    names1<-hold_names
    if("INFO" %in% hold_names) cat(print(paste("Interpreting the INFO column as the INFO column.")),file=log.file,sep="\n",append=TRUE)
    hold_names[hold_names %in%c("INFO")] <- "INFO"
    if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the INFO column.")),file=log.file,sep="\n",append=TRUE)
    
    names1<-hold_names
    if("P" %in% hold_names) cat(print(paste("Interpreting the P column as the P column.")),file=log.file,sep="\n",append=TRUE)
    hold_names[hold_names %in%c("P","PVALUE","PVAL","P_VALUE","P-VALUE","P.VALUE","P_VAL","GC_PVALUE","WALD_P")] <- "P"
    if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the P column.")),file=log.file,sep="\n",append=TRUE)
    
    names1<-hold_names
    if("N" %in% hold_names) cat(print(paste("Interpreting the N column as the N (sample size) column.")),file=log.file,sep="\n",append=TRUE)
    hold_names[hold_names %in%c("N","WEIGHT","NCOMPLETESAMPLES", "TOTALSAMPLESIZE", "TOTALN", "TOTAL_N","N_COMPLETE_SAMPLES" )] <- "N"
    if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the ", setdiff(names1, hold_names), " column as the N (sample size) column.")),file=log.file,sep="\n",append=TRUE)
    
    names1<-hold_names
    if("N_CAS" %in% hold_names) cat(print(paste("Interpreting the N_CAS column as the N_CAS (sample size for cases) column.")),file=log.file,sep="\n",append=TRUE)
    hold_names[hold_names %in%c("NCASE","N_CASE","N_CASES","N_CAS", "NCAS", "NCA")] <- "N_CAS"
    if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the ", setdiff(names1, hold_names), " column as the N_CAS (sample size for cases) column.")),file=log.file,sep="\n",append=TRUE)
    
    names1<-hold_names
    if("N_CON" %in% hold_names) cat(print(paste("Interpreting the N_CON column as the N_CON (sample size for controls) column.")),file=log.file,sep="\n",append=TRUE)
    hold_names[hold_names %in%c("NCONTROL","N_CONTROL","N_CONTROLS","N_CON","CONTROLS_N", "NCON", "NCO")] <- "N_CON"
    if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the ", setdiff(names1, hold_names), " column as the N_CON (sample size for controls) column.")),file=log.file,sep="\n",append=TRUE)
    
    # Print a message for misisng P value, rs, effect or allele columns
    if(sum(hold_names %in% "P") == 0) cat(print(paste0('Cannot find P-value column, try renaming it to P in the summary statistics file for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(hold_names %in% "A1") == 0) cat(print(paste0('Cannot find effect allele column, try renaming it to A1 in the summary statistics file for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(hold_names %in% "A2") == 0) cat(print(paste0('Cannot find other allele column, try renaming it to A2 in the summary statistics file for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(hold_names %in% "effect") == 0) cat(print(paste0('Cannot find beta or effect column, try renaming it to effect in the summary statistics file for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(hold_names %in% "SNP") == 0) cat(print(paste0('Cannot find rs-id column, try renaming it to SNP in the summary statistics file for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
    
    # Print a warning message when multiple columns interprets as P value, rs, effect or allele columns
    if(sum(hold_names %in% "P") > 1) cat(print(paste0('Multiple columns are being interpreted as the P-value column. Try renaming the column you dont want interpreted as P to P2 for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(hold_names %in% "A1") > 1) cat(print(paste0('Multiple columns are being interpreted as the effect allele column. Try renaming the column you dont want interpreted as effect allele column to A1_2 for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(hold_names %in% "A2") > 1) cat(print(paste0('Multiple columns are being interpreted as the other allele column. Try renaming the column you dont want interpreted as the other allele column to A2_2 for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(hold_names %in% "effect") > 1) cat(print(paste0('Multiple columns are being interpreted as the beta or effect column. Try renaming the column you dont want interpreted as the beta or effect column to effect2 for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
    if(sum(hold_names %in% "SNP") > 1) cat(print(paste0('Multiple columns are being interpreted as the rs-id column. Try renaming the column you dont want interpreted as rs-id to SNP2 for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
    
    # Throw warnings for misisng P valuue, rs, effect or allele columns
    if(sum(hold_names %in% "P") == 0) warning(paste0('Cannot find P-value column, try renaming it P in the summary statistics file for:',filenames[i]))
    if(sum(hold_names %in% "A1") == 0) warning(paste0('Cannot find effect allele column, try renaming it A1 in the summary statistics file for:',filenames[i]))
    if(sum(hold_names %in% "A2") == 0) warning(paste0('Cannot find other allele column, try renaming it A2 in the summary statistics file for:',filenames[i]))
    if(sum(hold_names %in% "effect") == 0) warning(paste0('Cannot find beta or effect column, try renaming it effect in the summary statistics file for:',filenames[i]))
    if(sum(hold_names %in% "SNP") == 0) warning(paste0('Cannot rs-id column, try renaming it SNP in the summary statistics file for:',filenames[i]))
    
    ##rename common MAF labels
    names1<-hold_names
    if("MAF" %in% hold_names) cat(print(paste("Interpreting the MAF column as the MAF (minor allele frequency) column.")),file=log.file,sep="\n",append=TRUE)
    hold_names[hold_names %in% toupper(c("MAF", "CEUaf", "Freq1", "EAF", "Freq1.Hapmap", "FreqAllele1HapMapCEU", "Freq.Allele1.HapMapCEU", "EFFECT_ALLELE_FREQ", "Freq.A1"))] <- "MAF"
    if(length(setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the ", setdiff(names1, hold_names), " column as the MAF (minor allele frequency) column.")),file=log.file,sep="\n",append=TRUE)
    
    #storing the original names
    names1<-names(files[[i]])
    #Replace the origonal names
    names(files[[i]]) <- hold_names
    
    if("MAF" %in% colnames(files[[i]])) {
      ##make sure MAF is actually MAF (i.e., max value is .5 or less)
      if(any(files[[i]]$MAF>1)){
        warning(paste0('Discovered MAF values above 1. This could indicate that the MAF column is not interpreted correctly.',filenames[i]))
      }else if(any(files[[i]]$MAF>1)){
        warning(paste0('Discovered MAF values above .5. Converting these to 1-VALUE.',filenames[i]))
        files[[i]]$MAF<-ifelse(files[[i]]$MAF <= .5, files[[i]]$MAF, (1-files[[i]]$MAF))
      }
    } else {
      #Add MAF here to not misinterpret any additional MAF column as the GWAS MAF
      files[[i]]$MAF<-NA_real_
    }
    
    # Compute N is N cases and N control is reported:
    if("N_CAS" %in% colnames(files[[i]])) {
      files[[i]]$N <- files[[i]]$N_CAS + files[[i]]$N_CON
      cat(print(paste("As the file includes both N_CAS and N_CON columns, the summation of these two columns will be used as the total sample size")),file=log.file,sep="\n",append=TRUE)
    }
    
    ##make sure all alleles are upper case for matching to reference file
    files[[i]]$A1 <- factor(toupper(files[[i]]$A1), c("A", "C", "G", "T"))
    files[[i]]$A2 <- factor(toupper(files[[i]]$A2), c("A", "C", "G", "T"))
    
    ##merge with ref file
    cat(print(paste("Merging file:", filenames[i], "with the reference file:", hm3)),file=log.file,sep="\n",append=TRUE)
    b<-nrow(files[[i]])
    cat(print(paste(b, "rows present in the full", filenames[i], "summary statistics file.")),file=log.file,sep="\n",append=TRUE)
    
    
    ##Addition: produce statistics on non-overlapping SNPs with reference
    if(doStatistics) {
      gwas.snps.not.in.ref_filepath<-file.path(path.dir.output,paste0(trait.names[i],".snps.not.in.ref.txt"))
      gwas.snps.not.in.ref<-merge(ref,files[[i]],by="SNP",all.x=F,all.y=T)
      gwas.snps.not.in.ref<-gwas.snps.not.in.ref[which(is.na(gwas.snps.not.in.ref$A1.x)),c("SNP")]
      write.table(x = gwas.snps.not.in.ref,file = gwas.snps.not.in.ref_filepath,sep="\t", quote = FALSE, col.names=F, row.names = F)
      ref.snps.not.in.gwas_filepath<-file.path(path.dir.output,paste0("ref.snps.not.in.",trait.names[i],".txt"))
      ref.snps.not.in.gwas<-merge(ref,files[[i]],by="SNP",all.x=T,all.y=F)
      ref.snps.not.in.gwas<-ref.snps.not.in.gwas[which(is.na(ref.snps.not.in.gwas$A1.y)),c("SNP")]
      write.table(x = ref.snps.not.in.gwas,file = ref.snps.not.in.gwas_filepath,sep="\t", quote = FALSE, col.names=F, row.names = F)
    }
    
    files[[i]] <- merge(ref,files[[i]],by="SNP",all.x=F,all.y=F)
    cat(print(paste((b-nrow(files[[i]])), "rows were removed from the", filenames[i], "summary statistics file as the rs-ids for these rows were not present in the reference file.")),file=log.file,sep="\n",append=TRUE)
    
    #immediate rename intruding columns from ref
    colnames(files[[i]])[colnames(files[[i]])=="MAF.x"]<-"MAF"
    
    ##remove any rows with missing p-values
    b<-nrow(files[[i]])
    if("P" %in% colnames(files[[i]])) {
      files[[i]]<-subset(files[[i]], !(is.na(files[[i]]$P)))
    }
    if(b-nrow(files[[i]]) > 0) cat(print(paste(b-nrow(files[[i]]), "rows were removed from the", filenames[i], "summary statistics file due to missing values in the P-value column")),file=log.file,sep="\n",append=TRUE)
    
    ##remove any rows with missing effects
    b<-nrow(files[[i]])
    if("effect" %in% colnames(files[[i]])) {
      files[[i]]<-subset(files[[i]], !(is.na(files[[i]]$effect)))
    }
    if(b-nrow(files[[i]]) > 0) cat(print(paste(b-nrow(files[[i]]), "rows were removed from the", filenames[i], "summary statistics file due to missing values in the effect column")),file=log.file,sep="\n",append=TRUE)
    
    ##determine whether it is OR or logistic/continuous effect based on median effect size 
    if("effect" %in% colnames(files[[i]])) {
      a1<-files[[i]]$effect[[1]]
      files[[i]]$effect<-ifelse(rep(round(median(files[[i]]$effect,na.rm=T)) == 1,nrow(files[[i]])), log(files[[i]]$effect),files[[i]]$effect)
      a2<-files[[i]]$effect[[1]]
      if(a1 != a2) cat(print(paste("The effect column was determined to be coded as an odds ratio (OR) for the", filenames[i], "summary statistics file. Please ensure this is correct.")),file=log.file,sep="\n",append=TRUE)
      
      # Flip effect to match ordering in ref file
      files[[i]]$effect<-ifelse(files[[i]]$A1.x != (files[[i]]$A1.y) & files[[i]]$A1.x == (files[[i]]$A2.y),files[[i]]$effect*-1,files[[i]]$effect)
    }
    
    ##remove SNPs that don't match A1 OR A2 in reference file.
    b<-nrow(files[[i]])
    files[[i]]<-subset(files[[i]], !(files[[i]]$A1.x != (files[[i]]$A1.y)  & files[[i]]$A1.x != (files[[i]]$A2.y)))
    if(b-nrow(files[[i]]) > 0) cat(print(paste(b-nrow(files[[i]]), "row(s) were removed from the", filenames[i], "summary statistics file due to the effect allele (A1) column not matching A1 or A2 in the reference file.")),file=log.file,sep="\n",append=TRUE)
    
    b<-nrow(files[[i]])
    files[[i]]<-subset(files[[i]], !(files[[i]]$A2.x != (files[[i]]$A2.y)  & files[[i]]$A2.x !=  (files[[i]]$A1.y)))
    if(b-nrow(files[[i]]) > 0) cat(print(paste(b-nrow(files[[i]]), "row(s) were removed from the", filenames[i], "summary statistics file due to the other allele (A2) column not matching A1 or A2 in the reference file.")),file=log.file,sep="\n",append=TRUE)
    
    ####VALIDITY CHECKS#####
    
    #Check that p-value column does not contain an excess of 1s/0s
    if((sum(files[[i]]$P > 1) + sum(files[[i]]$P < 0)) > 100){
      cat(print("In excess of 100 SNPs have P val above 1 or below 0. The P column may be mislabled!"),file=log.file,sep="\n",append=TRUE)
    }
    
    #Compute Z score
    if("effect" %in% colnames(files[[i]])) {
      files[[i]]$Z <- sign(files[[i]]$effect) * sqrt(qchisq(files[[i]]$P,1,lower=F))
    }
    
    ##filter on INFO column at designated threshold provided for the info.filter argument (default = 0.9)
    if("INFO" %in% colnames(files[[i]])) {
      b<-nrow(files[[i]])
      files[[i]] <- files[[i]][files[[i]]$INFO >= info.filter,]
      cat(print(paste(b-nrow(files[[i]]), "rows were removed from the", filenames[i], "summary statistics file due to INFO values below the designated threshold of", info.filter)),file=log.file,sep="\n",append=TRUE)
    }else{cat(print("No INFO column, cannot filter on INFO, which may influence results"),file=log.file,sep="\n",append=TRUE)}
    
    ##filter on MAF filter at designated threshold provided for the maf.filter argument (default = 0.01)
    if("MAF" %in% colnames(files[[i]])) {
      files[[i]]$MAF<-as.numeric(as.character(files[[i]]$MAF))
      b<-nrow(files[[i]])
      files[[i]] <- files[[i]][files[[i]]$MAF >= maf.filter,]
      files[[i]]<-subset(files[[i]], !(is.na(files[[i]]$MAF)))
      cat(print(paste(b-nrow(files[[i]]), "rows were removed from the", filenames[i], "summary statistics file due to missing MAF information or MAFs below the designated threshold of", maf.filter)),file=log.file,sep="\n",append=TRUE)
    }else{cat(print("No MAF column, cannot filter on MAF, which may influence results"),file=log.file,sep="\n",append=TRUE)}
    
    if(!("N" %in% colnames(files[[i]]))) {
      files[[i]]$N<-N[i]
    }
    
    #final rename columns
    colnames(files[[i]])[colnames(files[[i]])=="A1.x"]<-"A1"
    colnames(files[[i]])[colnames(files[[i]])=="A2.x"]<-"A2"
    
    colnames(files[[i]])[colnames(files[[i]])=="CHR.x"]<-"CHR"
    
    #final remove columns
    files[[i]] <- files[[i]][,!(colnames(files[[i]]) %in% c("effect","A1.y","A2.y","CHR.y","MAF.y"))]
    
    #output.colnames<- c("SNP","N","Z","A1","A2")
    output.colnames<- c("SNP")
    if("N" %in% colnames(files[[i]])) {
      output.colnames<- c(output.colnames,"N")
    }
    if("Z" %in% colnames(files[[i]])) {
      output.colnames<- c(output.colnames,"Z")
    }
    output.colnames<- c(output.colnames,c("A1","A2"))
    output.colnames.more<-colnames(files[[i]])[!(colnames(files[[i]]) %in% output.colnames)]
    
    
    output<-files[[i]][,c(output.colnames,output.colnames.more)]
    #original output
    #output <- cbind.data.frame(files[[i]]$SNP,files[[i]]$N,files[[i]]$Z,files[[i]]$A1.x,files[[i]]$A2.x)
    
    if(!("N" %in% names(files[[i]])) & (exists("N") == FALSE)) cat(warning(paste0('Cannot find sample size column for',filenames[i], " and a sample size was not provided for the N argument. Please either provide a total sample size to the N argument or try changing the name of the sample size column to N.")),file=log.file,sep="\n",append=TRUE)
    
    cat(print(paste(nrow(output), "SNPs are left in the summary statistics file", filenames[i], "after QC.")),file=log.file,sep="\n",append=TRUE)
    
    #remove spaces in trait.names file to avoid errors with fread functionality used for s_ldsc
    trait.names[i]<-str_replace_all(trait.names[i], fixed(" "), "") 
    
    nfilepath<-file.path(path.dir.output,trait.names[i])
    write.table(x = output,file = nfilepath,sep="\t", quote = FALSE, row.names = F)
    nfilepath.gzip<-gzip(nfilepath)
    cat(print(paste("I am done munging file:", filenames[i])),file=log.file,sep="\n",append=TRUE)
    cat(print(paste("The file is saved as", nfilepath.gzip, "in the current working directory.")),file=log.file,sep="\n",append=TRUE)
    
    #addition: producing per-chromosome files in a folder, as RAISS columns
    if(doChrSplit) {
      if("CHR" %in% colnames(files[[i]])){
        dir.create(paste0(nfilepath,".chr"), showWarnings = FALSE)
        validChromosomes<-c(1:22,"X","Y","XY","MT") #as per Plink standard
        for(chr in validChromosomes){
          output.chr<-output[which(output$CHR==chr),c("SNP","ORIGBP","A1","A2","Z")]
          colnames(output.chr)<-c("rsID","pos","A0","A1","Z")
          write.table(x = output.chr,file = file.path(paste0(nfilepath,".chr"), chr),sep="\t", quote = FALSE, row.names = F)
        }
      }
    }
    
  }
  
  end.time <- Sys.time()
  
  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- total.time-mins*60
  
  cat(paste("     "),file=log.file,sep="\n",append=TRUE)
  cat(print(paste0("Munging was completed at ",end.time), sep = ""),file=log.file,sep="\n",append=TRUE)
  
  cat(print(paste0("The munging of all files took ",mins," minutes and ",secs," seconds"), sep = ""),file=log.file,sep="\n",append=TRUE)
  cat(print(paste("Please check the log file", paste0(log2, "_munge.log"), "to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files")),file=log.file,sep="\n",append=TRUE)
  
  flush(log.file)
  close(log.file)
  
}
