# Functions for testing functions in IFFMathR/computeNonBinomialRobust.R
## Copyright 2011 Webtrends Inc. All Rights Reserved.
## Webtrends PROPRIETARY/CONFIDENTIAL. Use is subject to license terms.
## Functions in this file:


#source(file.path(Sys.getenv("R_USER"),"auxiliaryLib.R"))
#source(file.path(Sys.getenv("R_USER"),"IFFMathR","experimentDesignInteraction.R"))
#source(file.path(Sys.getenv("R_USER"),"IFFMathR","computeNonBinomialRobust.R"))
#source(file.path(Sys.getenv("R_USER"),"IFFMathR","callComputeAnovaTableInteraction.R"))
require(Rfit)

call.compute.NonBinomial.AllSaints <- function(plots=FALSE)
## Testing the East Coast data. See email from Phil on 6/12/12.
## 6/15/2012
{
printout("Testing the All Saints data")
    
    indir = file.path(Sys.getenv("R_DATA"),"AllSaints")
    file.csv = file.path(indir,"6-12-12","UKDat_SumOfOrderValue-530-601.csv")
    cnames.file=scan(file = file.csv, what = "character", nmax = -1, sep = ",",skip = 0, nlines = 1,quiet = TRUE) #read column names from the file
    
    factor.col = c("experimentName") # columns factors
    numeric.col = c("SumOfOrderValues") # numeric columns
    
    factor.col.ind = is.element(cnames.file,factor.col)
    numeric.col.ind = is.element(cnames.file,numeric.col)

    col.classes = vector("character",length=length(cnames.file)) # column classes
    col.classes[factor.col.ind] <- "factor"
    col.classes[numeric.col.ind] <- "numeric"
    dat = read.csv(file.csv,header=TRUE,colClasses=col.classes)

    fnames = c("experimentName","SumOfOrderValues") # ExperimentName and SumOfOrderValues
    colnames(dat) = fnames
    nf = length(fnames) # number of factors
    N = nrow(dat) # total number of observations
    
    treat = fnames[1] # treatment name
printout(treat)    
    lev = factor(sort(unique(dat[,treat])))# treatment levels
    n.treat = length(lev) # number of treatments
    n.i.treat = sapply(lev,function(x) sum(is.element(dat[,treat],x))) # number of replications for each treatment
    dat.exp.list = lapply(lev,function(x) dat$SumOfOrderValues[is.element(dat[,treat],x)]) # list dat organized by experiments
    names(dat.exp.list) = lev

 printout(head(dat))
 printout(lev)
 printout(n.i.treat)
 
 ### Debugging using the Mucociliary efficiency data, Hollander and Wolfe, p.116
#    X=c(2.9,3,2.5,2.6,3.2)
#   Y=c(3.8,2.7,4,2.4)
#   Z=c(2.8,3.4,3.7,2.2,2)
#   
#    dat.exp.list=list(E1=X,E2=Y,E3=Z)
    
#printout(dat.exp.list)
#

   r = rank(unlist(dat.exp.list),ties.method="average") # vector of ranks in the original order, i.e. (x,y)
printout(head(r))
    numn = sapply(dat.exp.list,length)
printout(numn)
    
    n.treat = length(dat.exp.list) # number of treatments (experiments)
    
    g = unlist(sapply(1:length(numn),function(ind) rep(names(dat.exp.list)[ind],numn[ind])))
    rank.list = split(r,g)
    rank.sum = lapply(rank.list,sum)

    ## Compute ties
    ties.size = table(r)[table(r)>1]

    if(length(ties.size) == 0)
    {
         ties.size = rep(1,sum(numn))   
    }
    

  res=kruskal.test(dat.exp.list)
printout(res) 
#printout(rank.sum)
#printout(numn)
#printout(ties.size)
  
  res = compute.ABn.robust(ranks=rank.sum,numn=numn,ties=ties.size,alpha=0.05,by.pvalue=FALSE)
#printout(res)

    ## Compute estimates of the location parameters, see Step 5 in the RobustAB doc
    est.list = vector("list",n.treat)
    for(e in 1:n.treat) # for each experiment
    {
        obs = dat.exp.list[[e]]
        est.list[[e]] = compute.HLestimates(obs,alpha=0.05)

    }
printout(est.list)
stop("debug")
    
    ## Contrasts
    #printout(contrasts(dat$treat))
    #contrasts(dat$treat) <- contr.sum
    #printout(contrasts(dat[,treat]))

    kpi = tail(fnames,1) # name of the kpi (last column)
    
    
    ## Graphical exploration of the data
    form = SumOfOrderValues ~ experimentName
    
    if(plots)
    {
        boxplot(form,dat) # boxplot, 
    
        for(i in 1:length(lev))
        {
            if(dev.cur() >= 1) dev.new()
            hist(dat[is.element(dat[,1],lev[i]),2],freq=F,main = paste("Histogram of" , "AverageRevenue"),xlab="AverageRevenue")
        }
    }
    ## Rfit
    fit <- rfit(form,dat,intercept=FALSE) # R fit
    
    #if(dev.cur() >= 1) dev.new()
    #qqnorm(rstudent(rfit(ldl~treat,quail))) # studentized residual normal q-q plot, see Fig. 4.2.1, p. 296
    
#    mf <- model.frame(form,dat)
#    #printout(mt)
#    mt <- attr(mf, "terms")
#    attr(mt, "intercept") <- 0
#    x <- model.matrix(mt, data = mf)
##    printout(x)
#    
    printout("Coefficients")
    printout(fit$coefficients)
    printout("tauhat")
    printout(fit$tauhat)
    printout("taushat")
    printout(fit$taushat)
    printout("betahat")
    printout(fit$betahat)
#    printout("fitted.values")
#    printout(fit$fitted.values)
    printout("median of the residuals")
    printout(median(fit$residuals))
    
    printout("Standard Error")
    se = sqrt(n.i.treat^(-1)*fit$tauhat^2+(fit$taushat^2-fit$tauhat^2)/N)
    #printout(se)
    
    fit.df = data.frame(fit$coefficients,se)
    colnames(fit.df) = c("Coeff","Standard Error")
    rownames(fit.df) = names(fit$coefficients)
    printout(fit.df)
    
#    printout("One way anova")
#    fit.one <- oneway.rfit(dat$ldl,dat$treat)
#    printout(fit.one$est)
#    printout(fit.one$se)
#    printout(fit.one$p.value)
#    printout(summary(fit.one))
    
    ## LS fit
    printout("LS fit")
    form.ls = prints_per_session ~ experimentName - 1
    fit.ls = lm(form.ls,dat)
    printout(summary(fit.ls))
    
    if(dev.cur() >= 1) dev.new()
    qqnorm(fit.ls$residuals)
    
    ## Rank-Based Tests of H0:mu1=...=muk, see p. 296 in Hettmansperger and McKean, 2011
    printout("Drop dispersion test")
    test.drop = drop.test(fit) # Rfit package
    test.drop.out = data.frame(test.drop$F,test.drop$tauhat,test.drop$df1,test.drop$df2,test.drop$p.value)
    colnames(test.drop.out) = c("Test Statistics","Scale","df1","df2","p-value")
    printout(test.drop.out)
    #printout(test.drop$F)
    #printout(test.drop$tauhat)
    #printout(c(test.drop$df1,test.drop$df2))
    #printout(test.drop$p.value)
#fit$residuals
}


call.compute.NonBinomial.EastCoast <- function(plots=FALSE,dataset=4)
  ## Testing the East Coast data. See email from Nura on 6/12/12.
  ## 6/12/2012
{
  printout("Testing the East Coast data")
  
  indir = file.path(Sys.getenv("R_DATA"),"EastCoast")
  
  if(dataset==1)
  {
    file.csv = file.path(indir,"TotalRevenue_061212.csv")
    cnames.file=scan(file = file.csv, what = "character", nmax = -1, sep = ",",skip = 0, nlines = 1,quiet = TRUE) #read column names from the file
    
    factor.col = c("experimentName","experimentID","conversionPoint") # columns factors
    numeric.col = c("TotalRevenue","Passengers") # numeric columns
    
    factor.col.ind = is.element(cnames.file,factor.col)
    numeric.col.ind = is.element(cnames.file,numeric.col)
    
    col.classes = vector("character",length=length(cnames.file)) # column classes
    col.classes[factor.col.ind] <- "factor"
    col.classes[numeric.col.ind] <- "numeric"
    dat1 = read.csv(file.csv,header=TRUE,colClasses=col.classes)
    
    dat = data.frame(dat1$experimentName,round(dat1$TotalRevenue/dat1$Passengers,0))
  }
  
  if(dataset == 2) ## see email from Tom 7/16/12
  {
    
    file.csv = file.path(indir,"RT13119_part2.csv")
printout(paste("read file",file.csv))
    cnames.file=scan(file = file.csv, what = "character", nmax = -1, sep = ",",skip = 0, nlines = 1,quiet = TRUE) #read column names from the file
    
    factor.col = c("conversionPoint","Passengers") # columns factors
    numeric.col = c("experimentID","UID","TotalRevenue","Upgrades") # numeric columns
    factor.col.ind = is.element(cnames.file,factor.col)
    numeric.col.ind = is.element(cnames.file,numeric.col)
    
    col.classes = vector("character",length=length(cnames.file)) # column classes
    col.classes[factor.col.ind] <- "factor"
    col.classes[numeric.col.ind] <- "numeric"
    
    dat.all = read.csv(file.csv,header=TRUE,colClasses=col.classes,na.strings = "NULL")
    
    # Need to compute the sum of total revenue for each UID
    
    
    exp.names = sort(unique(dat.all$experimentID))
    ##dat.exp1 = dat.all[dat.all$experimentID == exp.names[1], ] # Data for experiment 1 only
    ##dat.exp2 = dat.all[dat.all$experimentID == exp.names[2], ] # Data for experiment 2 only
    
#     for(exp in 1:2)
#     {
#       dat = dat.all[dat.all$experimentID == exp.names[exp], ]
#       
#       UID.unique = unique(dat$UID)
#       treat.vec = sapply(UID.unique,function(x) sum(dat$TotalRevenue[is.element(dat$UID,x)]) ) # sum of total revenue
#       
#       sumrev.df = data.frame(experimentID=UID.unique,SumTotalRevenue=treat.vec)
#       write.csv(sumrev.df,file=file.path(indir,paste("RT13119_part2_SumRevenue_",exp,".csv")))
#       
#       
#     }
#     
    
    
  } # end if(dataset == 2)
  
  if(dataset == 3) ## see email from Tom 7/16/12
  {
    file.csv = file.path(indir,"exp2.csv")
    cnames.file=scan(file = file.csv, what = "character", nmax = -1, sep = ",",skip = 0, nlines = 1,quiet = TRUE) #read column names from the file
    
    col.classes <- rep("numeric",length(cnames.file))
    dat = read.csv(file.csv,header=TRUE,colClasses=col.classes)
    
    treat.vec = dat$Sum.of.TotalRevenue  # Sum of TotalRevenue
    printout(length(treat.vec))
ptt <- proc.time()  
    out = sort.pair(treat.vec,ind=289859138,descreasing=FALSE)
    #out = compute.HLestimates(treat.vec,alpha=0.05,fast=TRUE)
printout(proc.time()-ptt)
printout(out)

stop("debug")
    
    #treat.vec.unique = sort(unique(treat.vec)) # sorted unique vector
    #printout(head(treat.vec.unique))
    
    #num.unique = sapply(treat.vec.unique,function(x) sum(treat.vec == x)) # number of each uniques in the original vector
    
    printout(head(num.unique))        
  }
  
  if(dataset == 4) ## see email from Tom 7/16/12
  {
    file.csv = file.path(indir,"RT13119_part2_SumRevenue_2.csv")
    cnames.file=scan(file = file.csv, what = "character", nmax = -1, sep = ",",skip = 0, nlines = 1,quiet = TRUE) #read column names from the file
 
    col.classes <- c("character", "numeric","numeric")
    dat = read.csv(file.csv,header=TRUE,row.names=1,colClasses=col.classes)
    
    treat.name = colnames(dat)[length(colnames(dat))]
    treat.vec = as.numeric(dat[,treat.name])  # Sum of TotalRevenue
printout(head(treat.vec))
printout(median(treat.vec))
    ptt <- proc.time()  
    out = sort.pair(treat.vec,ind=289859138,descreasing=FALSE)
    #out = compute.HLestimates(treat.vec,alpha=0.05,fast=TRUE)
    printout(proc.time()-ptt)
    printout(out)
    return(out)
  } # end if(dataset == 4)
  
  
  stop("debug")
  
  fnames = c("experimentName","AverageRevenue") # ExperimentName and TotalRevenue
  colnames(dat) = fnames
  nf = length(fnames) # number of factors
  N = nrow(dat) # total number of observations
  
  treat = fnames[1] # treatment name
  printout(treat)    
  lev = factor(unique(dat[,treat]),levels=paste("Experiment",1:16)) # treatment levels
  n.treat = length(lev) # number of treatments
  n.i.treat = sapply(lev,function(x) sum(is.element(dat[,treat],x))) # number of replications for each treatment
  
  write.csv(dat,file=file.path(indir,"AverageRevenue_R.csv"))
  
  printout(head(dat))
  printout(lev)
  printout(n.i.treat)
  
  
  ## Contrasts
  #printout(contrasts(dat$treat))
  #contrasts(dat$treat) <- contr.sum
  #printout(contrasts(dat[,treat]))
  
  kpi = tail(fnames,1) # name of the kpi (last column)
  
  
  ## Graphical exploration of the data
  form = AverageRevenue ~ experimentName
  
  if(plots)
  {
    boxplot(form,dat) # boxplot, 
    
    for(i in 1:length(lev))
    {
      if(dev.cur() >= 1) dev.new()
      hist(dat[is.element(dat[,1],lev[i]),2],freq=F,main = paste("Histogram of" , "AverageRevenue"),xlab="AverageRevenue")
    }
  }
  ## Rfit
  fit <- rfit(form,dat,intercept=FALSE) # R fit
  
  #if(dev.cur() >= 1) dev.new()
  #qqnorm(rstudent(rfit(ldl~treat,quail))) # studentized residual normal q-q plot, see Fig. 4.2.1, p. 296
  
  #    mf <- model.frame(form,dat)
  #    #printout(mt)
  #    mt <- attr(mf, "terms")
  #    attr(mt, "intercept") <- 0
  #    x <- model.matrix(mt, data = mf)
  ##    printout(x)
  #    
  printout("Coefficients")
  printout(fit$coefficients)
  printout("tauhat")
  printout(fit$tauhat)
  printout("taushat")
  printout(fit$taushat)
  printout("betahat")
  printout(fit$betahat)
  #    printout("fitted.values")
  #    printout(fit$fitted.values)
  printout("median of the residuals")
  printout(median(fit$residuals))
  
  printout("Standard Error")
  se = sqrt(n.i.treat^(-1)*fit$tauhat^2+(fit$taushat^2-fit$tauhat^2)/N)
  #printout(se)
  
  fit.df = data.frame(fit$coefficients,se)
  colnames(fit.df) = c("Coeff","Standard Error")
  rownames(fit.df) = names(fit$coefficients)
  printout(fit.df)
  
  #    printout("One way anova")
  #    fit.one <- oneway.rfit(dat$ldl,dat$treat)
  #    printout(fit.one$est)
  #    printout(fit.one$se)
  #    printout(fit.one$p.value)
  #    printout(summary(fit.one))
  
  ## LS fit
  printout("LS fit")
  form.ls = prints_per_session ~ experimentName - 1
  fit.ls = lm(form.ls,dat)
  printout(summary(fit.ls))
  
  if(dev.cur() >= 1) dev.new()
  qqnorm(fit.ls$residuals)
  
  ## Rank-Based Tests of H0:mu1=...=muk, see p. 296 in Hettmansperger and McKean, 2011
  printout("Drop dispersion test")
  test.drop = drop.test(fit) # Rfit package
  test.drop.out = data.frame(test.drop$F,test.drop$tauhat,test.drop$df1,test.drop$df2,test.drop$p.value)
  colnames(test.drop.out) = c("Test Statistics","Scale","df1","df2","p-value")
  printout(test.drop.out)
  #printout(test.drop$F)
  #printout(test.drop$tauhat)
  #printout(c(test.drop$df1,test.drop$df2))
  #printout(test.drop$p.value)
  #fit$residuals
}



call.compute.NonBinomial.Coupon <- function(plots=FALSE)
## Testing the Coupon data
## 3/26/2012
{
printout("Testing the Coupon data")
    
    indir = file.path(Sys.getenv("R_DATA"),"Coupons")
    file.csv = file.path(indir,"coupons_simple.csv")
    cnames.file=scan(file = file.csv, what = "character", nmax = -1, sep = ",",skip = 0, nlines = 1,quiet = TRUE) #read column names from the file
    
    col.classes = c(rep("factor",1),"numeric") # column classes
    dat = read.csv(file.csv,header=TRUE,colClasses=col.classes)
    
    ## if 7way: file format with experiments as columns
    #file.csv = file.path(indir,"coupons_7waysplit_data.csv")
    file.csv = file.path(indir,"Coupons_Print_Device_prints_per_session.csv")
    dat.file = read.csv(file.csv,header=TRUE,colClasses="numeric")
    exp.names = colnames(dat.file)
    n.exp = apply(dat.file,2,function(x) sum(!is.na(x)))
    dat.vec = unlist(sapply(1:length(n.exp),function(x) dat.file[1:n.exp[x],x]))
    dat.exp = unlist(sapply(1:length(n.exp),function(x) rep(exp.names[x],n.exp[x])))
    dat = data.frame(dat.exp,dat.vec)
    colnames(dat)=cnames.file
printout(dim(dat))    
    #dat = dat[dat[,2]<2,]
    #dat[,2] = log(dat[,2]) 
#printout(dim(dat))    
#printout(tail(dat.file))
#printout(n.exp)
#printout(head(dat))
#printout(sapply(dat,data.class))
#stop("debug")
    fnames = colnames(dat)
    nf = length(fnames) # number of factors
    N = nrow(dat) # total number of observations
    
    treat = fnames[1] # treatment name
printout(treat)    
    lev = unique(dat[,treat]) # treatment levels
    n.treat = length(lev) # number of treatments
    n.i.treat = sapply(lev,function(x) sum(is.element(dat[,treat],x))) # number of replications for each treatment
 
 printout(head(dat))
 printout(lev)
 printout(n.i.treat)
    
    ## Contrasts
    #printout(contrasts(dat$treat))
    #contrasts(dat$treat) <- contr.sum
    #printout(contrasts(dat[,treat]))

    kpi = tail(fnames,1) # name of the kpi (last column)
 
    #col.classes = c(rep("factor",nf),"numeric") # column classes
    
    ## Graphical exploration of the data
    form = prints_per_session ~ experimentName
    
    if(plots)
    {
        boxplot(form,dat) # boxplot, 
    
        for(i in 1:length(lev))
        {
            if(dev.cur() >= 1) dev.new()
            hist(dat[is.element(dat[,1],lev[i]),2],freq=F,main = paste("Histogram of" , "prints_per_session"),xlab="prints_per_session")
        }
    }
    ## Rfit
    fit <- rfit(form,dat,intercept=FALSE) # R fit
    
    #if(dev.cur() >= 1) dev.new()
    #qqnorm(rstudent(rfit(ldl~treat,quail))) # studentized residual normal q-q plot, see Fig. 4.2.1, p. 296
    
#    mf <- model.frame(form,dat)
#    #printout(mt)
#    mt <- attr(mf, "terms")
#    attr(mt, "intercept") <- 0
#    x <- model.matrix(mt, data = mf)
##    printout(x)
#    
    printout("Coefficients")
    printout(fit$coefficients)
    printout("tauhat")
    printout(fit$tauhat)
    printout("taushat")
    printout(fit$taushat)
    printout("betahat")
    printout(fit$betahat)
#    printout("fitted.values")
#    printout(fit$fitted.values)
    printout("median of the residuals")
    printout(median(fit$residuals))
    
    printout("Standard Error")
    se = sqrt(n.i.treat^(-1)*fit$tauhat^2+(fit$taushat^2-fit$tauhat^2)/N)
    #printout(se)
    
    fit.df = data.frame(fit$coefficients,se)
    colnames(fit.df) = c("Coeff","Standard Error")
    rownames(fit.df) = names(fit$coefficients)
    printout(fit.df)
    
#    printout("One way anova")
#    fit.one <- oneway.rfit(dat$ldl,dat$treat)
#    printout(fit.one$est)
#    printout(fit.one$se)
#    printout(fit.one$p.value)
#    printout(summary(fit.one))
    
    ## LS fit
    printout("LS fit")
    form.ls = prints_per_session ~ experimentName - 1
    fit.ls = lm(form.ls,dat)
    printout(summary(fit.ls))
    
    if(dev.cur() >= 1) dev.new()
    qqnorm(fit.ls$residuals)
    
    ## Rank-Based Tests of H0:mu1=...=muk, see p. 296 in Hettmansperger and McKean, 2011
    printout("Drop dispersion test")
    test.drop = drop.test(fit) # Rfit package
    test.drop.out = data.frame(test.drop$F,test.drop$tauhat,test.drop$df1,test.drop$df2,test.drop$p.value)
    colnames(test.drop.out) = c("Test Statistics","Scale","df1","df2","p-value")
    printout(test.drop.out)
    #printout(test.drop$F)
    #printout(test.drop$tauhat)
    #printout(c(test.drop$df1,test.drop$df2))
    #printout(test.drop$p.value)
#fit$residuals
}


#################################################

# 
#   obs=exp1[,3]
#   
#   res.full=compute.HLestimates(obs,alpha=0.05,fast=TRUE) 
#   
#   n=length(exp1[,3])
#   
#   n.25=floor(n*0.25)
# 
#   n.50=floor(n*0.5)
# 
#   n.75=floor(n*0.75)
# 
#  
#   res.sum = c(0,0,0)
# 
#  count.index = 0
#   
# #for(i in 10)
# 
# n.iter = 10
# 
# n.sample = n.75
# 
# while(count.index<n.iter)
#  {
#   
#   obs.sample=sample(obs) 
#   
#   obs.sub=obs.sample[1:n.sample]
#   
#   res.HLestimates = compute.HLestimates(obs.sub,alpha=0.05,fast=TRUE)
#   
#   res = c(res.HLestimates[[2]][1],res.HLestimates[[1]],res.HLestimates[[2]][2])
#   
#   printout(paste("res.HLestimates"))
#   
#   printout(res)
#   
#   res.sum = res.sum + res
#   
#   count.index = count.index + 1
#     
#  }
#   
#   res = res.sum/n.iter
# 
# 
# 
#   
# #### use normal approximations in wilcox test 
# 
# #obs = rnorm(10000)
# 
# #obs = rt(10000, 3, 0)
# 
# obs = rf(10000, 8, 3, 0)
# 
# #obs = rexp(10000, rate = 1)
# 
# n=length(obs)
# 
# 
# #res.sum = c(0,0,0)
# 
# count.index = 0
# 
# #for(i in 10)
# 
# n.iter = 10
# 
# res=matrix(0,nrow=n.iter,ncol=3)
# 
# while(count.index<n.iter)
# {
#   
#   wilcox.test.res = wilcox.test(obs,conf.int = T)
#   
#   res[count.index+1,] = c(wilcox.test.res$conf.int[1],wilcox.test.res$estimate,wilcox.test.res$conf.int[2])
#   
# #  printout(paste("wilcox.test.res"))
#   
# #  printout(wilcox.test.res)
#   
#   count.index = count.index + 1
#   
# }
# 
# #res = res.sum/n.iter
# 
# LC = mean(res[,1])
# 
# median = mean(res[,2])
# 
# UC = mean(res[,3])
# 
# 
# printout(paste("The mean result of the normal approximation for Wilcox test"))
# 
# printout(c(LC,median,UC))
# 

###############################################################
# 
# res=matrix(0,nrow=n.iter,ncol=3)
# 
#   res.HLestimates = compute.HLestimates(obs.sub,alpha=0.05,fast=TRUE)
#   
#   res = c(res.HLestimates[[2]][1],res.HLestimates[[1]],res.HLestimates[[2]][2])
#   
# 
# printout(paste("The result of HLestimates"))
# 
# printout(res)
# 
#   


### repeated the Boxcox data several times until it reaches a large dataset, for testing the results of Ranova vs Anova of psuedo observatinos
### oroginal Boxcox data: logSurv Poison Treatment, n=48

call.compute.NonBinomial.BoxCox <- function()
{
#  set.seed(123)
#  edf <- data.frame(sex = c(rep("Male", 300), rep("Female", 300), rep("Unknown", 300)),
#  head_length = c(1.2 * c(170:179 + rnorm(300)), 0.8 * c(150:159 + rnorm(300)), c(160:169 + rnorm(300)))/10,
#  body_length = c(c(170:179 + rnorm(300)), c(150:159 + rnorm(300)), c(160:169 + rnorm(300)))) 
#  edf$sex <- factor(as.character(edf$sex))
  
  printout("Testing the BoxCox Data.")
  indir = file.path("X:/projects","Boxcox")
 # file.csv = file.path(indir,"motors.csv")
  inter.array.file = file.path(indir,"interaction_array.csv")
  #cnames.file=scan(file = file.csv, what = "character", nmax = -1, sep = ",",skip = 0, nlines = 1,quiet = TRUE) #read column names from the file
  
  #fnames = head(cnames.file,-1)
  #nf = length(fnames) # number of factors
  
  #kpi = tail(cnames.file,1) # name of the kpi (last column)
  
  #col.classes = c(rep("factor",nf),"numeric") # column classes
  
  ## Read the data 
  #dat.file = read.csv(file.csv,header=TRUE,colClasses=col.classes)
  #printout(head(dat.file))
  #printout(dat.file[1,1:nf]) 

  data(BoxCox)
  printout(BoxCox)
  ### repeat the block 20 times, so n=48*20=960
  BoxCox1<-as.matrix(BoxCox[,2:3])
  
  rBoxCox<-NULL
  for(i in 1:20)
  {
    rBoxCox<-rbind(rBoxCox,BoxCox1)
  }
  
  
  n<-dim(rBoxCox)[1]
  
  Surv<-rnorm(n, mean = 5, sd = 2)
  rBoxCox<-cbind(Surv,rBoxCox)
  
  ### add some outliers
  #outlier.index<-rbern(n, prob=0.05)
  outlier.index<-rbinom(n, size=1, prob=0.05)
  outlier.size<-rnorm(n, mean = 3, sd = 2)
  outlier<-outlier.index*outlier.size
  rBoxCox[,1]<-Surv+outlier

  dat.file<-rBoxCox
  
  nf = dim(rBoxCox)[2] # number of factors
  nr =dim(rBoxCox)[1] # number of rows
  nc = dim(rBoxCox)[2] # number of columns
  lev = as.numeric(apply(rBoxCox[,2:nf],2,max))
  
  #lograte = log(dat.file[,nc]) # log of the surv
  #dat.file[,nc] <- lograte
  printout(head(rBoxCox))
  
  ## Generate the design array
  des.array = gen.comb.mtx(lev,ord="lex")+1
  printout(des.array)
  
  n.exp = dim(des.array)[1]
  exp.vec = 1:n.exp
  experiments = paste("Experiment",exp.vec)
  
  ## Read the interaction array
  
  inter.array = read.csv(inter.array.file,header=FALSE)
  
  ## Generate name.map
  name.map = lapply(lev,function(x) paste("L",1:x,sep=""))
  names(name.map) = fnames
  
  #kpi = rBoxCox[,1]
  
  
  ## Generate dat.list
  
  dat.exp  <-  function(exp.ind){
    out = NULL
    
    
    for(i in 1:nr){                          
      if(identical(des.array[exp.ind,],as.numeric(as.matrix(dat.file[i,2:nf])))) {
        #out = c(out,dat[i,kpi])
        out = c(out,dat.file[i,1])
      }      
    }
    
    return(out)
  }
  
  dat.list = lapply(exp.vec,dat.exp) # list of the kpi by experiment 
  
  names(dat.list) = experiments
  printout(dat.list)
  
  out = computeNonBinomialRobust(name.map,design.array=des.array,dat.list,interact.array=inter.array,x.mtx.main=NULL,
                                 comborder="lex",alpha=.05,thresh.outlier=5,by.pvalue=FALSE,err.threshold=1/5)
  
  return(out)
  
}




call.compute.NonBinomial.Quail <- function()
  ## Testing the LDL Cholesterol of Quail data, example 4.2.1 in Hettmansperger and McKean, 2011.
  ## 3/26/2012
{
  printout("Testing the LDL Cholesterol of Quail data, example 4.2.1 in Hettmansperger and McKean, 2011.")
  
  #     indir = file.path("X:/projects","Quail")
  #     file.csv = file.path(indir,"quail.csv")
  #     inter.array.file = file.path(indir,"interaction_array.csv")
  #     cnames.file=scan(file = file.csv, what = "character", nmax = -1, sep = ",",skip = 0, nlines = 1,quiet = TRUE) #read column names from the file
  # 
  #     fnames = head(cnames.file,-1)
  #     nf = length(fnames) # number of factors
  # 
  #     kpi = tail(cnames.file,1) # name of the kpi (last column)
  #     #kpi = dat[,2]
  # 
  #     col.classes = c(rep("factor",nf),"numeric") # column classes
  # 
  
  
  ## Read the data 
  #    dat.file = read.csv(file.csv,header=TRUE,colClasses=col.classes)
  
  data(quail)
  dat = quail
  fnames = colnames(dat)
  N = nrow(dat) # total number of observations
  
  lev = unique(dat$treat) # treatment levels
  n.treat = length(lev) # number of treatments
  n.i.treat = sapply(lev,function(x) sum(is.element(dat$treat,x))) # number of replications for each treatment
  kpi = tail(fnames,1) # name of the kpi (last column)
  
  
  # 
  # nr =dim(dat.file)[1] # number of rows
  # nc = dim(dat.file)[2] # number of columns
  # #lev = as.numeric(apply(dat.file[,1:nf],2,max))
  # lev=length(unique(dat.file$treat))
  # 
  # printout(head(dat.file))
  # 
  
  
  #   n.treat = length(lev) # number of treatments
  #    n.i.treat = sapply(lev,function(x) sum(is.element(dat$treat,x))) # number of replications for each treatment
  
  ## Contrasts
  #printout(contrasts(dat$treat))
  #contrasts(dat$treat) <- contr.sum
  #printout(contrasts(dat$treat))
  
  
  #kpi = tail(fnames,1) # name of the kpi (last column)
  
  
  #col.classes = c(rep("factor",nf),"numeric") # column classes
  
  form = ldl ~ treat
  #boxplot(form,dat) # boxplot, see Fig. 4.2.1, p. 296
  
  ## Rfit
  fit <- rfit.default(form,dat,intercept=FALSE) # R fit
  
  #if(dev.cur() >= 1) dev.new()
  #qqnorm(rstudent(rfit(ldl~treat,quail))) # studentized residual normal q-q plot, see Fig. 4.2.1, p. 296
  #    
  #    mf <- model.frame(form,dat)
  #    #printout(mt)
  #    mt <- attr(mf, "terms")
  #    attr(mt, "intercept") <- 0
  #    x <- model.matrix(mt, data = mf)
  #    printout(x)
  #    
  #   printout("Coefficients")
  #   printout(fit$coefficients)
  #printout("tauhat")
  #printout(fit$tauhat)
  #printout("taushat")
  #printout(fit$taushat)
  #printout("betahat")
  #printout(fit$betahat)
  #printout("fitted.values")
  #printout(fit$fitted.values)
  #   printout("median of the residuals")
  #   printout(median(fit$residuals))
  #    
  #  printout("Standard Error")
  se = sqrt(n.i.treat^(-1)*fit$tauhat^2+(fit$taushat^2-fit$tauhat^2)/N)
  #  printout(se)
  #    
  fit.df = data.frame(fit$coefficients,se)
  colnames(fit.df) = c("Coeff","Standard Error")
  rownames(fit.df) = names(fit$coefficients)
  printout(fit.df)
  
  printout("One way anova")
  fit.one <- oneway.rfit(dat$ldl,dat$treat)
  printout(fit.one$est)
  printout(fit.one$se)
  printout(fit.one$p.value)
  printout(summary(fit.one))
  
  #    ## ww package
  #     printout("Using the ww package")
  #    res.fit = wwfit(x, y=dat[,kpi], bij="WIL", center=F) 
  #printout(res.fit$tmp1$coefficients)
  #    
  #    
  #    ## LS fit
  #    printout("LS fit")
  #    form.ls = ldl ~ treat - 1
  #    fit.ls = lm(form.ls,dat)
  #    printout(summary(fit.ls))
  #    
  #    
  #    
  #    ## Rank-Based Tests of H0:mu1=...=muk, see p. 296 in Hettmansperger and McKean, 2011
  #    printout("Drop dispersion test")
  #    test.drop = drop.test(fit) # Rfit package
  #    test.drop.out = data.frame(test.drop$F,test.drop$tauhat,test.drop$df1,test.drop$df2,test.drop$p.value)
  #    colnames(test.drop.out) = c("Test Statistics","Scale","df1","df2","p-value")
  #    printout(test.drop.out)
  #    #printout(test.drop$F)
  #    #printout(test.drop$tauhat)
  #    #printout(c(test.drop$df1,test.drop$df2))
  #    #printout(test.drop$p.value)
  
  printout("computeNonBinomialRobust")
  # ## Generate the design array
  des.array = gen.comb.mtx(length(lev),ord="lex")+1
  printout("des.array")
  printout(des.array)
  
  n.exp = nrow(des.array)
  exp.vec = 1:n.exp
  experiments = paste("Experiment",exp.vec)
  nr = N
  nf = 1
  
  # 
  # ## Read the interaction array
  # 
  # #inter.array = read.csv(inter.array.file,header=FALSE)
  inter.array =NULL
  # 
  # ## Generate name.map
  name.map = list(treat  = paste("L",1:n.exp,sep=""))
  # name.map = lapply(lev,function(x) paste("L",1:x,sep=""))
  # names(name.map) = fnames
  
  ## Generate dat.list
  
  dat.exp  <-  function(exp.ind){
    out = NULL
    
    
    for(i in 1:nr){                          
      # if(identical(des.array[exp.ind,],as.numeric(as.matrix(dat.file[i,1:nf])))) {
      if(identical(as.numeric(des.array[exp.ind,]),as.numeric(as.matrix(dat[i,1:nf])))) {
        out = c(out,dat[i,kpi])
        #out = c(out,dat$kpi[i])
      }      
    }
    
    return(out)
  }
  
  dat.list = lapply(exp.vec,dat.exp) # list of the kpi by experiment 
  
  names(dat.list) = experiments
  printout(dat.list)
  
  out = computeNonBinomialRobust(name.map,design.array=des.array,dat.list,interact.array=inter.array,x.mtx.main=NULL,
                                 comborder="lex",alpha=.05,thresh.outlier=5,by.pvalue=FALSE,err.threshold=1/5)
  
  return(out)
  
}




call.compute.NonBinomial.Motors <- function()
## Testing the motor lifetime data, example 4.4.1 in Hettmansperger and McKean, 1998.
## 3/10/2011
{
printout("Testing the motor lifetime data, example 4.4.1 in Hettmansperger and McKean, 1998.")
    #indir = file.path("X:/projects/R/data","motors")
    indir = file.path("X:/projects","motors")
    file.csv = file.path(indir,"motors.csv")
    inter.array.file = file.path(indir,"interaction_array.csv")
    cnames.file=scan(file = file.csv, what = "character", nmax = -1, sep = ",",skip = 0, nlines = 1,quiet = TRUE) #read column names from the file

 fnames = head(cnames.file,-1)
 nf = length(fnames) # number of factors
 
 kpi = tail(cnames.file,1) # name of the kpi (last column)
 
 col.classes = c(rep("factor",nf),"numeric") # column classes
 
 ## Read the data 
 dat.file = read.csv(file.csv,header=TRUE,colClasses=col.classes)
#printout(head(dat.file))
#printout(dat.file[1,1:nf]) 
 nr =dim(dat.file)[1] # number of rows
 nc = dim(dat.file)[2] # number of columns
 lev = as.numeric(apply(dat.file[,1:nf],2,max))+1
 
 lograte = log(dat.file[,nc]) # log of the lifetime of motors
 dat.file[,nc] <- lograte
printout(head(dat.file))
 
 ## Generate the design array
 des.array = gen.comb.mtx(lev,ord="lex")
printout(des.array)
 
 n.exp = dim(des.array)[1]
 exp.vec = 1:n.exp
 experiments = paste("Experiment",exp.vec)
 
 ## Read the interaction array
 
 inter.array = read.csv(inter.array.file,header=FALSE)
#inter.array = NULL

 ## Generate name.map
 name.map = lapply(lev,function(x) paste("L",1:x,sep=""))
 names(name.map) = fnames
#printout(name.map)
 
 ## Generate dat.list
 
 dat.exp  <-  function(exp.ind){
                                out = NULL
                                
                               for(i in 1:nr){                          
                                  if(identical(des.array[exp.ind,],as.numeric(as.matrix(dat.file[i,1:nf])))) {
                                      out = c(out,dat.file[i,kpi])
                                  }      
                               }
                               return(out)
                               }
 
 dat.list = lapply(exp.vec,dat.exp) # list of the kpi by experiment 
 names(dat.list) = experiments
printout(dat.list)

    out = computeNonBinomialRobust(name.map,design.array=des.array,dat.list,interact.array=inter.array,x.mtx.main=NULL,
                             comborder="lex",alpha=.05,thresh.outlier=5,by.pvalue=FALSE,err.threshold=1/5)
    ## For debugging                         
#    out = computeNonBinomial(name.map,design.array=des.array,dat.list,interact.array=inter.array,x.mtx.main=NULL,
#                             comborder="lex",alpha=.05,thresh.outlier=5,by.pvalue=FALSE,err.threshold=1/5)
return(out) 
    
}


call.compute.NonBinomial.Marketing <- function()
## Testing the marketing data, example 4.6.1 in Hettmansperger and McKean, 1998.

{
 printout("Testing the marketing data, example 4.6.1 in Hettmansperger and McKean, 1998.")
 
 #indir = file.path("X:\\projects\\R","data","marketing")
 indir = file.path("X:\\projects","marketing")
 
 file.csv=file.path(indir,"marketing.csv") # data file
 
 inter.array.file = file.path(indir,"interaction_array.csv")
 
 
 
 cnames.file=scan(file = file.csv, what = "character", nmax = -1, sep = ",",skip = 0, nlines = 1,quiet = TRUE) #read column names from the file

 fnames = head(cnames.file,-1)
 nf = length(fnames) # number of factors
 
 kpi = tail(cnames.file,1) # name of the kpi (last column)
 
 col.classes = c(rep("factor",nf),"numeric") # column classes
 
 ## Read the data 
 dat.file = read.csv(file.csv,header=TRUE,colClasses=col.classes)
 printout(head(dat.file))
#printout(contrasts(dat.file$Fee)) 
 
 printout(dat.file[1,1:nf]) 
 nr =dim(dat.file)[1] # number of rows
 nc = dim(dat.file)[2] # number of columns
 lev = as.numeric(apply(dat.file[,1:nf],2,max))+1

 
 ## Generate the design array
 des.array = gen.comb.mtx(lev,ord="lex")
printout(des.array)
 
 n.exp = dim(des.array)[1]
 exp.vec = 1:n.exp
 experiments = paste("Experiment",exp.vec)
 
 ## Read the interaction array
 
 inter.array = read.csv(inter.array.file,header=FALSE)
 
 #inter.array = NULL
 
 ## Generate name.map
 name.map = lapply(lev,function(x) paste("L",1:x,sep=""))
 names(name.map) = fnames
printout(name.map)
 
 ## Generate dat.list
 
 dat.exp  <-  function(exp.ind){
                                out = NULL
                                
                               for(i in 1:nr){                          
                                  if(identical(des.array[exp.ind,],as.numeric(as.matrix(dat.file[i,1:nf])))) {
                                      out = c(out,dat.file[i,kpi])
                                  }      
                               }
                               return(out)
                               }
 
 dat.list = lapply(exp.vec,dat.exp) # list of the kpi by experiment 
 names(dat.list) = experiments
printout(dat.list)
 
 #rfit$tauhat

    out = computeNonBinomialRobust(name.map,design.array=des.array,dat.list,interact.array=inter.array,x.mtx.main=NULL,
                             comborder="lex",alpha=.05,thresh.outlier=5,by.pvalue=FALSE,err.threshold=1/5)
 


return(out) 
}




call.compute.NonBinomial.Eastcoast.two.methods <- function()
  ## Testing the east coast data given by EMEA.
  
{
  printout("Testing the east coast data")
  
  #indir = file.path("X:\\projects\\R","data","marketing")
  #indir = file.path("X:\\projects","EastCoast")
  indir = file.path("Y:\\widemileData\\EastCoast")
  
  ### generate the working matrix (datamatrix.scv), by combining the data matrix and the design matrix
  
  #data.file.csv=file.path(indir,"datamatrix.csv") # data file
  #design.file.csv=file.path(indir,"e16 (2 2 2 2 2 4 4).csv")
  
  ### read the raw dataset
  #results_61212 <- read.csv("X:/projects/EastCoast/results_61212.csv")
  #fn <- file.path(indir,"cc4290f1-e82f-4e22-8620-340990c36493.csv")
  fn = file.path(indir,"East Coast test_0","cc4290f1-e82f-4e22-8620-340990c36493.csv")
  results_61212 <- read.csv(fn)

  ### extract the tests with "step7_custom_data"
  #Revenue.new=results_61212$TotalRevenue[which(results_61212$conversionPoint == "step7_custom_data")]
  Revenue.new=results_61212$TotalRevenue[which(results_61212$conversionEvent == "step7_custom_data")]
  head(Revenue.new)
  
  ### convert from factor to numeric, 
  #Revenue.new = as.numeric(levels(Revenue.new))[Revenue.new] 
  #head(Revenue.new)
  
  Revenue.new = as.numeric(Revenue.new)
  head(Revenue.new)
  
  ### extract the experiments
  #experimentName.new=results_61212$experimentName[which(results_61212$conversionPoint == "step7_custom_data")]
  experimentName.new=results_61212$experimentName[which(results_61212$conversionEvent == "step7_custom_data")]
  head(experimentName.new)
  
  experiment = substr(experimentName.new,12,13)
  experiment = as.numeric(experiment)
  
  ### initialize a working matrix
  n=length(experimentName.new)
  #working.data = matrix(0,nrow=n,ncol=9)
  working.data = matrix(0,nrow=n,ncol=5)
  
  ### combine the experiment and revenue into the working dataset
  #working.data[,1] = experiment
  #working.data[,9] = Revenue.new
  
  working.data[,1] = experiment
  working.data[,5] = Revenue.new
  
  ### read the design matrix
  #e16 <- read.csv("X:/projects/EastCoast/e16 (2 2 2 2 2 4 4).csv", header=F) 
  #e16 <- matrix(unlist(e16),nrow = 16,byrow=F)
  
  #e16 <- read.csv("X:/projects/EastCoast/e16 (2 2 4).csv", header=F) 
  #e16 <- matrix(unlist(e16),nrow = 16,byrow=F)
  
  e16 <- gen.comb.mtx(lev=c(2,2,4),ord="standard")
  
  ### assign the design matrix row into the data matrix based on the experiment number
  for(i in 1:n)
  {
      #working.data[i,2:8] = e16[working.data[i,1],] 
      working.data[i,2:4] = e16[working.data[i,1],]
  }
  head(working.data)
  
  working.data = working.data[1:1000,]
   
  ### write the data into a csv file
  write.csv(working.data, file.path(indir,"working.data.csv"), row.names=FALSE)
  
  ### sample it, to get a representative sampling, including all the 16 experiments
  #working.data = working.data[1:1000,]
  
  ### in the working data, delete the first column of experiments
  working.data1 = working.data[,-1]
  head(working.data1)
  
  ### write the data into a csv file
  write.csv(working.data1, file.path(indir,"working.data1.csv"), row.names=FALSE)
  
  ### read the generated data matrix
  file.csv=file.path(indir,"working.data1.csv")
  
  inter.array.file = file.path(indir,"interaction_array.csv")

  cnames.file=scan(file = file.csv, what = "character", nmax = -1, sep = ",",skip = 0, nlines = 1,quiet = TRUE) #read column names from the file
  
  #fnames = head(cnames.file,-c(1,9))
  fnames = head(cnames.file,-1)
  #fnames = tail(fnames,-1)
  #fnames = head(cnames.file)
  
  nf = length(fnames) # number of factors
  
  kpi = tail(cnames.file,1) # name of the kpi (last column)
  
  col.classes = c(rep("factor",nf),"numeric") # column classes
  
  ## Read the data 
  dat.file = read.csv(file.csv,header=TRUE,colClasses=col.classes)
  class(dat.file)
  printout(head(dat.file))
  #printout(data.class(dat.file[,8]))
  printout(data.class(dat.file[,4]))
  
  ### convert from data.frame into matrix
  ### notice: need to convert from factor to numeric
  
  
  #dat.file = data.matrix(dat.file)
  #class(dat.file) 
  
  printout(head(dat.file))
  printout(dat.file[1,1:(nf+1)]) 
  nr =dim(dat.file)[1] # number of rows
  nc = dim(dat.file)[2] # number of columns
  
  lev = as.numeric(apply(dat.file[,1:nf],2,max))+1
  
  
  ## Generate the design array
  #des.array = gen.comb.mtx(lev,ord="lex")
  des.array = e16
  printout(des.array)
  
  
  n.exp = dim(des.array)[1]
  exp.vec = 1:n.exp
  experiments = paste("Experiment",exp.vec)
  
  ## Read the interaction array
  
  #inter.array = read.csv(inter.array.file,header=FALSE)
  
  inter.array = NULL
  
  ## Generate name.map
  name.map = lapply(lev,function(x) paste("L",1:x,sep=""))
  names(name.map) = fnames
  printout(name.map)
  
  ## Generate dat.list
  
  dat.exp  <-  function(exp.ind){
    out = NULL
    
    for(i in 1:nr){                          
      if(identical(as.numeric(des.array[exp.ind,]),as.numeric(as.matrix(dat.file[i,1:nf])))) {
        out = c(out,dat.file[i,kpi])
      }      
    }
    return(out)
  }

  
  dat.list = lapply(exp.vec,dat.exp) # list of the kpi by experiment 
  names(dat.list) = experiments
  printout(dat.list)
  
  
  out = computeNonBinomialRobust(name.map,design.array=des.array,dat.list,interact.array=inter.array,x.mtx.main=NULL,
                                 comborder="lex",alpha=.05,thresh.outlier=5,by.pvalue=FALSE,err.threshold=1/5)
  
  return(out) 
}




#################################################################
##### Purpose: check whether our function get correct results

# n = 1000
# mu = 2
# lev = c(2, 2, 3)
# 
# # leveff = c(1, 2, 1.5, 0.6, 0.8, 0.7, 0.1)
# leveff[[1]] = c(-1, 1)
# leveff[[2]] = c(-0.5, 0.5)
# leveff[[3]] = c(1, -3, 2)

self.check.leveleffect<-function(n, mu, lev, m, inter.array, des.array, leveff)
{ 
  
  ## Input: 
  ### n: the number of observations in our simulation
  ### mu: the overall mean 
  ### lev: the factor level combination vector
  ### m: number of effects we are interested in
  ### leveff: level effects we assigned
  
  n = 1000
  mu = 2
  lev = c(2, 2, 3)
  
  # intercation array, and it determines which effects we are interested in
  # inter.array = NULL
  # inter.array = matrix(ncol=length(lev),byrow=TRUE,c(1,1,0,1,0,1,0,1,1,1,1,1))
  # inter.array = matrix(ncol=length(lev),byrow=TRUE,c(1,1,0,1,0,1))
  inter.array = matrix(ncol=length(lev),byrow=TRUE,c(1,1,0))
  printout(inter.array)
  
  # if only consider the main effect
  if(inter.array == NULL)
  {
    m = length(lev)   
    
    #leveff = list()
    #length(leveff) <- m
    leveff = vector("list", m)
    
    for(i in 1:m)
    {
      leveff[[i]] = rep(0,lev[i])
    }
    
    # leveff = c(1, 2, 1.5, 0.6, 0.8, 0.7, 0.1)
    leveff[[1]] = c(-1, 1)
    leveff[[2]] = c(-0.5, 0.5)
    leveff[[3]] = c(1, -3, 2)
  }
  
  # based on lev and intercation array, get the mapping matrix
  map.mtx = matrix("list", dim(inter.array)[1])
    
    
    
  # based on the design array row and mapping matrix row, derive how many level effects are included,and the level effects list
    
  # based on lev and mapping matrix row, find how many interaction effects
  n.interaction = 0

  
  # simulate the random terms 
  epsilon = rnorm(n) 
  
  # kpi
  y = rep(0,n)
  
  # experiment
  exp.col = rep(0,n)
  
  # design array
  des.array = gen.comb.mtx(lev,ord="lex")
  printout(des.array)
  
  
  
  # how many experiments
  exp.number = dim(des.array)[1]
  #exp = c(1:exp.number)
  
  # among all the observations, how many times each experiment appears on average
  # n1 = floor(n/exp.number)
  
  # assign experiments to each observation randomly
  exp.col = sample(1:exp.number,n,replace=T)
  
  
  n.exp = dim(des.array)[1]
  exp.vec = 1:n.exp
  experiments = paste("Experiment",exp.vec)
  
  
  ## Generate name.map
  name.map = lapply(lev,function(x) paste("L",1:x,sep=""))
  names(name.map) = fnames
  printout(name.map)

  # based on the level, we choose the level effects from leveff list
  pick.leveff <- function(level,leveff)
  {
    out = NULL
    leveff.certain = rep(0,m)
    
    for(i in 1:m)
     {                          
        leveff.certain[i] = leveff[[i]][level[i]+1]     
     }
    #printout(leveff.certain)
    
    # calculate the sum effect of those level effects
    out = sum(leveff.certain)

    return(out)
  }
  
  
  # simulate kpi based on experiment number
  for(i in 1:n)
  {
    # the row in the design array corresponding to the experiment number
    level = des.array[exp.col[i],]
    #printout(level)
     
    a = pick.leveff(level,leveff)
    #printout(a)
   
    y[i] = mu + a + epsilon[i]
    
  }
  
  #####################
  printout("Testing the simulation data.")
  
  #indir = file.path("X:\\projects\\R","data","marketing")
  indir = file.path("X:\\projects","simulation")
    
  working.data = matrix(0,nrow=n,ncol=5)
  
  ### combine the experiment and revenue into the working dataset
  
  working.data[,1] = exp.col
  working.data[,5] = y
  
  ### assign the design matrix row into the data matrix based on the experiment number
  for(i in 1:n)
  {
    working.data[i,2:4] = des.array[working.data[i,1],]
  }
  head(working.data)
  
  
  ### write the data into a csv file
  write.csv(working.data, file.path(indir,"working.data.csv"), row.names=FALSE)
  
  ### in the working data, delete the first column of experiments
  working.data1 = working.data[,-1]
  head(working.data1)
  
  ### write the data into a csv file
  write.csv(working.data1, file.path(indir,"working.data1.csv"), row.names=FALSE)
  
  ### read the generated data matrix
  file.csv=file.path(indir,"working.data1.csv")
  
  #inter.array.file = file.path(indir,"interaction_array.csv")
  inter.array.file = NULL
  
  cnames.file=scan(file = file.csv, what = "character", nmax = -1, sep = ",",skip = 0, nlines = 1,quiet = TRUE) #read column names from the file
  
  fnames = head(cnames.file,-1)
  
  nf = length(fnames) # number of factors
  
  kpi = tail(cnames.file,1) # name of the kpi (last column)
  
  col.classes = c(rep("factor",nf),"numeric") # column classes
  
  ### read data file
  dat.file = read.csv(file.csv,header=TRUE,colClasses=col.classes)
  head(dat.file)
  
  
  ## Generate dat.list
  
  dat.exp  <-  function(exp.ind){
    out = NULL
    
    for(i in 1:n){                          
      if(identical(as.numeric(des.array[exp.ind,]),as.numeric(as.matrix(dat.file[i,1:nf])))) {
        out = c(out,dat.file[i,kpi])
      }      
    }
    return(out)
  }
  
  dat.list = lapply(exp.vec,dat.exp) # list of the kpi by experiment 
  names(dat.list) = experiments
  printout(dat.list)
  
  out = computeNonBinomialRobust(name.map,design.array=des.array,dat.list,interact.array=inter.array,x.mtx.main=NULL,
                                 comborder="lex",alpha=.05,thresh.outlier=5,by.pvalue=FALSE,err.threshold=1/5)
  
  return(out)  
  
}



# indir = file.path("X:\\projects","Barclaycard")
# design.file.csv=file.path(indir,"e12 (2 2 2 3).csv")
# design.array=read.csv(design.file.csv,header=FALSE)
# lev = c(2,2,2,3)
# fnames=c("V1","V2","V3","V4")
# name.map = lapply(lev,function(x) paste("L",1:x,sep=""))
# names(name.map) = fnames
# printout(name.map)



call.compute.NonBinomial.Hamilton <- function()
  ## Testing the one sample location problem Hollander&Wolfe(1999), p.51-57 
  ## Example 3.1 and 3.3 Hamilton Depression Scale Factor
{
  
  printout("Testing the one sample location problem")
#  source(file.path("X:\\projects\\R\\data\\hamilton\\hamilton.txt"))
#  Z = sort(hamilton$Y-hamilton$X)
   Z=c(-1.022, -0.952, -0.620, -0.590, -0.490, -0.430, -0.010,  0.080,  0.147)
  nobs = length(Z)
  printout(Z)  
  ## Testing the repetitions
  #Zrep = sort(c(Z,rep(Z[3],2),rep(Z[5],3))) # for testing the repetitions
  #Z = Zrep
  #nobs = length(Z)
  #####################################
  
#   ### Testing the abs differences
#   indij = t(combn(nobs,2))
#   diff.arr = get.abdord(vec.z=Z,indij=indij)
#   
#   diff.arr.tbl = matrix(ncol=nobs,nrow=nobs)
#   counter = 1
#   for(i in 1:nrow(indij))
#   {
#     diff.arr.tbl[indij[i,1], indij[i,2]] = diff.arr[counter]
#     counter = counter + 1
#   }
#   
# printout(diff.arr.tbl)
  
  ######################################################################

  ztbl = t(sapply(1:(nobs),function(x) c(rep(NA,(x-1)),(Z[(x):nobs]+Z[x]))))
printout(ztbl)
#return(ztbl)
#ind =56
printout(paste("index",ind))
  out = sort.pair(Z,ind,descreasing=FALSE)
  
  
  #out = compute.HLestimates(Z,alpha=0.04,fast=F)

}


# gen.bounds<-function(n,m)
# ## Generates bounds for an nxm rectangular array
# {
#   gridij = expand.grid(i=1:n,j=1:m)
#   lb = apply(gridij,1,function(x) x[1]*x[2]) 
#   lb.mtx = matrix(lb,nrow=n,ncol=m)
# printout(lb.mtx)
#   ub = apply(gridij,1,function(x) n*m+1-(n-x[1]+1)*(m-x[2]+1)) 
#   ub.mtx = matrix(ub,nrow=n,ncol=m)
# printout(ub.mtx)
# }
  
# test.sort<-function(sz)
# ## Tests timing of sorting algorithms
# {
#   #nrep = 1000
#   #x = rep(rnorm(floor(sz/nrep)),nrep)
#   x = rnorm(sz)
# printout(length(x)) 
#   #sort(x, decreasing = FALSE, na.last = NA,  method="quick")
#   printout(system.time(sort(x, method = "quick")))
#   
#   printout(system.time(sort(x, method = "shell")))
# }
# 
# 
# compare.rows <- function(x,mtx)
# {
#   apply(x,1,function(i){ any(apply(mtx,1,identical,i)) }) 
#   
# }


callComputeAnovaTableCompareRobust <- function(
        name.file,
        design.array.file,
        conversions.report.file,
        interact.array.file=NULL,
        comborder="lex",
        alpha = .05)
{
    
    options(mathdebug=TRUE)
    inputdata <- getInputData(name.file,design.array.file,conversions.report.file)
    
    fnames <- names(inputdata$name.map)
    con.rate.v <- inputdata$conversions$Rate  
    name.map=inputdata$name.map
    design.array = inputdata$design.array
    
    #model.formula = "rates ~ ."
    ## Interaction array
    interact.array <- NULL
    if(!is.null(interact.array.file)) # No interact.array.file
    {
        interact.array <- read.table(interact.array.file,header = FALSE,sep = ",",dec=".",fill = TRUE, strip.white = TRUE,comment.char="",col.names=fnames)
        #formulastr=interactArray2Formula(colnames(design.array),as.matrix(interact.array)) 
        
        #model.formula = formulastr
    }
  
    cat("\n\n\n")
    
    #call old version of computeAnovaTable
    printout("computeAnovaTable")
printout(name.map)
    startime=proc.time()
    retvalues = computeAnovaTable(name.map,design.array,con.rate.v,interact.array,x.mtx.main=NULL,comborder,alpha)
    
    endtime = proc.time()
    printout("elapsed time old method ")
    printout(endtime-startime)
    
    
    ### call  new version of anova calculations

    design.array.factors = designArrayAsFactors(name.map,design.array)
    startime=proc.time()
    #retvaluesAOV = computeAOV(model.formula,design.array.factors,con.rate.v)
    #retvaluesAOV = computeAOV(inputdata$name.map,interact.array,inputdata$design.array,inputdata$conversions$Rate)#return(list(sig.ind=sig.ind, aov=aovout,aov.revised=revisedAOV))
    
    experiments = rownames(inputdata$design.array)
    #dat.list = lapply(experiments,function(x) dat.file[dat.file[,exp.name]==x,kpi]) # list of the kpi by experiment
    dat.list = lapply(con.rate.v,identity)
    names(dat.list) = rownames(inputdata$design.array)
    retvaluesAOV <- computeNonBinomial(inputdata$name.map,inputdata$design.array,dat.list,interact.array=interact.array,x.mtx.main=NULL,
                             comborder="lex",alpha=.05,thresh.outlier=50,by.pvalue=FALSE,err.threshold=1/5)
    endtime=proc.time()
    printout("elapsed time new ")
    printout(endtime-startime)
    printout("-----------------------------------------------")
    #compare results just for level influences and intial anova table
    
    #for(i in 1:(length(retvaluesAOV$influences$tables))){
    for(i in 1:(length(retvalues$levelEffects))){
        cat("\nold: ")
        cat(fnames[i])
        cat(": ")
        print(retvalues$levelEffects[[i]])
        
        cat("\nnew: ")
        cat(fnames[i])
        cat(": ")
        #print(retvaluesAOV$influences$tables[[i]])
        print(retvaluesAOV$levelEffects[[i]])
        cat("\n")
    }
    cat("\n----anova old\n")
    printout(retvalues$initialANOVATable$anovatable[,1:3])
    
    cat("\n-----anova new\n")
    #printout(summary(retvaluesAOV$anova))
    printout(retvaluesAOV$initialANOVATable$anovatable[,1:3])
    
    cat("\n-----old sig ind ")
    cat(retvalues$initialANOVATable$sig.ind)
    
    cat("\n---new sig ind ")
    #cat(retvaluesAOV$sig.ind)
    cat(retvaluesAOV$initialANOVATable$sig.ind)
    
    cat("\n--- old revised anova\n")
    #don't print out factor influence or significance columns for easier comparison
    printout(retvalues$revisedANOVATable$anovatable[c(1,2,3,4,6)])
    
    cat("\n--- new revised anova\n")
    printout(retvaluesAOV$revisedANOVATable$anovatable[c(1,2,3,4,6)])
 #   
#    if(!is.null(retvaluesAOV$revisedAOV)){
#        cat("\n--- new revised anova\n")
#        printout(anova(retvaluesAOV$revisedAOV))
#    }
    return(list(retvalues = retvalues,
                    retvaluesAOV=retvaluesAOV,
                    name.map=name.map,
                    inputdata=design.array.factors))
    
} # end of getANOVATable
