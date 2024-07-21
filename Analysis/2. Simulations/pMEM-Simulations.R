##
## Development of predictive Moran's eigenvector maps (pMEM)
## Monte-Carlo simulations analysis script
## 
## rm(list=ls())

library(magrittr)
library(knitr)
library(RPostgreSQL)
library(dplyr)
library(pMEM)
library(yaml)

source("../pMEM-aux.R")

## Obtain the simulation data from the database
if(FALSE) {
  
  ## Confidentials:
  creds <- yaml.load_file("../../Credentials/Credentials.yml")
  
  conn <- sshPortRedirect(creds$hostname, creds$username, 5432, 5434)
  
  dbConnect(
    RPostgreSQL::PostgreSQL(),
    host="localhost",
    port="5434",
    user=creds$username,
    password=creds$password,
    dbname=creds$dbname
  ) -> sim
  
  ## showProcessingJobs(sim)
  
  ## All the simulation results:
  showCompletedJobs(sim) %>%
    as_tibble -> jobs
  
  dbDisconnect(sim)
  sshPortShutdown(conn)
  
  rm(creds,conn,sim)
  
  jobs %>%
    select(!worker) %>%
    mutate(qof = log10(msd) - log10(mse)) %>%
    select(!msd & !mse) -> simRes
  
  rm(jobs)
  
  simRes$map %<>% as.factor
  simRes$sample %<>% as.factor
  simRes$fun %<>% as.factor
  
  saveRDS(simRes, "../../Data/simRes.rds")
  
  ## Simulations with the best DWF
  best <- c()
  
  for(i in unique(simRes$map))
    for(j in unique(simRes$sample))
      for(k in unique(simRes$n))
        which(
          (simRes$map==i)&
            (simRes$sample==j)&
            (simRes$n==k)
        ) %>%
    {.[which.max(simRes$qof[.])]} %>%
    c(best,.) -> best
  
  simResBest <- simRes[best,]
  
  rm(best,i,j,k)
  
  saveRDS(simResBest, "../../Data/simResBest.rds")
  
} else {
  simRes <- readRDS("../../Data/simRes.rds")
  simResBest <- readRDS("../../Data/simResBest.rds")
}

## Analysis of the simulation results - step 1
if(FALSE) {
  
  simRes %>%
    aov(
      qof ~ (log10(n) + fun + map + sample)^3,
      data=.,
      contrast = list(
        fun = contr.treatment,
        map = contr.helmert,
        sample = contr.helmert
      )
    ) %>%
    anova -> aovAll
  
  simResBest %>%
    aov(
      qof ~ log10(n) * map * sample,
      data=.,
      contrast = list(
        map = contr.helmert,
        sample = contr.helmert
      )
    ) %>%
    anova -> aovBest
  
  aovtab <- list(all = aovAll)
  
  class(aovtab$all) <- class(aovtab$all)[2L]
  
  c(
    "$\\log_{10}n$","$DWF$","$Map$","$Sample$","$\\log_{10}n \\times DWF$",
    "$\\log_{10}n \\times Map$", "$\\log_{10}n \\times Sample$",
    "$DWF \\times Map$","$DWF \\times Sample$","$Map \\times Sample$",
    "$\\log_{10}n \\times DWF \\times Map$",
    "$\\log_{10}n \\times DWF \\times Sample$",
    "$\\log_{10}n \\times Map \\times Sample$",
    "$DWF \\times Map \\times Sample$","$\\mathrm{Residuals}$"
  ) -> rownames(aovtab$all)
  aovtab$all[["Sum Sq"]] <- NULL
  aovtab$all[["Mean Sq"]] <- NULL
  aovtab$all[["F value"]] %<>% signif(4) %>% as.character
  aovtab$all[["F value"]][is.na(aovtab$all[["F value"]])] <- ""
  
  aovtab$all[["Pr(>F)"]] %<>%
    symnum(
      corr = FALSE, na = FALSE, 
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("$< 0.000\\,1$", "$< 0.001$", "$< 0.01$", "$< 0.05$",
                  "$> 0.05$")
    )
  
  aovtab$all[["Pr(>F)"]][is.na(aovtab$all[["Pr(>F)"]])] <- ""
  
  kable(
    aovtab$all,
    col.names = c("$\\nu$","$F_{\\nu,\\nu_{res}}$","$P$")
  )
  
  aovtab$best <- aovBest
  
  class(aovtab$best) <- class(aovtab$best)[2L]
  
  c(
    "$\\log_{10}n$", "$Map$", "$Sample$",
    "$\\log_{10}n \\times Map$", "$\\log_{10}n \\times Sample$",
    "$Map \\times Sample$", "$\\log_{10}n \\times Map \\times Sample$",
    "$\\mathrm{Residuals}$"
  ) -> rownames(aovtab$best)
  
  aovtab$best[["Sum Sq"]] <- NULL
  aovtab$best[["Mean Sq"]] <- NULL
  aovtab$best[["F value"]] %<>% signif(4) %>% as.character
  aovtab$best[["F value"]][is.na(aovtab$best[["F value"]])] <- ""
  
  aovtab$best[["Pr(>F)"]] %<>%
    symnum(
      corr = FALSE, na = FALSE, 
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("$< 0.000\\,1$", "$< 0.001$", "$< 0.01$", "$< 0.05$",
                  "$> 0.05$")
    )
  
  aovtab$best[["Pr(>F)"]][is.na(aovtab$best[["Pr(>F)"]])] <- ""
  
  kable(
    aovtab$best,
    col.names = c("$\\nu$","$F_{\\nu,\\nu_{res}}$","$P$")
  )
  
  save(aovAll, aovBest, aovtab, file="../../Data/sim_AOV_tab.rda")
  
} else load(file="../../Data/sim_AOV_tab.rda")

## Analysis of the simulation results - step 2
if(FALSE) {
  
  simResBest %>%
    lm(
      qof ~ log10(n) * map * sample,
      data=.,
      contrast = list(
        map = contr.helmert,
        sample = contr.helmert
      )
    ) -> lmSim
  
  cbind(1,log10(c(10,20,50,100,200,500))) %*%
    coefficients(lmSim)[1L:2L] %>%
    {1 - 10^-.} -> NPred
  ## -log10(1-NPred)
  
  dimnames(NPred) <- list(c(10,20,50,100,200,500),"Mean P")
  
  saveRDS(NPred, "../../Data/NPred.rds")
  
} else NPred <- readRDS("../../Data/NPred.rds")

## simResBest$fun %>% summary %>% sort(decreasing = TRUE)

## Analysis of the simulation results - step 3
if(FALSE) {
  
  mapPred <- c()
  
  for(i in levels(simResBest$map))
    simResBest %>%
    .[.$map==i,] %>%
    lm(
      qof ~ log10(n) * sample,
      data=.,
      contrast = list(
        sample = contr.helmert
      )
    ) %>%
    coefficients %>%
    .[1L:2L] %>%
    {cbind(1,log10(c(10,20,50,100,200,500))) %*% .} %>%
    {1 - 10^-.} %>%
    as.numeric %>%
    rbind(mapPred,.) -> mapPred
  
  list(
    levels(simResBest$map),
    c(10,20,50,100,200,500)
  ) -> dimnames(mapPred)
  
  ## for(i in 2L:ncol(mapPred))
  ##   mapPred[,i] <- mapPred[,i] - rowSums(mapPred[,1L:(i - 1L),drop=FALSE])
  
  saveRDS(mapPred, "../../Data/mapPred.rds")
  
} else mapPred <- readRDS("../../Data/mapPred.rds")

## Analysis of the simulation results - step 4
if(FALSE) {
  
  samplePred <- c()
  
  for(i in levels(simResBest$sample))
    simResBest %>%
    .[.$sample==i,] %>%
    lm(
      qof ~ log10(n) * map,
      data=.,
      contrast = list(
        map = contr.helmert
      )
    ) %>%
    coefficients %>%
    .[1L:2L] %>%
    {cbind(1,log10(c(10,20,50,100,200,500))) %*% .} %>%
    {1 - 10^-.} %>%
    as.numeric %>%
    rbind(samplePred,.) -> samplePred
  
  list(
    levels(simResBest$sample),
    c(10,20,50,100,200,500)
  ) -> dimnames(samplePred)
  
  ## for(i in 2L:ncol(samplePred))
  ##   samplePred[,i] <- samplePred[,i] - rowSums(samplePred[,1L:(i - 1L),drop=FALSE])
  
  saveRDS(samplePred, "../../Data/samplePred.rds")
  
} else samplePred <- readRDS("../../Data/samplePred.rds")

## Analysis of the simulation results - step 5
if(FALSE) {
  
  funPred <- c()
  
  for(i in levels(simRes$fun))
    simRes %>%
    .[.$fun==i,] %>%
    lm(
      qof ~ log10(n) * map * sample,
      data=.,
      contrast = list(
        map = contr.helmert,
        sample = contr.helmert
      )
    ) %>%
    coefficients %>%
    .[1L:2L] %>%
    {cbind(1,log10(c(10,20,50,100,200,500))) %*% .} %>%
    {1 - 10^-.} %>%
    as.numeric %>%
    rbind(funPred,.) -> funPred
  
  list(
    levels(simRes$fun),
    c(10,20,50,100,200,500)
  ) -> dimnames(funPred)
  
  ## for(i in 2L:ncol(funPred))
  ##   funPred[,i] <- funPred[,i] - rowSums(funPred[,1L:(i - 1L),drop=FALSE])
  ## funPred %<>% .[order(rowSums(.)),]
  
  saveRDS(funPred, "../../Data/funPred.rds")
  
} else funPred <- readRDS("../../Data/funPred.rds")

## Analysis of the simulation results - step 6
if(FALSE) {
  
  simResBest$fun %>%
    tapply(.,.,length) %>%
    sort(decreasing = TRUE) -> dwf
  
  saveRDS(dwf, "../../Data/bestDwf.rds")
  
} else dwf <- readRDS("../../Data/bestDwf.rds")

## Plotting the results: simulations per maps
if(FALSE) {
  
  png(file="../../Image/Simulations-map.png", width = 800, height = 400)
  par(mfrow=c(1L,1L), mar=c(3.0,4.6,1.0,4.6))
  
  mapPred %>% scale(scale=FALSE) %>% svd(nu = 1L) %>% .$u %>% order -> ord
  
  plot(
    NA, xlim=c(1,length(ord)), ylim=c(-0.05,1), axes=FALSE,
    ylab=expression(paste("Mean", ~ italic(P)^2))
  )
  axis(1L, 1:length(ord), labels=rownames(mapPred)[ord], las=2L)
  axis(2L, las=2L)
  col <- grey(seq(0.3, 0.9, length.out=ncol(mapPred)))
  for(i in 1L:ncol(mapPred)) {
    ## lines(x=1:length(ord), y=mapPred[ord,i])
    points(x=1:length(ord), y=mapPred[ord,i], pch=22L, bg=col[i], cex=2.5)
  }
  
  legend(
    x = length(ord) + 0.75, y = 0.7, pch=22L, xpd=TRUE, pt.cex=2, box.lwd = 0,
    legend = parse(text=sprintf("italic(n) == %s",c(500,200,100,50,20,10))),
    pt.bg = rev(col)
  )
  
  dev.off()
  
  rm(ord,i,col)
  
}

## Plotting the results: simulations per subsets
if(FALSE) {
  
  png(file="../../Image/Simulations-subset.png", width = 800, height = 400)
  par(mfrow=c(1L,1L), mar=c(3.0,4.6,2.0,4.6))
  
  samplePred %>% scale(scale=FALSE) %>% svd(nu = 1L) %>% .$u %>% order -> ord
  
  plot(
    NA, xlim=c(1,length(ord)), ylim=c(-0.05,1), axes=FALSE,
    ylab=expression(paste("Mean", ~ italic(P)^2))
  )
  axis(1L, 1:length(ord), labels=rownames(samplePred)[ord], las=2L)
  axis(2L, las=2L)
  col <- grey(seq(0.3, 0.9, length.out=ncol(samplePred)))
  for(i in 1L:ncol(samplePred)) {
    ## lines(x=1:length(ord), y=samplePred[ord,i])
    points(x=1:length(ord), y=samplePred[ord,i], pch=22L, bg=col[i], cex=2.5)
  }
  
  legend(
    x = length(ord) + 0.75, y = 0.7, pch=22L, xpd=TRUE, pt.cex=2, box.lwd = 0,
    legend = parse(text=sprintf("italic(n) == %s",c(500,200,100,50,20,10))),
    pt.bg = rev(col)
  )
  
  dev.off()
  
  rm(ord,i,col)
  
}

labSwap <- c(linear="Linear", power="Power", hyperbolic="Hyperbolic",
             spherical="Spherical", exponential="Exponential",
             Gaussian="Gaussian", hole_effect="Hole effect")

## Plotting the results: simulations per best distance-weighting function
if(FALSE) {
  
  png(file="../../Image/Simulations-function.png", width = 400, height = 400)
  par(mar=c(7,4.1,0.6,5))
  
  funPred %>% scale(scale=FALSE) %>% svd(nu = 1L) %>% .$u %>% order -> ord
  
  plot(
    NA, xlim=c(0.5,length(ord) + 0.5), ylim=c(-0.05,1), axes=FALSE,
    ylab=expression(paste("Mean", ~ italic(P)^2)), xlab=""
  )
  axis(1L, 1:length(ord), labels=labSwap[rownames(funPred)[ord]], las=2L)
  axis(2L, las=2L)
  col <- grey(seq(0.3, 0.9, length.out=ncol(funPred)))
  for(i in 1L:ncol(funPred)) {
    ## lines(x=1:length(ord), y=funPred[ord,i])
    points(x=1:length(ord), y=funPred[ord,i], pch=22L, bg=col[i], cex=2.5)
  }
  
  legend(
    x = length(ord) + 0.75, y = 0.7, pch=22L, xpd=TRUE, pt.cex=2, box.lwd = 0,
    legend = parse(text=sprintf("italic(n) == %s",c(500,200,100,50,20,10))),
    pt.bg = rev(col)
  )
  
  dev.off()
  
  rm(ord,i,col)
  
}
