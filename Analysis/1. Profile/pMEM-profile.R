##
## Development of predictive Moran's eigenvector maps (pMEM)
## Eigenfunction profile script
## 
## rm(list=ls())

library(magrittr)
library(pMEM)

source("../pMEM-aux.R")

## Distance-weighting functions figure
if(FALSE) {
  
  ## X11(width=4, height=12)
  png(file="../../Image/Common weighting functions.png", width = 400, height = 1200)
  par(mfrow=c(7,1),mar=c(5.1,5.1,0.6,0.6))
  
  d <- seq(0,5,0.001)
  
  plot(x = d, y = genDWF("linear",1)(d), type = "l", ylim = c(0,1), las = 1L,
       xlab = "", ylab = "", cex.axis = 2, lwd = 2)
  lines(x = d, y = genDWF("linear",2)(d), col = "red", lwd = 2)
  lines(x = d, y = genDWF("linear",3)(d), col = "blue", lwd = 2)
  text(x = 1.5, y = 0.85, label = "Linear", adj = 0,
       cex = 2)
  legend(x = 3.2, y = 0.85, lwd=1, col=c("black","red","blue"), cex=1.25,
         legend = expression(
           italic(d[max]) == 1,
           italic(d[max]) == 2,
           italic(d[max]) == 3))
  
  plot(x = d, y = genDWF("power",1,1)(d), type = "l", ylim = c(0,1),
       las = 1L, xlab = "", ylab = "", cex.axis = 2, lwd = 2)
  lines(x = d, y = genDWF("power",2,0.5)(d), col = "red", lwd = 2)
  lines(x = d, y = genDWF("power",3,0.5)(d), col = "blue", lwd = 2)
  text(x = 1.5, y = 0.85, label = "Power", adj = 0,
       cex = 2)
  legend(x = 3.2, y = 0.85, lwd=1, col=c("black","red","blue"), cex=1.25,
         legend = expression(
           paste(italic(d[max]) == 1, ", ",  alpha == 1),
           paste(italic(d[max]) == 2, ", ",  alpha == 0.5),
           paste(italic(d[max]) == 3, ", ",  alpha == 0.5)))
  
  plot(x = d, y = genDWF("hyperbolic",1,1)(d), type = "l", ylim = c(0,1),
       las = 1L, xlab = "", ylab = "", cex.axis = 2, lwd = 2)
  lines(x = d, y = genDWF("hyperbolic",2,0.5)(d), col = "red", lwd = 2)
  lines(x = d, y = genDWF("hyperbolic",3,0.5)(d), col = "blue", lwd = 2)
  text(x = 1.5, y = 0.85, cex = 2,
       label = "Hyperbolic", adj=0)
  legend(x = 3.2, y = 0.85, lwd=1, col=c("black","red","blue"), cex=1.25,
         legend = expression(
           paste(italic(d[max]) == 1, ", ",  alpha == 1),
           paste(italic(d[max]) == 2, ", ",  alpha == 0.5),
           paste(italic(d[max]) == 3, ", ",  alpha == 0.5)))
  
  plot(x = d, y = genDWF("spherical",1)(d), type = "l", ylim = c(0,1), las = 1L,
       xlab = "", ylab = "", cex.axis = 2, lwd = 2)
  lines(x = d, y = genDWF("spherical",2)(d), col = "red", lwd = 2)
  lines(x = d, y = genDWF("spherical",3)(d), col = "blue", lwd = 2)
  text(x = 1.5, y = 0.85, label = "Spherical", adj = 0, cex = 2)
  legend(x = 3.2, y = 0.85, lwd=1, col=c("black","red","blue"), cex=1.25,
         legend = expression(
           italic(d[max]) == 1,
           italic(d[max]) == 2,
           italic(d[max]) == 3))
  
  plot(x = d, y = genDWF("exponential",1)(d), type = "l", ylim = c(0,1),
       las = 1L, xlab = "", ylab = "", cex.axis = 2, lwd = 2)
  lines(x = d, y = genDWF("exponential",0.5)(d), col = "red", lwd = 2)
  lines(x = d, y = genDWF("exponential",2)(d), col = "blue", lwd = 2)
  text(x = 1.5, y = 0.85, label = "Exponential", adj = 0, cex = 2)
  legend(x = 3.2, y = 0.85, lwd=1, col=c("black","red","blue"), cex=1.25,
         legend = expression(
           italic(d[max]) == 1,
           italic(d[max]) == 0.5,
           italic(d[max]) == 2))
  
  plot(x = d, y = genDWF("Gaussian",1)(d), type = "l", ylim = c(0,1),
       las = 1L, xlab = "", ylab = "", cex.axis = 2, lwd = 2)
  lines(x = d, y = genDWF("Gaussian",0.5)(d), col = "red", lwd = 2)
  lines(x = d, y = genDWF("Gaussian",2)(d), col = "blue", lwd = 2)
  text(x = 1.5, y = 0.85, label = "Gaussian", adj = 0, cex = 2)
  legend(x = 3.2, y = 0.85, lwd=1, col=c("black","red","blue"), cex=1.25,
         legend = expression(
           italic(d[max]) == 1,
           italic(d[max]) == 0.5,
           italic(d[max]) == 2))
  
  plot(x = d, y = genDWF("hole_effect",1)(d), type = "l", ylim = c(-0.2,1),
       las = 1L, xlab = "", ylab = "", cex.axis = 2, lwd = 2)
  lines(x = d, y = genDWF("hole_effect",0.5)(d), col = "red", lwd = 2)
  lines(x = d, y = genDWF("hole_effect",2)(d), col = "blue", lwd = 2)
  abline(h = 0, lty = 3L)
  text(x = 1.5, y = 0.85, label = "Hole effect", adj = 0, cex = 2)
  legend(x = 3.2, y = 0.85, lwd=1, col=c("black","red","blue"), cex=1.25,
         legend = expression(
           italic(d[max]) == 1,
           italic(d[max]) == 0.5,
           italic(d[max]) == 2))
  
  dev.off()
  
  if(FALSE) {
    par(mar=c(5.1,5.1,0.6,0.6))
    d <- seq(0,5,0.001)
    plot(x = d, y = genDWF("hyperbolic",1,1)(d), type = "l", ylim = c(0,1),
         las = 1L, xlab = "", ylab = "", cex.axis = 2, lwd = 2)
    lines(x = d, y = genDWF("hyperbolic",2,0.5)(d), col = "red", lwd = 2)
    lines(x = d, y = genDWF("hyperbolic",0.5,2)(d), col = "blue", lwd = 2)
    text(x = 2, y = 0.9, cex = 2,
         label = "Hyperbolic", adj=0)
    rm(d)
  }
  
  rm(d)
}

## One-dimensional spatial eigenfunction profile
if(FALSE) {
  
  n <- 11
  
  profile <- list()
  
  if(!require(parallel)) break
  cl <- makeForkCluster(4L)
  
  ## i = "equidistant"
  ## i = "random"
  for(i in c("equidistant","random")) {
    if(i == "equidistant") {
      x <- (n - 1)*seq(0, 1, length.out=n)
    } else if(i == "random") {
      x <- round(c(0, cumsum(runif(n - 1L, 0.5, 1.5))),2)
      x <- (n - 1)*x/max(x)
    }
    
    ## x <- profile[[i]]$x
    profile[[i]] <- list(x = x)
    
    ## j = "hyperbolic"
    for(j in c("linear","power","hyperbolic","spherical","exponential",
               "Gaussian","hole_effect")) {
      switch(
        j,
        linear = draw(
          x, genDistMetric(), genDWF(j, 5), 0.1, 0.001
        ),
        power = draw(
          x, genDistMetric(), genDWF(j, 5, 0.5), 0.1, 0.001
        ),
        hyperbolic = draw(
          x, genDistMetric(), genDWF(j, 5, 1), 0.1, 0.001
        ),
        spherical = draw(
          x, genDistMetric(), genDWF(j, 5), 0.1, 0.001
        ),
        exponential = draw(
          x, genDistMetric(), genDWF(j, 5), 0.1, 0.001
        ),
        Gaussian = draw(
          x, genDistMetric(), genDWF(j, 5), 0.1, 0.001
        ),
        hole_effect = draw(
          x, genDistMetric(), genDWF(j, 5), 0.1, 0.001
        )
      ) -> profile[[i]][[j]]$drw
      getDerivatives(profile[[i]][[j]]$drw, cl) -> profile[[i]][[j]]$drv
    }
  }
  
  stopCluster(cl)
  
  rm(n,i,j,x,cl)
  
  for(i in names(profile$equidistant))
    sprintf("../../Data/profile-EQ-%s.rds",i) %>%
    saveRDS(profile$equidistant[[i]], file=.)
  
  for(i in names(profile$random))
    sprintf("../../Data/profile-RD-%s.rds",i) %>%
    saveRDS(profile$random[[i]], file=.)
  
  rm(i)
  
} else {
  
  profile <- list(equidistant = list(), random = list())
  
  for(i in c("x","linear","power","hyperbolic","spherical","exponential",
             "Gaussian","hole_effect")) {
    
    sprintf("../../Data/profile-EQ-%s.rds",i) %>%
      readRDS -> profile$equidistant[[i]]
    
    sprintf("../../Data/profile-RD-%s.rds",i) %>%
      readRDS -> profile$random[[i]]
    
  }
  
  rm(i)
  
}

## Code to display the spatial eigenfunctions four-by-four:
if(FALSE) {
  
  if(dev.cur() > 1) dev.off()
  
  pts <- "equidistant"  ## pts <- "random"  ##
  
  fun <- "hyperbolic"
  
  n <- length(profile[[pts]][[fun]]$drw$sef$getLambda())
  
  for(i in 1L:n) {
    m <- i + (0L:3L)
    m[m>n] - n -> m[m>n]
    plotDraw1D(drw=profile[[pts]][[fun]]$drw, m=m)
  }
  
  rm(pts,fun,n,i,m)
  
}

## Code to display each eigenfunction with its first three derivatives:
if(FALSE) {
  
  if(dev.cur() > 1) dev.off()
  
  pts <- "equidistant"  ## pts <- "random"
  
  fun <- "hyperbolic"
  
  n <- length(profile[[pts]][[fun]]$drw$sef$getLambda())
  
  for(i in 1L:n)
    plotDerivatives(
      drw = profile[[pts]][[fun]]$drw,
      drv = profile[[pts]][[fun]]$drv,
      j=i
    )
  
  rm(pts,fun,n,i)
  
}

rm(profile)

## Three-point one-dimensional spatial eigenfunction profiles
if(FALSE) {
  
  profile3pts <- list()
  
  if(!require(parallel)) break
  
  cl <- makeForkCluster(4L)
  
  ## j="hyperbolic"
  for(j in c("linear","power","hyperbolic","spherical","exponential",
             "Gaussian","hole_effect")) {
    list(
      x = seq(0.04,0.96,0.02),
      drw = list(),
      drv = list()
    ) -> profile3pts[[j]]
    
    ## i=1L
    for(i in 1L:length(profile3pts[[j]]$x)) {
      switch(
        j,
        linear = draw(
          c(0,profile3pts[[j]]$x[i],1), genDistMetric(), genDWF(j, 1),
          0.1, 0.001
        ),
        power = draw(
          c(0,profile3pts[[j]]$x[i],1), genDistMetric(), genDWF(j, 1, 0.5),
          0.1, 0.001
        ),
        hyperbolic = draw(
          c(0,profile3pts[[j]]$x[i],1), genDistMetric(), genDWF(j, 1, 1),
          0.1, 0.001
        ),
        spherical = draw(
          c(0,profile3pts[[j]]$x[i],1), genDistMetric(), genDWF(j, 1),
          0.1, 0.001
        ),
        exponential = draw(
          c(0,profile3pts[[j]]$x[i],1), genDistMetric(), genDWF(j, 1),
          0.1, 0.001
        ),
        Gaussian = draw(
          c(0,profile3pts[[j]]$x[i],1), genDistMetric(), genDWF(j, 1),
          0.1, 0.001
        ),
        hole_effect = draw(
          c(0,profile3pts[[j]]$x[i],1), genDistMetric(), genDWF(j, 1),
          0.1, 0.001
        )
      ) -> profile3pts[[j]]$drw[[i]]
      getDerivatives(profile3pts[[j]]$drw[[i]], cl) -> profile3pts[[j]]$drv[[i]]
    }
  }
  
  stopCluster(cl)
  
  rm(cl,i,j)
  
  for(i in names(profile3pts))
    sprintf("../../Data/profile3pts-%s.rds",i) %>%
    saveRDS(profile3pts[[i]], file=.)
  
  rm(i)
  
} else {
  
  profile3pts <- list()
  
  for(i in c("linear","power","hyperbolic","spherical","exponential",
             "Gaussian","hole_effect"))
    sprintf("../../Data/profile3pts-%s.rds",i) %>%
    readRDS -> profile3pts[[i]]
  
  rm(i)
  
}

## Three-point one-dimensional spatial eigenfunction profile plotting script:
if(FALSE) {
  
  if(dev.cur() > 1) dev.off()
  
  fun <- "hyperbolic"
  
  n <- length(profile3pts[[fun]]$x)
  
  for(i in 1L:n) {
    plotDerivatives3Pts(
      drw = profile3pts[[fun]]$drw[[i]],
      drv = profile3pts[[fun]]$drv[[i]]
    )
  }
  
  rm(fun,n,i)
  
}

rm(profile3pts)

## two-dimensional spatial eigenfunction profiles
if(FALSE) {
  
  profile2D <- list()
  
  c( 0.4172758, 0.03786392,
     0.1963081,-0.32355223,
     0.5910233, 0.08128149,
    -0.1712318,-0.86602540,
     0.7605802, 0.86602540,
    -1.0000000,-0.51367787,
     1.0000000,-0.06088640) %>%
    matrix(
      ncol=2L, byrow=TRUE,
      dimnames=list(NULL,c("x","y"))
    ) -> rnd.coords
  
  if(!require(parallel)) break
  
  cl <- makeForkCluster(2L)
  
  pids <- unlist(clusterCall(cl, Sys.getpid))  ## If one needs to kill them
  
  ## i="equidistant"
  ## i="random"
  for(i in c("equidistant","random")) {
    cat(i,'\n')
    if(i == "equidistant") {
      cbind(
        x = c(-0.5,0.5,-1,0,1,-0.5,0.5),
        y = c(rep(sqrt(3)/2,2L),rep(0,3L),rep(-sqrt(3)/2,2L))
      ) -> coords
    } else if(i == "random") {
      coords <- rnd.coords
    }
    
    profile2D[[i]] <- list(coords = coords)
    ## coords <- profile2D[[i]]$coords
    
    ## j="hyperbolic"
    for(j in c("linear","power","hyperbolic","spherical","exponential",
               "Gaussian","hole_effect")) {
      cat(j,"- ")
      switch(
        j,
        linear = draw2D(
          x=coords, m=genDistMetric(), f=genDWF(j,2), ext=0.25, by=0.025
        ),
        power = draw2D(
          x=coords, m=genDistMetric(), f=genDWF(j,2,0.5), 0.25, 0.025
        ),
        hyperbolic = draw2D(
          x=coords, m=genDistMetric(), f=genDWF(j,2,1), 0.25, 0.025
        ),
        spherical = draw2D(
          x=coords, m=genDistMetric(), f=genDWF(j,2), 0.25, 0.025
        ),
        exponential = draw2D(
          x=coords, m=genDistMetric(), f=genDWF(j,2), 0.25, 0.025
        ),
        Gaussian = draw2D(
          x=coords, m=genDistMetric(), f=genDWF(j,2), 0.25, 0.025
        ),
        hole_effect = draw2D(
          x=coords, m=genDistMetric(), f=genDWF(j,2), 0.25, 0.025
        )
      ) -> profile2D[[i]][[j]]$drw
      getDerivatives2D(profile2D[[i]][[j]]$drw, cl) -> profile2D[[i]][[j]]$drv
    }
  }
  
  stopCluster(cl)
  
  rm(rnd.coords,i,j,coords,cl)
  
  ## system(sprintf("kill -9 %s",paste(pids,collapse=" ")))
  ## rm(pids)
  
  for(i in names(profile2D$equidistant))
    sprintf("../../Data/profile2D-EQ-%s.rds",i) %>%
    saveRDS(profile2D$equidistant[[i]], file=.)
  
  for(i in names(profile2D$random))
    sprintf("../../Data/profile2D-RD-%s.rds",i) %>%
    saveRDS(profile2D$random[[i]], file=.)
    
  rm(i)
  
} else {
  
  profile2D <- list(equidistant = list(), random = list())
  
  for(i in c("coords","linear","power","hyperbolic","spherical",
             "exponential","Gaussian","hole_effect")) {
    sprintf("../../Data/profile2D-EQ-%s.rds",i) %>%
      readRDS -> profile2D$equidistant[[i]]
    sprintf("../../Data/profile2D-RD-%s.rds",i) %>%
      readRDS -> profile2D$random[[i]]
    
  }
  
  rm(i)
  
}

if(FALSE) {
  
  if(dev.cur() > 1) dev.off()
  
  pts <- "random"
  
  fun <- "hyperbolic"
  
  n <- ncol(profile2D[[pts]][[fun]]$drw$scr)
  
  for(i in 1L:n)
    plotDerivatives2D(
      drw = profile2D[[pts]][[fun]]$drw,
      drv = profile2D[[pts]][[fun]]$drv,
      j = i
    )
  
  rm(pts,fun,n,i)
  
}

rm(profile2D)

kable(
  data.frame(
    Name = c("Linear","Concave up","Concave down"),
    Definition = c(
      "$a_{i,j} = 1-\\frac{d_{i,j}}{d_{max}}$",
      "$a_{i,j} = 1-\\left(\\frac{d_{i,j}}{d_{max}}\\right)^\\alpha$",
      "$a_{i,j} = \\frac{1}{d_{i,j}^\\alpha}$"
    ),
    Ref = c("(T1 1)","(T1 2)","(T1 3)")
  ),
  align = "llr",col.names = c("Name","Definition","")
)

kable(
  data.frame(
    Name = sprintf("$\\mathrm{%s}$",c("Spherical","Exponential","Gaussian","Hole effect")),
    Definition = c(
      "$w_i = \\begin{cases}
          d_i < d_{max}, 1 - 1.5\\left(\\frac{d_i}{d_{max}}\\right) +  0.5\\left(\\frac{d_i}{d_{max}}\\right)^3 \\\\
          d_i \\geq d_{max}, 0
        \\end{cases}$",
      "$w_i = \\mathrm{e}^{-\\frac{d_i}{d_{max}}}$",
      "$w_i = \\mathrm{e}^{-\\left(\\frac{d_i}{d_{max}}\\right)^2}$",
      "$w_i = \\begin{cases}
          d_i = 0, 1 \\\\
          d_i > 0, \\frac{d_{max}}{\\pi d_i}\\sin \\frac{\\pi d_i}{d_{max}}
        \\end{cases}.$"
    ),
    Ref = sprintf("$\\mathrm{%s}$",c("(T2 1)","(T2 2)","(T2 3)","(T2 4)"))
  ),
  align = "llr",col.names = c("Name","Definition","")
)
