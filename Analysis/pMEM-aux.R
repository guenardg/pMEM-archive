## 
##  Auxiliary functions
##
draw <- function(x, m, f, ext, by, tol = .Machine$double.eps^0.5) {
  
  sef <- genSEF(x, m, f, tol)
  
  xx <- seq(min(x) - ext, max(x) + ext, by)
  
  scr <- sef$getPredictor(xx)
  
  list(x = x, sef = sef, by = by, xx = xx, scr = scr)
  
}

getDerivatives <- function(drw, cl, w = 13L, ...) {
  
  mm <- seq(0, (w - 1L)*drw$by, drw$by)
  
  mm <- mm - mean(mm)
  
  mm <- cbind(intercept=1, x=mm, x2=mm^2, x3=mm^3)
  
  res <- list()
  for(j in 1L:ncol(drw$scr)) {
    parSapply(
      cl = cl,
      X = 1L:(length(drw$xx) - w + 1L),
      FUN = function(i, y, mm, w, ...) {
        lm.fit(x = mm, y = y[i + (1L:w) - 1L], ...) -> lm0
        c(
          lm0$coefficients[2L],
          2*lm0$coefficients[3L],
          6*lm0$coefficients[4L]
        )
      },
      y = drw$scr[,j] * sign(drw$scr[1L,j]),
      mm = mm,
      w = w,
      ...
    ) -> res[[j]]
    
    as.data.frame(
      rbind(
        matrix(NA,w %/% 2L, 3L),
        t(res[[j]]),
        matrix(NA,w %/% 2L, 3L)
      )
    ) -> res[[j]]
    
  }
  
  res
  
}

plotDraw1D <- function(drw, m, ylab, xlab = "Location", ...) {
  
  par(no.readonly = TRUE) -> tmp
  
  ev <- drw$sef$getIMoran()
  
  drw$sef$getSEF(m)
  
  ylim <- max(abs(drw$sef$getSEF(m)),abs(drw$scr[,m])) * c(-1,1)
  
  if(missing(ylab))
    ylab <- sprintf("MEM %02d",m)
  
  par(mfrow=c(4,1), mar=c(4.1,4.1,0.5,2.0))
  for(i in 1L:length(m)) {
    sign(drw$scr[1,m[i]]) -> s
    plot(x=drw$xx, y=s*drw$scr[,m[i]], type="l", ylim=ylim, las=1L,
         ylab=ylab[i], xlab="", ...)
    abline(h=0, lty=3L)
    drw$sef$getSEF(m[i])
    points(x=drw$x, y=s*drw$sef$getSEF(m[i]))
    mtext(round(ev[m[i]],4L), 4L, 0.5, FALSE, 0)
  }
  
  mtext(text=xlab, side=1L, line=-1.5, outer=TRUE)
  
  par(tmp)
  
}

plotDerivatives <- function(drw, drv, xlab = "Location", ..., j=1L) {
  
  par(no.readonly = TRUE) -> tmp
  
  par(mfrow=c(4,1), mar=c(4.1,4.1,0.5,2.0))
  
  sign(drw$scr[1,j]) -> s
  
  plot(x = drw$xx, y = s*drw$scr[,j], type = "l", las = 1L, xlab="",
       ylab=sprintf("MEM %02d", j), ...)
  
  abline(v = drw$x, lty = 3L)
  
  abline(h = 0, lty = 2L)
  
  plot(x = drw$xx, y = drv[[j]][,1L], type="l", las = 1L, xlab="",
       ylab = expression(paste(1^{st}," derivative")),
       ylim = max(abs(drv[[j]][,1L]),na.rm=TRUE)*c(-1,1), ...)
  
  abline(v = drw$x, lty = 3L)
  
  abline(h = 0, lty = 2L)
  
  plot(x = drw$xx, y = drv[[j]][,2L], type="l", las = 1L, xlab="",
       ylab = expression(paste(2^{nd}," derivative")),
       ylim = max(abs(drv[[j]][,2L]),na.rm=TRUE)*c(-1,1), ...)
  
  abline(v = drw$x, lty = 3L)
  
  abline(h = 0, lty = 2L)
  
  plot(x = drw$xx, y = drv[[j]][,3L], type="l", las = 1L, xlab="",
       ylab = expression(paste(3^{rd}," derivative")),
       ylim = max(abs(drv[[j]][,3L]),na.rm=TRUE)*c(-1,1), ...)
  
  abline(v = drw$x, lty = 3L)
  
  abline(h = 0, lty = 2L)
  
  mtext(text=xlab, side=1L, line=-1.5, outer=TRUE)
  
  par(tmp)
  
  invisible(NULL)
  
}

plotDerivatives3Pts <- function(drw, drv, xlab = "Location", ...) {
  
  par(no.readonly = TRUE) -> tmp
  par(mfrow=c(4,2), mar=c(4.1,5.1,0.5,0.5))
  
  for(j in 1L:2L) {
    
    sign(drw$scr[1,j]) -> s
    
    plot(x = drw$xx, y = s*drw$scr[,j], type = "l", las = 1L, xlab="",
         ylab = if(j==1L) "MEM 1" else "MEM 2", ...)
    
    abline(v = drw$x, lty = 3L)
    
    abline(h = 0, lty = 2L)
    
  }
  
  for(j in 1L:2L) {
    
    plot(x = drw$xx, y = drv[[j]][,1L], type="l", las = 1L, xlab="",
         ylab = if(j==1L) "First derivative\n" else "",
         ylim = max(abs(drv[[j]][,1L]),na.rm=TRUE)*c(-1,1), ...)
    
    abline(v = drw$x, lty = 3L)
    
    abline(h = 0, lty = 2L)
    
  }
  for(j in 1L:2L) {
    
    plot(x = drw$xx, y = drv[[j]][,2L], type="l", las = 1L, xlab="",
         ylab = if(j==1L) "Second derivative\n" else "",
         ylim = max(abs(drv[[j]][,2L]),na.rm=TRUE)*c(-1,1), ...)
    
    abline(v = drw$x, lty = 3L)
    
    abline(h = 0, lty = 2L)
    
  }
  for(j in 1L:2L) {
    
    plot(x = drw$xx, y = drv[[j]][,3L], type="l", las = 1L, xlab="",
         ylab = if(j==1L) "Third derivative\n" else "",
         ylim = max(abs(drv[[j]][,3L]),na.rm=TRUE)*c(-1,1), ...)
    
    abline(v = drw$x, lty = 3L)
    
    abline(h = 0, lty = 2L)
    
  }
  
  mtext(text=xlab, side=1L, line=-1.5, outer=TRUE)
  
  par(tmp)
  
  invisible(NULL)
  
}

rerange <- function(x, min, max)
  (x - min(x)) * (max - min) / (max(x) - min(x)) + min

draw2D <- function(x, m, f, ext, by, tol = .Machine$double.eps^0.5) {
  
  sef <- genSEF(x, m, f, tol)
  
  seq(min(x[,1L]) - ext, max(x[,1L]) + ext, by) -> xx
  
  seq(min(x[,2L]) - ext, max(x[,2L]) + ext, by) -> yy
  
  list(
    x = xx,
    y = yy,
    coords = cbind(
      x = rep(xx, length(yy)),
      y = rep(yy, each = length(xx))
    )
  ) -> ss
  
  scr <- sef$getPredictor(ss$coords)
  
  list(x = x, sef = sef, by = by, ss = ss, scr = scr)
  
}

getDerivatives2D <- function(drw, cl, w = 13L, ..., verbose = TRUE) {
  
  mm <- seq(0, (w - 1L)*drw$by, drw$by)
  
  mm <- mm - mean(mm)
  
  cbind(
    x = rep(mm, each = length(mm)),
    y = rep(mm, length(mm))
  ) -> mm
  
  cbind(
    intercept=1, x=mm[,1L], y=mm[,2L], xx=mm[,1L]^2, yy=mm[,2L]^2,
    xy=mm[,1L]*mm[,2L], xxx=mm[,1L]^3, yyy = mm[,2L]^3,
    xxy=mm[,1L]^2*mm[,2L], xyy=mm[,1L]*mm[,2L]^2
  ) -> mm
  
  array(
    NA,
    c(ncol(drw$scr),length(drw$ss$x),length(drw$ss$y),6L),
    list(
      sprintf("SEV%02d",1L:ncol(drw$scr)),
      sprintf("X%d",1L:length(drw$ss$x)),
      sprintf("Y%d",1L:length(drw$ss$y)),
      c("dx","dy","d2x","d2y","d3x","d3y")
    )
  ) -> res
  
  ii <- (1L:(length(drw$ss$x) - w + 1L)) + w %/% 2
  
  for(k in 1L:ncol(drw$scr)) {
    
    matrix(
      drw$scr[,k] * sign(drw$scr[1L,k]),
      length(drw$ss$x),
      length(drw$ss$y)
    ) -> Y
    
    for(j in 1L:(length(drw$ss$y) - w + 1L)) {
      
      parSapply(
        cl = cl,
        X = 1L:(length(drw$ss$x) - w + 1L),
        FUN = function(i, y, mm, w, ...) {
          lm.fit(x = mm, y = as.numeric(y[i + (1L:w) - 1L,]), ...) -> lm0
          c(
            lm0$coefficients[2L],
            lm0$coefficients[3L],
            2*lm0$coefficients[4L],
            2*lm0$coefficients[5L],
            6*lm0$coefficients[7L],
            6*lm0$coefficients[8L]
          )
        },
        y = Y[,j + (1L:w) - 1L],
        mm = mm,
        w = w,
        ...
      ) -> tmp
      
      res[k,ii,j + w %/% 2,] <- t(tmp)
      
    }
    
    if(verbose) cat(k,' ')
    
  }
  
  if(verbose) cat('\n')
  
  res
  
}

plotDerivatives2D <- function(drw, drv, ..., j=1L,
                              col=head(rainbow(1200L),1000L)) {
  
  par(no.readonly = TRUE) -> tmp
  
  par(mfrow = c(2L,2L), mar=c(3.6,3.6,1.1,1.1))
  image(
    z=matrix(sign(drw$scr[1L,j])*drw$scr[,j], length(drw$ss$x),
             length(drw$ss$y)),
    x=drw$ss$x, y=drw$ss$y, xlab="", ylab="", asp=1, las=1L,
    col=col
  )
  
  points(drw$x)
  
  text(x=-1.15, y=1.0, label=parse(text=sprintf("bold(u)[%d]",j)),
       font=2L, xpd=TRUE, cex=1.25, adj=0)
  
  image(
    z=matrix((drv[j,,,1L]^2 + drv[j,,,2L]^2)^0.5, length(drw$ss$x),
             length(drw$ss$y)),
    x=drw$ss$x, y=drw$ss$y, xlab="", ylab="", asp=1, las=1L,
    col=col
  )
  
  points(drw$x)
  
  text(x=-1.15, y=1.0, label="First", font=2L, xpd=TRUE, adj=0)
  
  image(
    z=matrix((drv[j,,,3L]^2 + drv[j,,,4L]^2)^0.5, length(drw$ss$x),
             length(drw$ss$y)),
    x=drw$ss$x, y=drw$ss$y, xlab="", ylab="", asp=1, las=1L,
    col=col
  )
  
  points(drw$x)
  
  text(x=-1.15, y=1.0, label="Second", font=2L, xpd=TRUE, adj=0)
  
  image(
    z=matrix((drv[j,,,5L]^2 + drv[j,,,6L]^2)^0.5, length(drw$ss$x),
             length(drw$ss$y)),
    x=drw$ss$x, y=drw$ss$y, xlab="", ylab="", asp=1, las=1L,
    col=col
  )
  
  points(drw$x)
  
  text(x=-1.15, y=1.0, label="Third", font=2L, xpd=TRUE, adj=0)
  
  par(tmp)
  
  invisible(NULL)
  
}

genBrownianSurface <- function(n, k, mu0 = 0, sigma0 = 2, sigma = 0.01) {
  
  if(!require(spdep))
    stop("This function requires package 'spdep'")
  
  sqrt32 <- sqrt(3)/2
  
  cbind(
    x = rep(seq(-(n - 1L)/2, (n - 1L)/2,1), n) +
      rep(rep(c(-0.25,0.25), each = n), length.out=n^2),
    y = rep(seq(-(n - 1L)*sqrt32/2, (n - 1L)*sqrt32/2, sqrt32),
            each=n),
    z = NA
  ) -> grid
  
  lnk <- nb2listw(tri2nb(grid[,c("x","y")]), style="B")$neighbours
  
  attributes(lnk) <- NULL
  
  grid[sample(1L:nrow(grid), k, FALSE),"z"] <- rnorm(k, mu0, sigma0)
  
  while(any(is.na(grid[,"z"]))) {
    which(
      sapply(
        X = 1L:nrow(grid),
        FUN = function(x,y,z) is.na(y[x,"z"]) && any(!is.na(y[z[[x]],"z"])),
        y = grid,
        z = lnk
      )
    ) -> wh
    if(length(wh) > 1L)
      wh <- sample(wh, 1L)
    rnorm(
      n = 1L,
      mean = mean(grid[lnk[[wh]][!is.na(grid[lnk[[wh]],"z"])],"z"]),
      sd = sigma
    ) -> grid[wh,"z"]
  }
  
  grid
  
}

mapKernelSmoother <- function(map, sigma = c(1,10,100,1000)) {
  
  res <- matrix(NA, NROW(map), length(sigma))
  
  d <- as.matrix(dist(map[,1L:2L]))
  
  for(i in 1L:length(sigma)) {
    w <- exp(-d^2/sigma[i])
    apply(
      w,
      1L,
      function(x, y)
        sum(y*x)/sum(x),
      y = map[,"z"]
    ) -> res[,i]
  }
  
  res
  
}

getSEFpred <- function(x, xx, m, f, tol = .Machine$double.eps^0.5) {
  
  sef <- genSEF(x, m, f, tol)
  
  if(!missing(xx))
    scr <- sef$getPredictor(xx)
  
  list(U = sef$getSEF(), Up = if(!missing(xx)) scr)
  
}

objfsim <- function(par, xy, map, train, fun, getMod = FALSE) {
  
  switch(
    fun,
    linear = getSEFpred(x = xy[train,], xx = xy[-train,], m = genDistMetric(),
                        f = genDWF(fun,par[1L])),
    power = getSEFpred(x = xy[train,], xx = xy[-train,],
                       m = genDistMetric(),
                       f = genDWF(fun,par[1L],par[2L])),
    hyperbolic = getSEFpred(x = xy[train,], xx = xy[-train,],
                            m = genDistMetric(),
                            f = genDWF(fun,par[1L],par[2L])),
    spherical = getSEFpred(x = xy[train,], xx = xy[-train,],
                           m = genDistMetric(), f = genDWF(fun,par[1L])),
    exponential = getSEFpred(x = xy[train,], xx = xy[-train,],
                             m = genDistMetric(), f = genDWF(fun,par[1L])),
    Gaussian = getSEFpred(x = xy[train,], xx = xy[-train,], m = genDistMetric(),
                          f = genDWF(fun,par[1L])),
    hole_effect = getSEFpred(x = xy[train,], xx = xy[-train,],
                             m = genDistMetric(), f = genDWF(fun,par[1L]))
  ) -> tmp
  
  getMinMSE(
    U = tmp$U, y = map[train],
    Up = tmp$Up, yy = map[-train],
    complete = FALSE
  ) -> mod
  
  if(getMod) mod else mod$mse
  
}

getSEFpred2 <- function(x, xx, m, f, tol = .Machine$double.eps^0.5) {
  
  sef <- genSEF(x, m, f, tol)
  
  wh <- which(sef$getIMoran() > 0)
  
  if(!missing(xx))
    scr <- sef$getPredictor(xx,wh)
  
  list(U = sef$getSEF(wh), Up = if(!missing(xx)) scr)
  
}

objfsim2 <- function(par, xy, map, train, fun, getMod = FALSE) {
  
  switch(
    fun,
    linear=getSEFpred2(x=xy[train,], xx=xy[-train,], m=genDistMetric(),
                       f=genDWF(fun,par[1L])),
    power=getSEFpred2(x=xy[train,], xx=xy[-train,], m=genDistMetric(),
                      f=genDWF(fun,par[1L],par[2L])),
    hyperbolic=getSEFpred2(x=xy[train,], xx=xy[-train,], m=genDistMetric(),
                           f=genDWF(fun,par[1L],par[2L])),
    spherical=getSEFpred2(x=xy[train,], xx=xy[-train,], m=genDistMetric(),
                          f=genDWF(fun,par[1L])),
    exponential=getSEFpred2(x=xy[train,], xx=xy[-train,], m=genDistMetric(),
                            f=genDWF(fun,par[1L])),
    Gaussian=getSEFpred2(x=xy[train,], xx=xy[-train,], m=genDistMetric(),
                         f=genDWF(fun,par[1L])),
    hole_effect=getSEFpred2(x=xy[train,], xx=xy[-train,], m=genDistMetric(),
                            f=genDWF(fun,par[1L]))
  ) -> tmp
  
  if(ncol(tmp$U)) {
    getMinMSE(
      U = tmp$U, y = map[train],
      Up = tmp$Up, yy = map[-train],
      complete = FALSE
    ) -> mod
    
  } else {
    list(
      betasq = 1,
      mse = mean((mean(map[train]) - map[-train])^2)
    ) -> mod
  }
  
  if(getMod) mod else mod$mse
  
}

sshPortRedirect <- function(host, user, port, port2) {
  
  system(
    sprintf(
      "ssh -L %d:127.0.0.1:%d -N %s@%s &",
      port2, port, user, host
    )
  )
  
  list(host=host, user=user, port=port, port2=port2)
  
}

sshPortShutdown <- function(conn)
  with(conn,{
    proc <- system("ps aux | grep ssh",TRUE)
    unlist(
      lapply(
        strsplit(proc," +"),
        function(x,h,u,p1,p2)
          (x[11L] == "ssh") &&
          (x[12L] == "-L") &&
          (x[13L] == sprintf("%d:127.0.0.1:%d",p2,p1)) &&
          (x[14L] == "-N") &&
          (x[15L] == sprintf("%s@%s",user,host)),
        h = host, u = user, p1 = port, p2 = port2
      )
    ) -> wh
    
    if(sum(wh)) {
      unlist(
        lapply(
          strsplit(proc[wh]," +"),
          function(x) x[2L]
        )
      ) -> pid
      system(paste("kill", paste(pid, collapse=" ")))
    } else
      warning("No such ssh TCP port redirection process was found\n\nUser: ",
              user,"\nHost: ",host,"\nHost port: ",port,"\nRecipient port: ",
              port2,"\n")
  })

showProcessingJobs <- function(sim)
  dbGetQuery(
    sim,
    "SELECT id, map, sample, n, fun, worker
     FROM results
     WHERE worker IS NOT NULL AND worker <> 'done'"
  )

showCompletedJobs <- function(sim) {
  
  unlist(
    dbGetQuery(
      sim,
      "SELECT COUNT(*)
       FROM results
       WHERE worker IS NOT NULL AND worker = 'done'"
    )
  ) -> done
  
  unlist(
    dbGetQuery(
      sim,
      "SELECT COUNT(*)
       FROM results"
    )
  ) -> total
  
  cat(
    sprintf(
      "%.2f%% done (%d of %d; remaining: %d)\n",
      100*done/total,done,total,total - done
    )
  )
  
  dbGetQuery(
    sim,
    "SELECT *
       FROM results
       WHERE worker IS NOT NULL AND worker = 'done'
       ORDER BY id"
  ) -> jobs
  
  invisible(jobs)
  
}

clearFailedJobs <- function(sim, worker) {
  
  dbGetQuery(
    sim,
    "SELECT id, worker
     FROM results
     WHERE worker IS NOT NULL and worker <> 'done'"
  ) -> tmp
  
  wh <- which(unlist(lapply(strsplit(tmp$worker,":"), head, n=1L)) == worker)
  
  dbExecute(
    sim,
    sprintf(
      "UPDATE results
       SET worker = null, msd = null, dmax = null, expn = null,
       betasq = null, mse = null
       WHERE id IN (%s)",
      paste(tmp$id[wh], collapse=",")
    )
  )
  
}

objf <- function(par, dwf, foldid, pts, y, z, ..., verbose = FALSE,
                 minMoran = 0) {
  
  range <- par[1L]
  
  shape <- par[2L]
  
  pMEM::genSEF(
    x = pts,
    m = pMEM::genDistMetric(),
    f = pMEM::genDWF(dwf, range=range, shape=shape)
  ) -> sef
  
  wh <- which(sef$getIMoran() > minMoran)
  
  cv.glmnet(
    x = cbind(if(!missing(z)) z, sef$getSEF(wh)),
    y = y,
    foldid = foldid,
    ...
  ) -> mod
  
  cvm <- mod$cvm[mod$index[1L,1L]]
  
  lambda <- mod$lambda[mod$index[1L,1L]]
  
  if(verbose)
    cat(
      sprintf(
        "dwf: %s, lambda: %f, range: %f, shape: %f, cvm: %f\n",
        dwf, lambda, range, shape, cvm
      )
    )
  
  cvm
  
}

getMod <- function(par, dwf, foldid, pts, y, z, ..., minMoran = 0) {
  
  range <- par[1L]
  
  shape <- par[2L]
  
  genSEF(
    x = pts,
    m = genDistMetric(),
    f = genDWF(dwf, range=range, shape=shape)
  ) -> sef
  
  wh <- which(sef$getIMoran() > minMoran)
  
  cv.glmnet(
    x = cbind(if(!missing(z)) z, sef$getSEF(wh)),
    y = y,
    foldid = foldid,
    ...
  ) -> mod
  
  lambda <- mod$lambda[mod$index[1L,1L]]
  
  list(
    sef = sef,
    model = mod,
    par = list(dwf=dwf, lambda=lambda, range=range, shape=shape, wh=wh)
  )
  
}

objfmv <- function(par, dwf, foldid, pts, Y, Z, ..., verbose = FALSE,
                   minMoran = 0) {
  
  range <- par[1L]
  
  shape <- par[2L]
  
  pMEM::genSEF(
    x = pts,
    m = pMEM::genDistMetric(),
    f = pMEM::genDWF(dwf, range=range, shape=shape)
  ) -> sef
  
  wh <- which(sef$getIMoran() > minMoran)
  
  DE <- cbind(const=1, if(!missing(Z)) Z, sef$getSEF(wh))
  
  SP <- cbind(1, contr.treatment(ncol(Y)))
  
  dimnames(SP) <- list(colnames(Y),colnames(Y))
  
  XX <- kronecker(X=SP, Y=DE, make.dimnames=TRUE)[,-1L]
  
  cv.glmnet(
    x = XX,
    y = as.numeric(Y),
    foldid = rep(foldid, NCOL(Y)),
    ...
  ) -> mod
  
  cvm <- mod$cvm[mod$index[1L,1L]]
  
  lambda <- mod$lambda[mod$index[1L,1L]]
  
  if(verbose)
    cat(
      sprintf(
        "dwf: %s, lambda: %f, range: %f, shape: %f, cvm: %f\n",
        dwf, lambda, range, shape, cvm
      )
    )
  
  cvm
  
}

getModmv <- function(par, dwf, foldid, pts, Y, Z, ..., verbose = FALSE,
                     minMoran = 0) {
  
  range <- par[1L]
  
  shape <- par[2L]
  
  pMEM::genSEF(
    x = pts,
    m = pMEM::genDistMetric(),
    f = pMEM::genDWF(dwf, range=range, shape=shape)
  ) -> sef
  
  wh <- which(sef$getIMoran() > minMoran)
  
  DE <- cbind(const=1, if(!missing(Z)) Z, sef$getSEF(wh))
  
  SP <- cbind(1, contr.treatment(ncol(Y)))
  
  dimnames(SP) <- list(colnames(Y),colnames(Y))
  
  XX <- kronecker(X=SP, Y=DE, make.dimnames=TRUE)[,-1L]
  
  cv.glmnet(
    x = XX,
    y = as.numeric(Y),
    foldid = rep(foldid, NCOL(Y)),
    ...
  ) -> mod
  
  lambda <- mod$lambda[mod$index[1L,1L]]
  
  predict(
    object = mod,
    newx = XX,
    s = "lambda.min",
    y = as.numeric(Y),
    type = "response",
    exact = TRUE
  ) -> fit
  
  dim(fit) <- dim(Y)
  
  dimnames(fit) <- dimnames(Y)
  
  list(
    sef = sef,
    model = mod,
    par = list(dwf=dwf, lambda=lambda, range=range, shape=shape, wh=wh),
    sp = SP,
    fit = fit
  )
  
}

plotMap <- function(conn, m, s, n, col, ..., smooth = 0, xlab="X", ylab="Y",
                    axes=FALSE, col.train = "black") {
  
  if(!require(base64enc))
    stop("Needs package 'base64enc'")
  
  if(missing(col))
    col <- rainbow(1200L)[1L:1000L]
  
  sprintf(
    "SELECT points.*, data.value AS value
     FROM data
     JOIN points ON data.point = points.id
     WHERE data.map = %d AND data.smooth = %d
     ORDER BY points.id",
    m, smooth
  ) %>%
    dbGetQuery(conn,.) -> dat
  
  if(nrow(dat)) {
    rng <- range(dat$value, na.rm=TRUE)
    dat$value[!is.na(dat$value)] %>%
      {1L + floor(999*(. - rng[1L])/(rng[2L] - rng[1L]))} -> idx
    plot(NA, xlim=range(dat[!is.na(dat$value),"x"]),
         ylim=range(dat[!is.na(dat$value),"y"]), asp=1, xlab=xlab,
         ylab=ylab, axes=axes, ...)
    seq(0.5*pi, 2.5*pi, length.out=7L) %>%
      {0.5*cbind(cos(.),sin(.))} -> pl
    #' i=1L
    for(i in 1L:nrow(dat))
      polygon(x=dat$x[i] + pl[,1L], y=dat$y[i] + pl[,2L],
              col=col[idx[i]], border=col[idx[i]])
    if(!missing(s) && !missing(n)) {
      sprintf(
        "SELECT encode(data,'base64') FROM samples WHERE id = %d",
        s
      ) %>%
        dbGetQuery(conn,.) %>%
        unlist %>%
        base64decode %>%
        unserialize -> smp
      if(n > 500L)
        stop("The maximum value of argument 'n' is 500.")
      train <- smp[1L:n]
      for(i in train)
        polygon(x=dat$x[i] + pl[,1L], y=dat$y[i] + pl[,2L],
                col="transparent", border=col.train, ...)
    }
  }
  
  invisible(
    sprintf(
      "SELECT * FROM results
       WHERE results.map = %d AND results.smooth = %d%s%s
       ORDER BY map, smooth, sample, n",
        m, smooth,
      if(!missing(s)) {
        sprintf("AND results.sample = %d", s)
      } else "",
      if(!missing(n)) {
        sprintf("AND results.n = %d", n)
      } else ""
      ) %>%
        dbGetQuery(conn,.)
  )
  
}

plotResults <- function(conn, r, col, ..., xlab="X", ylab="Y",
                        axes=FALSE, col.train = "black") {
  
  if(missing(col))
    col <- rainbow(1200L)[1L:1000L]
  
  sprintf(
    "SELECT *
     FROM results
     WHERE id = %d",
    r
  ) %>%
    dbGetQuery(sim,.) %>%
    as.list -> rr
  
  if(is.na(rr$msd))
    stop("No results available for r = ",r)
  
  sprintf(
    "SELECT points.x, points.y, data.value AS value
     FROM data
     JOIN points ON data.point = points.id
     WHERE data.map = %d AND data.smooth = %d
     ORDER BY points.id",
    rr$map, rr$smooth
  ) %>%
    dbGetQuery(conn,.) -> dat
  
  sprintf(
    "SELECT encode(data,'base64')
     FROM samples WHERE id = %d",
    rr$sample
  ) %>%
    dbGetQuery(conn,.) %>%
    unlist %>%
    base64decode %>%
    unserialize -> smp
  train <- smp[1L:rr$n]
  
  dat[train,1L:2L] %>%
    as.matrix %>%
    genSEF(
      m = genDistMetric(),
      f = genDWF(rr$fun, range=rr$dmax, shape=rr$expn)
    ) -> sef
  
  getPredsMinMSE <- function(x, y, xx, betasq) {
    ym <- mean(y)
    yn2 <- sum((y - ym)^2)
    b <- t(x) %*% (y - ym)
    b[(b^2/yn2) < betasq] <- 0
    ym + xx %*% b
  }
  
  getPredsMinMSE(
    x = sef$getSEF(),
    y = dat$value[train],
    xx = sef$getPredictor(as.matrix(dat[,1L:2L])),
    betasq = rr$betasq
  ) %>%
    as.numeric -> prd
  
  rng <- range(dat$value, prd, na.rm=TRUE)
  seq(0.5*pi, 2.5*pi, length.out=7L) %>%
    {0.5*cbind(cos(.),sin(.))} -> pl
  
  par(no.readonly = TRUE) -> parSafe
  par(mfrow = c(1L,3L), mar=c(4,4,2,2))
  
  dat$value[!is.na(dat$value)] %>%
    {1L + floor(999*(. - rng[1L])/(rng[2L] - rng[1L]))} -> idx
  dev.hold()
  
  plot(NA, xlim=range(dat[!is.na(dat$value),"x"]),
       ylim=range(dat[!is.na(dat$value),"y"]), asp=1, xlab="", ylab="",
       axes=axes,
       main = sprintf(
         "Map:%d, Sample:%d,\nSmooth:%d, N:%d",
         rr$map, rr$sample, rr$smooth, rr$n),
       ...)
  
  for(i in 1L:nrow(dat))
    polygon(x=dat$x[i] + pl[,1L], y=dat$y[i] + pl[,2L],
            col=col[idx[i]], border=col[idx[i]])
  
  for(i in train)
    polygon(x=dat$x[i] + pl[,1L], y=dat$y[i] + pl[,2L],
            col="transparent", border=col.train)
  
  prd %>%
    {1L + floor(999*(. - rng[1L])/(rng[2L] - rng[1L]))} -> idx
  
  plot(NA, xlim=range(dat[!is.na(dat$value),"x"]),
       ylim=range(dat[!is.na(dat$value),"y"]), asp=1, xlab="", ylab="",
       axes=axes, main=sprintf("Fun: %s\nDmax: %f",rr$fun, rr$dmax), ...)
  
  for(i in 1L:nrow(dat))
    polygon(x=dat$x[i] + pl[,1L], y=dat$y[i] + pl[,2L],
            col=col[idx[i]], border=col[idx[i]])
  
  res <- dat$value - prd
  
  rng <- range(res, na.rm=TRUE)
  
  res %>%
    {1L + floor(999*(. - rng[1L])/(rng[2L] - rng[1L]))} -> idx
  
  plot(NA, xlim=range(dat[!is.na(dat$value),"x"]),
       ylim=range(dat[!is.na(dat$value),"y"]), asp=1, xlab="", ylab="",
       axes=axes, ...)
  
  for(i in 1L:nrow(dat))
    polygon(x=dat$x[i] + pl[,1L], y=dat$y[i] + pl[,2L],
            col=col[idx[i]], border=col[idx[i]])
  dev.flush()
  
  mtext(xlab, 1L, 2, TRUE)
  mtext(ylab, 2L, 2, TRUE)
  
  par(parSafe)
  
  c(
    msd = mean((dat$value - mean(dat$value, na.rm=TRUE))^2,na.rm=TRUE),
    mse = mean((dat$value - prd)^2,na.rm=TRUE)
  )
  
}

getColors <- function(x, col, range)
  col[1L + floor((length(col) - 1L)*(x - range[1L])/(range[2L] - range[1L]))]

Psquare <- function(x, y) 1 - sum((x - y)^2)/sum((x - mean(x))^2)

prdOribatids <- function(
    sp,
    dwf = names(which.min(optOribatids$bestval)),
    desc_dwf = list(
      SubsDens=names(which.min(optSubsDens$bestval)),
      WatrCont=names(which.min(optWatrCont$bestval))
    ),
    s = "lambda.min",
    type = "response"
) {
  grid$Other %>%
    .[attr(.,"sf_column")] -> prd
  predict(
    object = modOribatids[[dwf]]$model,
    newx = kronecker(
      modOribatids[[dwf]]$sp[sp,,drop=FALSE],
      cbind(
        const = 1,
        grid$SubsDens[[desc_dwf$SubsDens]],
        grid$WatrCont[[desc_dwf$WatrCont]],
        as.matrix(st_drop_geometry(grid$Other)),
        modOribatids[[dwf]]$sef$getPredictor(
          st_coordinates(prd),
          modOribatids[[dwf]]$par$wh
        )
      ),
      make.dimnames = TRUE
    )[,-1L],
    s=s,
    type=type
  ) %>% as.numeric -> prd[[sp]]
  prd
}
