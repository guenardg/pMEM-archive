##
## Development of predictive Moran's eigenvector maps (pMEM)
## Oribatids example script
## 
## rm(list=ls())

library(magrittr)
library(raster)
library(sf)
library(pMEM)
library(glmnet)
library(parallel)
library(DEoptim)

source("../pMEM-aux.R")

data(geoMite)
attach(geoMite)

if(FALSE) {
  
  col <- list()
  col[["substrate"]] <- c(Sphagn1 = "#00ff00", Sphagn2 = "#fffb00",
                          Sphagn3 = "#774b00", Sphagn4 = "#ff8400",
                          Litter = "#ee00ff", Barepeat = "#ff0004")
  col[["water"]] <- c(Water = "#008cff", Flooded = "#ffffff00",
                      Predict = "#ffffff00")
  col[["shrub"]] <- c(None = "#dfdfdf", Few = "#a7a7a7", Many = "#5c5c5c")
  col[["topo"]] <- c(Blanket = "#74cd00", Hummock = "#bc9d00")
  
  ## X11(width=8.0, height=6.0)
  png(filename="../../Image/Oribatid_map.png", width=800, height=600)
  par(mar=c(1.1,2.1,1.1,0.1), mfrow=c(1L,4L))
  
  substrate %>%
    {plot(st_geometry(.), col=col[["substrate"]][.$Type], main="Substrate")}
  
  axis(1L, pos=0)
  
  axis(2L, pos=0)
  
  water[1L,] %>%
    {plot(st_geometry(.), col=col[["water"]][.$Type], add=TRUE)}
  
  water[-1L,] %>%
    {plot(st_geometry(.), col=col[["water"]][.$Type], lty=3L, add=TRUE)}
  
  core %>%
    {plot(st_geometry(.), pch = 21L, bg = "black", add=TRUE)}
  
  shrub %>%
    {plot(st_geometry(.), col = col[["shrub"]][.$Type], main="Shrubs")}
  
  axis(1L, pos=0)
  
  water[1L,] %>%
    {plot(st_geometry(.), col = col[["water"]][.$Type], add = TRUE)}
  
  water[-1L,] %>%
    {plot(st_geometry(.), col=col[["water"]][.$Type], lty=3L, add=TRUE)}
  
  core %>%
    {plot(st_geometry(.), pch = 21L, bg = "black", add=TRUE)}
  
  topo %>%
    {plot(st_geometry(.), col = col[["topo"]][.$Type], main="Topography")}
  
  axis(1L, pos=0)
  
  water[1L,] %>%
    {plot(st_geometry(.), col = col[["water"]][.$Type], add = TRUE)}
  
  water[-1L,] %>%
    {plot(st_geometry(.), col=col[["water"]][.$Type], lty=3L, add=TRUE)}
  
  core %>%
    {plot(st_geometry(.), pch = 21L, bg = "black", add=TRUE)}
  
  plot(NA, xlim=c(0,1), ylim=c(0,1), axes = FALSE)
  legend(x=0, y=0.9, pch=22L, pt.cex = 2.5, pt.bg=col[["substrate"]],
         box.lwd = 0, legend=names(col[["substrate"]]), title="Substrate")
  legend(x=0, y=0.6, pch=c(22L,NA,NA), pt.cex = 2.5, pt.bg=col[["water"]],
         box.lwd = 0, lty = c(0L,3L,2L),
         legend=c("Open water","Flooded area","Prediction polygon"))
  legend(x=0, y=0.4, pch=22L, pt.cex = 2.5, pt.bg=col[["shrub"]], box.lwd = 0,
         legend=names(col[["shrub"]]), title="Shrubs")
  legend(x=0, y=0.2, pch=22L, pt.cex = 2.5, pt.bg=col[["topo"]], box.lwd = 0,
         legend=names(col[["topo"]]), title="Topography")
  
  dev.off()
  
  rm(col)
}

## Make predictions for the two continuous variables
if(FALSE) {

  res <- c(x = 0.025, y = 0.025)
  
  list(
    x = seq(0,2.6,res["x"]),
    y = seq(0,10,res["y"])
  ) %>%
    {data.frame(
      x = rep(.$x, length(.$y)),
      y = rep(.$y, each = length(.$x))
    )} %>%
    st_as_sf(
      coords = c("x","y"),
      crs = "+proj=tmerc +units=m"
    ) %>%
    .[is.na(st_join(., water[1L,])$Type),] %>%
    {list(SubsDens=., WatrCont=., Other=.)} -> grid
  
  ## .[is.na(st_join(., water[1L,])$Type) &
  ##     !is.na(st_join(., hull[1L,])$Type),] %>%
  ##   {list(SubsDens=., WatrCont=., Other=.)} -> grid
  
  core %>%
    st_coordinates %>%
    dist %>%
    range
  
  FID <- rep(1L:7L, length.out=nrow(core))
  
  cl <- makeForkCluster(7L)
  ## stopCluster(cl) ; rm(cl)
  
  if(FALSE) {
    
    optSubsDens <- list(dwf = list())
    
    for(dwf in c("linear","power","hyperbolic","spherical","exponential",
                 "Gaussian","hole_effect"))
      DEoptim(
        fn = objf,
        lower = c(0.5,0.15),
        upper = c(5,1.85),
        dwf = dwf,
        foldid = FID,
        pts = st_coordinates(core),
        y = core$SubsDens,
        control = list(cluster=cl)
      ) -> optSubsDens$dwf[[dwf]]
    
    optSubsDens$dwf %>%
      lapply(function(x) x$optim$bestval) %>%
      unlist -> optSubsDens$bestval
    
    optSubsDens$dwf %>%
      lapply(function(x) x$optim$bestmem) %>%
      unlist %>%
      matrix(ncol=2L) -> optSubsDens$bestpar
    
    list(
      names(optSubsDens$dwf),
      c("range", "shape")
    ) -> dimnames(optSubsDens$bestpar)
    
    modSubsDens <- list()
    
    for(dwf in names(optSubsDens$dwf))
      getMod(
        par = optSubsDens$bestpar[dwf,],
        dwf = dwf,
        foldid = FID,
        pts = st_coordinates(core),
        y = core$SubsDens
      ) -> modSubsDens[[dwf]]
    
    rm(dwf)
    
    save(optSubsDens, modSubsDens, file="../../Data/SubsDens.rda")
    
  } else load(file="../../Data/SubsDens.rda")
  
  if(FALSE) {
    
    optWatrCont <- list(dwf = list())
    
    for(dwf in c("linear","power","hyperbolic","spherical","exponential",
                 "Gaussian","hole_effect"))
      DEoptim(
        fn = objf,
        lower = c(0.5,0.15),
        upper = c(5,1.85),
        dwf = dwf,
        foldid = FID,
        pts = st_coordinates(core),
        y = core$WatrCont,
        control = list(cluster=cl)
      ) -> optWatrCont$dwf[[dwf]]
    
    optWatrCont$dwf %>%
      lapply(function(x) x$optim$bestval) %>%
      unlist -> optWatrCont$bestval
    
    optWatrCont$bestval %>%
      which.min %>%
      names -> optWatrCont$bestdwf
    
    optWatrCont$dwf %>%
      lapply(function(x) x$optim$bestmem) %>%
      unlist %>%
      matrix(ncol=2L) -> optWatrCont$bestpar
    
    list(
      names(optWatrCont$dwf),
      c("range", "shape")
    ) -> dimnames(optWatrCont$bestpar)
    
    modWatrCont <- list()
    
    for(dwf in names(optWatrCont$dwf))
      getMod(
        par = optWatrCont$bestpar[dwf,],
        dwf = dwf,
        foldid = FID,
        pts = st_coordinates(core),
        y = core$WatrCont
      ) -> modWatrCont[[dwf]]
    
    rm(dwf)
    
    save(optWatrCont, modWatrCont, file="../../Data/WatrCont.rda")
    
  } else load(file="../../Data/WatrCont.rda")
  
  ## Generate the predictions: substrate density
  for(dwf in names(modSubsDens))
    predict(
      object = modSubsDens[[dwf]]$model,
      newx = modSubsDens[[dwf]]$sef$getPredictor(
        st_coordinates(grid$SubsDens),
        modSubsDens[[dwf]]$par$wh
      ),
      s = "lambda.min",
      type = "response"
    ) %>%
    as.numeric -> grid$SubsDens[[dwf]]
  
  ## Generate the predictions: water concentration
  for(dwf in names(modWatrCont))
    predict(
      object = modWatrCont[[dwf]]$model,
      newx = modWatrCont[[dwf]]$sef$getPredictor(
        st_coordinates(grid$WatrCont),
        modWatrCont[[dwf]]$par$wh
      ),
      s = "lambda.min",
      type = "response"
    ) %>%
    as.numeric -> grid$WatrCont[[dwf]]
  
  rm(dwf)
  
  colmap <- rainbow(1200L)[1L:1000L]
  
  ## X11(width=8, height=8)
  png("../../Image/SubsWatr.png", width=740, height=800)
  
  par(fig=c(0.0,0.2,0.0,1.0), mar=c(4,6.1,4,2))
  
  seq(
    min(core$SubsDens),
    max(core$SubsDens),
    length.out=1000L
  ) %>%
    matrix(1L,1000L) %>%
    image(
      zlim = range(core$SubsDens),
      col = colmap,
      axes = FALSE,
      ylab = expression(paste("Substrate density ",(g~L^{-1}))),
      cex.lab = 1.5
    )
  
  box()
  
  seq(20,80,10) %>%
    {axis(
      side = 2L,
      label = .,
      at = (. - min(core$SubsDens))/
        (max(core$SubsDens) - min(core$SubsDens)),
      cex.axis = 1.5
    )}
  
  par(fig=c(0.2,0.5,0.0,1.0), mar=c(2.5,1.5,2,0.5), new=TRUE)
  
  grid$SubsDens %>%
    {rasterFromXYZ(
      cbind(
        st_coordinates(.),
        .[[names(which.min(optSubsDens$bestval))]]
      ),
      crs = st_crs(.)
    )} %>%
    as.matrix %>% t %>% .[,ncol(.):1L] %>%
    image(
      x = seq(0,2.6,length.out=nrow(.)),
      y = seq(0,10,length.out=ncol(.)),
      zlim = range(core$SubsDens),
      col = colmap,
      axes = FALSE,
      asp = 1,
      xlab = "",
      ylab = ""
    )
  
  axis(1L)
  
  axis(2L)
  
  box()
  
  points(
    st_coordinates(core),
    pch=21L,
    bg = core$SubsDens %>%
      {(. - min(.))/(max(.) - min(.))} %>%
      {colmap[1L + round(999*.)]},
    cex=2
  )
  
  par(fig=c(0.5,0.7,0.0,1.0), mar=c(4,6.1,4,2), new=TRUE)
  
  seq(
    min(core$WatrCont),
    max(core$WatrCont),
    length.out=1000L
  ) %>%
    matrix(1L,1000L) %>%
    image(
      zlim = range(core$WatrCont),
      col = colmap,
      axes = FALSE,
      ylab = expression(paste("Water content ",(g~L^{-1}))),
      cex.lab = 1.5
    )
  
  box()
  
  seq(200,800,100) %>%
    {axis(
      side = 2L,
      label = .,
      at = (. - min(core$WatrCont))/
        (max(core$WatrCont) - min(core$WatrCont)),
      cex.axis = 1.5
    )}
  
  par(fig=c(0.7,1.0,0.0,1.0), mar=c(2.5,1.5,2,0.5), new=TRUE)
  
  grid$WatrCont %>%
    {rasterFromXYZ(
      cbind(
        st_coordinates(.),
        .[[names(which.min(optWatrCont$bestval))]]
      ),
      crs = st_crs(.)
    )} %>%
    as.matrix %>% t %>% .[,ncol(.):1L] %>%
    image(
      x = seq(0,2.6,length.out=nrow(.)),
      y = seq(0,10,length.out=ncol(.)),
      zlim = range(core$WatrCont),
      col = colmap,
      axes = FALSE,
      asp = 1,
      xlab = "",
      ylab = ""
    )
  
  axis(1L)
  
  box()
  
  points(
    st_coordinates(core),
    pch=21L,
    bg = core$WatrCont %>%
      {(. - min(.))/(max(.) - min(.))} %>%
      {colmap[1L + round(999*.)]},
    cex=2
  )
  
  dev.off()
  
  ## Plot one-by-one: substrate density
  grid$SubsDens %>%
    {rasterFromXYZ(
      cbind(
        st_coordinates(.),
        .[[names(which.min(optSubsDens$bestval))]]
      ),
      crs = st_crs(.)
    )} %>%
    plot(col=head(rainbow(1200L),1000L))
  
  ## Plot one-by-one: water content
  grid$WatrCont %>%
    {rasterFromXYZ(
      cbind(
        st_coordinates(.),
        .[[names(which.min(optWatrCont$bestval))]]
      ),
      crs = st_crs(.)
    )} %>%
    plot(col=head(rainbow(1200L),1000L))
  
}

if(FALSE) {
  
  dat <- list()
  
  core %>%
    .[,substr(colnames(.),1L,7L)=="Species"] -> dat$Y
  
  dat$Y %<>% st_set_geometry(NULL) %>% as.matrix
  
  core %>%
    .[,substr(colnames(.),1L,7L)!="Species"] -> dat$Z
  
  dat$Z %<>% st_set_geometry(NULL)
  
  dat$Z[,3L:8L] %<>% {./rowSums(.)}
  
  dat$shrubPoly <- poly(dat$Z$Shrub,2L)
  
  dat$Z %<>% cbind(dat$shrubPoly)
  
  dat$Z$Shrub <- NULL
  
  colnames(dat$Z)[11L:12L] <- c("Shrub^1","Shrub^2")
  
  dat$Z$Topo <- ifelse(
    dat$Z$Topo=="Blanket",
    1,
    -sum(dat$Z$Topo=="Blanket")/sum(dat$Z$Topo!="Blanket")
  )  ## sum(dat$Z$Topo)
  
  dat$Z$Flooded %<>% as.numeric
  
  dat$Z %<>% as.matrix
  
  saveRDS(dat, file="../../Data/Oribatids_dat.rds")
  
} else dat <- readRDS(file="../../Data/Oribatids_dat.rds")

## Cannot be run in RStudio; the PDE is not stable enough...
if(FALSE) {
  
  cl <- makeForkCluster(2L)
  ## stopCluster(cl) ; rm(cl)
  
  ## For testing only:
  objfmv(par=c(1,1), dwf="Gaussian", foldid=FID, pts=st_coordinates(core),
         Y=dat$Y, Z=dat$Z, family = "poisson", control=list(cluster=cl))
  
  optOribatids <- list(dwf = list())
  
  for(dwf in c("linear","power","hyperbolic","spherical","exponential",
               "Gaussian","hole_effect")) {
    cat(sprintf("DWF: %s\n",dwf))
    DEoptim(
      fn = objfmv,
      lower = c(0.5,0.15),
      upper = c(5,1.85),
      dwf = dwf,
      foldid = FID,
      pts = st_coordinates(core),
      Y = dat$Y,
      Z = dat$Z,
      family = "poisson",
      control = list(cluster=cl)
    ) -> optOribatids$dwf[[dwf]]
  }
  
  optOribatids$dwf %>%
    lapply(
      function(x)
        x$optim$bestval
    ) %>%
    unlist -> optOribatids$bestval
  
  optOribatids$dwf %>%
    lapply(
      function(x)
        x$optim$bestmem
    ) %>%
    unlist %>%
    matrix(ncol=2L, byrow=TRUE) -> optOribatids$bestmem
  
  list(
    names(optOribatids$dwf),
    c("range","shape")
  ) -> dimnames(optOribatids$bestmem)
  
  modOribatids <- list()
  
  for(dwf in names(optOribatids$dwf)) {
    getModmv(
      par = optOribatids$bestmem[dwf,],
      dwf = dwf,
      foldid = FID,
      pts = st_coordinates(core),
      Y = Y,
      Z = Z,
      family = "poisson",
      control = list(cluster=cl)
    ) -> modOribatids[[dwf]]
  }
  
  rm(dwf)
  
  stopCluster(cl)
  
  rm(cl)
  
  save(optOribatids, modOribatids, file="Data/optOribatids.rda")
  
} else load(file="../../Data/optOribatids.rda")

## 
if(FALSE) {
  
  st_intersects(grid$Other, substrate) %>%
    sapply(
      function(x,g) colMeans(g[x,-1L]),
      g = st_drop_geometry(substrate)
    ) %>%
    t -> tmp
  
  grid$Other[["Substrate.Sphagn1"]] <- tmp[,"Substrate.Sphagn1"]
  grid$Other[["Substrate.Sphagn2"]] <- tmp[,"Substrate.Sphagn2"]
  grid$Other[["Substrate.Sphagn3"]] <- tmp[,"Substrate.Sphagn3"]
  grid$Other[["Substrate.Sphagn4"]] <- tmp[,"Substrate.Sphagn4"]
  grid$Other[["Substrate.Litter"]] <- tmp[,"Substrate.Litter"]
  grid$Other[["Substrate.Barepeat"]] <- tmp[,"Substrate.Barepeat"]
  
  rm(tmp)
  
  st_intersects(grid$Other, topo) %>%
    sapply(
      function(x,g)
        if(length(x) == 1L) c(Blanket = 1, Hummock = -1.8)[g[x]] else 0,
      g = st_drop_geometry(topo)$Type
    ) %>%
    as.numeric -> grid$Other[["topo"]]
  
  st_intersects(grid$Other, water) %>%
    sapply(
      function(x) if(!length(x)) 0 else 1
    ) -> grid$Other[["Flooded"]]
  
  st_intersects(grid$Other, shrub) %>%
    sapply(
      function(x) x[1L]
    ) %>%
    st_drop_geometry(shrub)$Type[.] %>%
    predict(dat$shrubPoly,.) -> tmp
  grid$Other[["Shrub^1"]] <- tmp[,1L]
  grid$Other[["Shrub^2"]] <- tmp[,2L]
  
  rm(tmp)
  
  SpeciesPsquare <- list()
  
  dat$Y %>%
    apply(
      MARGIN = 2L,
      FUN = function(x)
        glm(y~1, data=data.frame(y=x), family = poisson())$coef
    ) -> SpeciesPsquare$lambda
  
  dat$Y %>%
    apply(
      MARGIN = 2L,
      FUN = function(x)
        -2*sum(
          dpois(
            x = x,
            lambda = exp(
              glm(y~1, data=data.frame(y=x), family = poisson())$coef
            ),
            log = TRUE
          )
        )
    ) -> SpeciesPsquare$null
  
  dat$Y %>%
    apply(
      MARGIN = 2L,
      FUN = function(x)
        -2*sum(
          dpois(
            x = x,
            lambda = x,
            log = TRUE
          )
        )
    ) -> SpeciesPsquare$max
  
  modOribatids %>%
    lapply(
      function(x, Y)
        -2*colSums(
          dpois(
            x = Y,
            lambda = x$fit,
            log=TRUE
          )
        ),
      Y = dat$Y
    ) -> SpeciesPsquare$fit
  
  SpeciesPsquare$fit %>%
    lapply(
      function(x, y)
        1 - (y$max - x)/(y$max - y$null),
      y = SpeciesPsquare
    ) -> SpeciesPsquare$R2
  
  SpeciesPsquare$R2 %>%
    lapply(
      function(x)
        -log10(1 - x)
    ) -> SpeciesPsquare$Q
  
  lmQL <- lm(Q.power~lambda, data = as.data.frame(SpeciesPsquare[c(1L:3L,6L)]))
  
  ## summary(lmQL)
  
  c(0.157,1,10,35.26) %>%
    log %>%
    data.frame(lambda = .) %>%
    predict(lmQL, newdata = .) %>%
    {1 - 10^-.}
  
  save(SpeciesPsquare, lmQL, file="../../Data/SpeciesPsquare.rda")
  
} else load(file="../../Data/SpeciesPsquare.rda")

if(FALSE) {
  rng <- c(0,1000)
  
  colmap <- rainbow(1200L)[1L:1000L]
  
  #' dwf = "linear"
  for(dwf in names(modOribatids)) {
    
    ## X11(width = 16, height = 24)
    png(
      sprintf("../../Image/Oribatid_predictions-%s.png",dwf),
      width = 800, height = 1600
    )
    
    par(fig = c(0,0.125,0,1), mar=c(2,5.3,3,0.5))
    
    seq(rng[1L], rng[2L], length.out=256L) %>%
      matrix(1L,length(.)) %>%
      image(
        zlim = rng,
        col = colmap,
        axes = FALSE,
        ylab = expression(paste("Abundance ",(ind~core^{-1}))),
        cex.lab = 2
      )
    
    box()
    
    c(0,1,3,10,30,100,300,1000) %>%
      {axis(
        side = 2L,
        label = .,
        at = log1p(.)/log1p(1000),
        cex.axis = 2
      )}
    
    par(mar=c(0.5,0.5,2.5,0.5))
    
    osx <- 1/8
    osy <- 1/5
    figx <- osx
    figy <- 1
    
    #' sp=rownames(modOribatids[[dwf]]$sp)[1L]
    for(sp in rownames(modOribatids[[dwf]]$sp)) {
      par(fig = c(figx,figx + osx,figy - osy,figy), new=TRUE)
      prdOribatids(sp = sp, dwf = dwf) %>%
        {rasterFromXYZ(
          cbind(
            st_coordinates(.),
            .[[sp]]
          ),
          crs = st_crs(.)
        )} %>%
        as.matrix %>%
        t %>%
        .[,ncol(.):1L] %>%
        log1p %>%
        image(
          x = seq(0,2.6,length.out=nrow(.)),
          y = seq(0,10,length.out=ncol(.)),
          zlim = log1p(rng),
          col = colmap,
          main = sprintf(
            "%s\n(%0.3f)",
            strsplit(sp,".",fixed=TRUE)[[1L]][2L],
            SpeciesPsquare$R2[[dwf]][sp]
          ),
          axes = FALSE,
          asp = 1,
          xlab = "",
          ylab = ""
        )
      points(
        st_coordinates(core),
        pch = ifelse(core[[sp]], 21L, 4L),
        bg = colmap[1L + round((length(colmap) - 1)*log1p(core[[sp]])/log1p(1000))]
      )
      figx <- figx + osx
      if(figx >= 1) {
        figx <- osx
        figy <- figy - osy
      }
    }
    
    dev.off()
    
  }
  
  rm(rng,colmap,dwf,sp,figx,figy,osx,osy)
  
}

## Get the model's coefficients
if(FALSE) {
  
  modOribatids$power$model %>%
    coef(s="lambda.min") %>%
    as.numeric %>%
    matrix(ncol = 35L) -> coefMat
  rownames(modOribatids$power$sp) -> colnames(coefMat)
  
  c(
    "cte",
    colnames(dat$Z),
    names(modOribatids$power$par$wh)
  ) -> rownames(coefMat)
  
  coefMat %>%
    apply(1L, function(x) sum(x != 0)) %>%
    {100*./ncol(coefMat)}
  
  coefMat %>%
    apply(2L, function(x) sum(x != 0)) %>%
    {100*./nrow(coefMat)}
  
  saveRDS(coefMat, file="../../Data/coefMat.rds")
  
} else coefMat <- readRDS(file="../../Data/coefMat.rds")

if(FALSE) {
  cbind(
    dat$Y %>% apply(2L, min),
    dat$Y %>%
      apply(
        MARGIN = 2L,
        FUN = function(x)
          glm(y~1, data=data.frame(y=x), family = poisson())$coef
      ) %>%
      exp,
    dat$Y %>% apply(2L, max),
    modOribatids$power$fit %>% apply(2L, min),
    modOribatids$power$fit %>% apply(2L, function(x) exp(mean(log(x)))),
    modOribatids$power$fit %>% apply(2L, max),
    colnames(dat$Y) %>%
      sapply(
        function(x) {
          tmp <- prdOribatids(x,"power")[[x]]
          c(min(tmp, na.rm=TRUE),
            exp(mean(log(tmp), na.rm=TRUE)),
            max(tmp, na.rm=TRUE))
        }
      ) %>% t
  ) -> rangeTable
  
  paste(
    rep(c("Observed","Model","Grid"),each=3L),
    rep(c("min","mean","max"),3L),
    sep="_"
  ) -> colnames(rangeTable)
  
  rangeTable %<>% cbind(R2 = SpeciesPsquare$R2$power)
  
  rownames(rangeTable) %<>%
    strsplit(".",fixed=TRUE) %>%
    lapply(tail, n=1L) %>%
    unlist
  
  data.frame(
    R2 = sprintf("%0.3f",rangeTable[,"R2"]),
    Observed = rangeTable %>%
      {sprintf("%0.2f [%0.2f, %0.2f]",.[,2L],.[,1L],.[,3L])},
    Model = rangeTable %>%
      {sprintf("%0.2f [%0.2f, %0.2f]",.[,5L],.[,4L],.[,6L])},
    Grid = rangeTable %>%
      {sprintf("%0.2f [%0.2f, %0.2f]",.[,8L],.[,7L],.[,9L])},
    row.names = rownames(rangeTable)
  ) -> rngTab
  
  saveRDS(rngTab, file="../../Data/rngTab.rds")
  
} else rngTab <- readRDS(file="../../Data/rngTab.rds")

## Beta diversity
if(FALSE) {
  
  BD <- list()
  
  BD$n <- nrow(dat$Y)
  
  dat$Y %>%
    {(. / rowSums(.))^0.5} -> BD$Y.tr
  
  BD$cmY <- colMeans(BD$Y.tr)
  
  BD$Y.tr %>%
    scale(center=BD$cmY, scale=FALSE) %>%
    {.^2} -> BD$S
  
  BD$SS_total <- sum(BD$S)
  
  BD$BD_total <- BD$SS_total/(BD$n - 1L)
  
  BD$var.i <- colSums(BD$S)/(BD$n - 1L)
  BD$var.j <- rowSums(BD$S)/(BD$n - 1L)
  
  BD$SCBD <- colSums(BD$S)/BD$SS_total
  BD$LCBD <- rowSums(BD$S)/BD$SS_total
  
  ## BD$BD_total
  ## quantile(BD$SCBD, c(0,0.5,1), na.rm=TRUE)
  ## quantile(BD$LCBD, c(0,0.5,1), na.rm=TRUE)
  
  BD$fit <- list()
  
  grid$LCBD <- grid$Other[,"geometry",drop=FALSE]
  
  #' dwf="linear"
  for(dwf in names(modOribatids)) {
    
    BD$fit[[dwf]] <- list()
    
    modOribatids[[dwf]]$fit %>%
      {(./rowSums(.))^0.5} -> tmp
    
    BD$fit[[dwf]]$cmY <- colMeans(tmp)
    
    tmp %>%
      scale(center=BD$fit[[dwf]]$cmY, scale=FALSE) %>%
      {.^2} -> BD$fit[[dwf]]$S
    
    BD$fit[[dwf]]$S %>%
      {rowSums(.)/sum(.)} -> BD$fit[[dwf]]$LCBD
    
    BD$fit[[dwf]]$S %>%
      {colSums(.)/sum(.)} -> BD$fit[[dwf]]$SCBD
    ## quantile(BD$fit[[dwf]]$LCBD, c(0,0.5,1), na.rm=TRUE)
    ## quantile(BD$fit[[dwf]]$SCBD, c(0,0.5,1), na.rm=TRUE)
    
    tmp <- NULL
    
    #' sp=rownames(modOribatids[[dwf]]$sp)[1L]
    for(sp in rownames(modOribatids[[dwf]]$sp))
      tmp %<>% cbind(prdOribatids(sp = sp, dwf = dwf)[[sp]])
    
    tmp %>%
      {(./rowSums(.))^0.5} %>%
      scale(center=BD$fit[[dwf]]$cmY, scale=FALSE) %>%
      {.^2} %>%
      {rowSums(.)/sum(BD$fit[[dwf]]$S)} -> grid$LCBD[[dwf]]
    
    ## quantile(grid$LCBD[[dwf]], c(0,0.5,1), na.rm=TRUE)
    
  }
  
  rm(dwf,tmp,sp)
  
}

if(FALSE) {
  
  png("../../Image/Oribatid_LCBD.png", width = 320, height = 320)
  par(mar=c(4.1,4.1,0.6,0.6))
  
  data.frame(
    fit = log10(BD$fit$power$LCBD),
    obs = log10(BD$LCBD)
  ) -> tmp
  
  plot(
    fit~obs,
    data = tmp,
    xlim = log10(c(0.003,0.05)),
    ylim = log10(c(0.003,0.05)),
    xlab = "LCBD from observed species abundances",
    ylab = "LCBD from fitted species abundances",
    axes = FALSE
  )
  
  c(0.002,0.005,0.01,0.02,0.05) %>%
    axis(1L, at=log10(.), label=.)
  
  c(0.002,0.005,0.01,0.02,0.05) %>%
    axis(2L, at=log10(.), label=.)
  
  tmp <- lm(fit~obs, data = tmp)
  
  abline(tmp)
  
  xx <- seq(log10(0.002), log10(0.05), length.out=200L)
  
  tmp %>%
    predict(
      newdata = data.frame(obs=xx),
      interval = "confidence"
    ) -> int
  lines(x = xx, y = int[,"lwr"], lty = 2L)
  lines(x = xx, y = int[,"upr"], lty = 2L)
  
  tmp %>%
    predict(
      newdata = data.frame(obs=xx),
      interval = "prediction"
    ) -> int
  lines(x = xx, y = int[,"lwr"], lty = 3L)
  lines(x = xx, y = int[,"upr"], lty = 3L)
  
  ss <- summary(tmp)
  tt <- (ss$coefficients[2L,1L] - 1)/ss$coefficients[2L,2L]
  ss$coefficients[1L,4L]
  
  text(
    labels = parse(
      text = sprintf(
        "italic(R)[italic(adj)]^2 == %0.3f",
        ss$adj.r.squared
      ),
    ),
    x = log10(0.003),
    y = log10(0.04),
    adj = 0
  )
  
  text(
    labels = parse(
      text = sprintf(
        "atop(italic(y) == alpha*italic(x)^beta,atop(italic(Pr)[alpha == 1] == %0.3f,italic(Pr)[beta == 1] == %0.3f))",
        ss$coefficients[1L,4L],
        2*pt(q = tt, df = ss$df[2L], lower.tail = FALSE)
      ),
    ),
    x = log10(0.003),
    y = log10(0.025),
    adj = 0
  )
  
  dev.off()
  
  rm(tmp, xx, ss, tt, int)
  
}

BD$fit$power$LCBD %>%
  quantile(na.rm=TRUE, prob=c(0,0.5,1))

BD$fit$power$SCBD %>%
  quantile(na.rm=TRUE, prob=c(0,0.5,1))

grid$LCBD %>%
  st_drop_geometry %>%
  apply(2L, quantile, na.rm=TRUE, prob=c(0,0.5,1)) %>% t

if(FALSE) {
  c(linear="Linear",power="Power",hyperbolic="Hyperbolic",
    spherical="Spherical",exponential="Exponential",
    Gaussian="Gaussian",hole_effect="Hole effect") -> labSwap
  
  saveRDS(labSwap, file="../../Data/labSwap.rds")
  saveRDS(grid, file="../../Data/grid.rds")
} else {
  labSwap <- readRDS(file="../../Data/labSwap.rds")
  grid <- readRDS(file="../../Data/grid.rds")
}

if(FALSE) {
  
  rng <- c(0,0.05)
  
  colmap <- rainbow(1200L)[1L:1000L]
  
  ## X11(width = 16.5, height = 6.40)
  png(
    "../../Image/Oribatid_Beta.png",
    width = 825, height = 320
  )
  
  figx <- 0.125
  
  par(fig = c(0,figx,0,1), mar=c(2,5.3,3,0.5))
  
  seq(rng[1L], rng[2L], length.out=256L) %>%
    matrix(1L,length(.)) %>%
    image(
      zlim = rng,
      col = colmap,
      axes = FALSE,
      ylab = expression(paste(italic(LCBD)," index")),
      cex.lab = 1.75
    )
  
  box()
  
  seq(rng[1L],rng[2L],0.01) %>%
    {axis(
      side = 2L,
      label = .,
      at = ./rng[2L],
      cex.axis = 1.25
    )}
  
  #' dev.off()
  
  par(mar=c(2.0,2.0,2,0.5))
  
  #' dwf="linear"
  for(dwf in names(modOribatids)) {
    
    par(fig = c(figx,figx + 0.125,0,1), new=TRUE)
    
    grid$LCBD %>%
      {rasterFromXYZ(
        cbind(
          st_coordinates(.),
          .[[dwf]]
        ),
        crs = st_crs(.)
      )} %>%
      as.matrix %>%
      t %>%
      .[,ncol(.):1L] %>%
      image(
        x = seq(0,2.6,length.out=nrow(.)),
        y = seq(0,10,length.out=ncol(.)),
        zlim = rng,
        col = colmap,
        main = labSwap[dwf],
        axes = FALSE,
        asp = 1,
        xlab = "",
        ylab = ""
      )
    
    points(
      st_coordinates(core),
      pch = 21L,
      bg = getColors(BD$fit[[dwf]]$LCBD, colmap, rng)
    )
    
    axis(1L)
    if(dwf=="linear") axis(2L)
    
    figx <- figx + 0.125
  }
  
  dev.off()
  
  ## For the best model only:
  
  rng <- c(0,0.05)
  
  colmap <- rainbow(1200L)[1L:1000L]
  
  ## X11(width = 6.0, height = 12.8)
  png(
    "../../Image/Oribatid_Beta_power.png",
    width = 300, height = 640
  )
  
  par(fig = c(0,0.35,0,1), mar=c(2,4.8,2,1.5))
  
  seq(rng[1L], rng[2L], length.out=256L) %>%
    matrix(1L,length(.)) %>%
    image(
      zlim = rng,
      col = colmap,
      axes = FALSE,
      ylab = expression(paste(italic(LCBD)," index")),
      cex.lab = 1.5
    )
  
  box()
  
  seq(rng[1L],rng[2L],0.01) %>%
    {axis(
      side = 2L,
      label = .,
      at = ./rng[2L],
      cex.axis = 1.5
    )}
  
  par(fig = c(0.35,1,0,1), mar=c(2,2,1,0.5), new=TRUE)
  
  dwf <- "power"
  
  grid$LCBD %>%
    {rasterFromXYZ(
      cbind(
        st_coordinates(.),
        .[[dwf]]
      ),
      crs = st_crs(.)
    )} %>%
    as.matrix %>%
    t %>%
    .[,ncol(.):1L] %>%
    image(
      x = seq(0,2.6,length.out=nrow(.)),
      y = seq(0,10,length.out=ncol(.)),
      zlim = rng,
      col = colmap,
      axes = FALSE,
      asp = 1,
      xlab = "",
      ylab = ""
    )
  
  points(
    st_coordinates(core),
    pch = 21L,
    ## bg = getColors(BP$LCBD, colmap, rng)
    bg = getColors(BD$fit[[dwf]]$LCBD, colmap, rng)
  )
  
  axis(1L)
  
  axis(2L)
  
  dev.off()
  
  rm(rng,colmap,dwf,figx)
  
}

if(FALSE) {
  
  ## Predictions for substrate density
  
  y <- numeric(nrow(core))
  
  dwf <- names(which.min(optSubsDens$bestval))
  
  ## i=1L
  for(i in unique(FID)) {
    core %>%
      st_coordinates %>%
      .[i != FID,] %>%
      genSEF(
        m = genDistMetric(),
        f = genDWF(
          fun = dwf,
          range = optSubsDens$bestpar[dwf,"range"],
          shape = optSubsDens$bestpar[dwf,"shape"]
        )
      ) -> tmp
    
    tmp %>%
      as.matrix %>%
      glmnet(
        y = core$SubsDens[i != FID],
        alpha = 1,
        lambda = modSubsDens[[dwf]]$model$lambda.min
      ) %>%
      predict(
        newx = predict(
          tmp,
          core %>%
            st_coordinates %>%
            .[i == FID,]
        )
      ) -> y[i == FID]
  }
  
  sqrt(var(core$SubsDens))
  
  mean((core$SubsDens - y)^2) %>% sqrt
  
  Psquare(core$SubsDens, y)
  
  ## Predictions for substrate water content
  
  y <- numeric(nrow(core))
  
  dwf <- names(which.min(optWatrCont$bestval))
  
  ## i=1L
  for(i in unique(FID)) {
    core %>%
      st_coordinates %>%
      .[i != FID,] %>%
      genSEF(
        m = genDistMetric(),
        f = genDWF(
          fun = dwf,
          range = optWatrCont$bestpar[dwf,"range"],
          shape = optWatrCont$bestpar[dwf,"shape"]
        )
      ) -> tmp
    
    tmp %>%
      as.matrix %>%
      glmnet(
        y = core$WatrCont[i != FID],
        alpha = 1,
        lambda = modWatrCont[[dwf]]$model$lambda.min
      ) %>%
      predict(
        newx = predict(
          tmp,
          core %>%
            st_coordinates %>%
            .[i == FID,]
        )
      ) -> y[i == FID]
  }
  
  sqrt(var(core$WatrCont))
  
  mean((core$WatrCont - y)^2) %>% sqrt
  
  Psquare(core$WatrCont, y)
  
  rm(y,dwf,i,tmp)
  
}

rmdwc::rmdcount(c("../../pMEM.Rmd","../../Supplementary/Appendices.Rmd"))
