##
## quick estimates 

tBinomTest <- function(x,trim=0.1) {
  np <- sum(x > trim,na.rm=TRUE)
  nm <- sum(x < (-1*trim),na.rm=TRUE)
  N <- np + nm
  if (N <= 0) {
    p$p.value = NA
    p$statistic = NA
  } else {
    p <- binom.test(np,N,p=0.5);
    p$statistic <- (2*np - N)/sqrt(N)	
  }
  return(p)
}

regmap <- function(m,scale=c(-6,6),na.color=par("bg"),...) {
  hold <- as.matrix(m[nrow(m):1,])
  colnames(hold) <- colnames(m)
  m <- hold
  rm(hold)
  di <- dim(m)
  nr <- di[1]
  nc <- di[2]
  
  if(!is.null(scale)) {
    m[m < scale[1]] <- scale[1]
    m[m > scale[2]] <- scale[2]      
    graph.scale <- rep(NA,length=ncol(m))
    gs <- scale[1]:scale[2]
    graph.scale[1:length(gs)] <- gs
    m <- rbind(graph.scale,m)
  }
  
  r.cex <- 0.2 + 1/log10(nr)
  ##   c.cex <- 0.2 + 1/log10(nc)
  c.cex <- r.cex
  if (c.cex > r.cex) 
    c.cex <- r.cex
  op <- par(no.readonly = TRUE)
  layout(matrix(c(0, 3, 2, 1), 2, 2, byrow = TRUE), widths = c(1,4), heights = c(1, 4), respect = FALSE)
  par(mar = c(5, 0, 0, 5))
  image(1:ncol(m), 1:nrow(m), t(m), axes = FALSE, xlim = c(0.5,ncol(m) + 0.5), ylim = c(0.5, nrow(m) + 0.5), xlab = ",",ylab = ",",...)
  if(!is.na(na.color) & any(is.na(m))) {
    na.m <- ifelse(is.na(m),1,NA)
    ##print(na.m)
    image(1:ncol(m), 1:nrow(m), t(na.m), axes = FALSE,, xlab = ",",ylab = ",",col=na.color,add=T)
  }
  axis(3, 1:ncol(m), las = 2, line = -0.5, tick = 0, labels = colnames(m), cex.axis = c.cex)
  axis(2, 1:nrow(m), las = 2, line = -0.5, tick = 0, labels = rownames(m), cex.axis = r.cex)
  par(op)
}

summarizeByRegion <- function (eset, genome, chrom = "ALL", ref = NULL, center = TRUE, aggrfun = NULL, p.value = 0.005, FUN = t.test, verbose = TRUE, explode = FALSE, ...) 
{
  if (chrom == "ALL") {
    chrom <- names(attr(genome, "chromInfo"))
    if (is.null(chrom) || is.na(chrom)) 
      stop("chromLoc object does not contain any chromInfo\n")
  }
  else if (chrom == "arms") {
    chrom <- genome@chromLocs$armList
  }
  else if (chrom == "bands") {
    chrom <- genome@chromLocs$bandList
  }
  else if (chrom == "mb") {
    chrom <- genome@chromLocs$mbList
  }
  if (class(eset) == "exprSet") {
    exprs <- eset@exprs
    .Deprecated(msg = Biobase:::EXPRSET_DEPR_MSG)
  }
  else if (class(eset) == "ExpressionSet") 
    exprs <- assayData(eset)$exprs
  else exprs <- eset
  if (!is.null(ref)) {
    if (!is.numeric(ref)) 
      stop("column index's required")
  }
  if (!is.null(ref)) {
    if(verbose) cat("Creating ratios...", "\n")
    ref_mean <- apply(exprs[, ref], 1, mean, na.rm = TRUE)
    exprs <- sweep(exprs, 1, ref_mean, "-")
  }
  if (center) 
    exprs <- scale(exprs, scale = F)
  chrom <- as.character(chrom)
  sum.statistic <- matrix(NA,ncol=ncol(exprs),nrow=length(chrom))
  rownames(sum.statistic) <- chrom
  colnames(sum.statistic) <- colnames(exprs)
  for (i in chrom) {
    if(verbose) 
      cat("Testing ",i,"\n")
    ag.list <- .usedChromExprs(exprs, genome, i, aggrfun)
    if (is.null(ag.list))
      next
    region.gx <- ag.list$exprs
    stat <- try(apply(region.gx, 2, FUN, ...), silent = !verbose)
    if (inherits(stat, "try-error"))
      next
    if (is.list(stat)) {
      ps <- unlist(lapply(stat, function(x) x$p.value))
      stat <- unlist(lapply(stat, function(x) x$statistic))
      if (!is.na(p.value)) {
        stat[ps > p.value] <- NA
      }
    }
    sum.statistic[i,] <- stat
  }
  if (explode) {
    nExprs <- exprs
    nExprs[1:nrow(nExprs), 1:ncol(nExprs)] <- NA
    if (verbose)
      cat("Exploding summary matrix...", "\n")
    for (i in chrom) for (j in 1:ncol(exprs)) nExprs[.usedChromExprs(exprs, 
                                                                     genome, i)$geneIDs, j] <- sum.statistic[i, j]
    return(nExprs)
  }
  return(sum.statistic)
}

cgma <- summarizeByRegion

##
## This returns a MATRIX or a VECTOR, depending on the summarize option
##
## smoothing code

movt <- function(v,span=NULL,summarize=mean){

   if (is.null(span)) {
        spanErr <- try(span <- seq(25, length(v) * 0.3, by = 5),
            silent = T)
        if (inherits(spanErr, "try-error"))
            return(NULL)
    }
    if (any(span < 1)) {
        span <- floor(length(v) * span)
    }
    if (length(v) <= 3)
        return(NULL)
    if (max(span) > length(v)) {
        warning("span is longer then data series. It will be truncated")
        span[span > length(v)] <- length(v)
    }
    if (is.null(span))
        stop("Invalid window span")


   finalMatrix <- matrix(data=NA,nrow=length(span),ncol=length(v))

   for(ABC in 1:length(span)){
		meanMatrix <- .genMeanMatrix(v,span[ABC],length(v))
		summary <- apply(meanMatrix,2,mean,na.rm=T)
		finalMatrix[ABC,] <- summary
	}

    if (!is.null(summarize)) finalMatrix <- apply(finalMatrix, 2, summarize)
    return(finalMatrix)
}



.genMeanMatrix <- function(test,curWindowSize,sampleSize) {
	meanMatrix <- matrix(data=NA,ncol=sampleSize,nrow=sampleSize)
	start <- 1
	end <- curWindowSize
	meanCounter <- 1
	while((meanCounter + curWindowSize - 1) <= sampleSize) {
	   vec1 <- test[c(start:end)]
	   m <- try(t.test(vec1))
	   if (inherits(m,"try-error")) {
	      m <- NA
	   } else {
	      m <- m$statistic
	   }
	   meanMatrix[meanCounter,c(meanCounter:(meanCounter+curWindowSize-1))] <- m
	   start <- (start + 1)
	   end <- (end + 1)
	   meanCounter <- (meanCounter + 1)
	}
	return(meanMatrix)
}




movbin  <- function(v,span=NULL,summarize=mean) {
  if(is.null(span)){
	spanErr <- try(span <- seq(25,length(v)*.3,by=5),silent=T)
	if(inherits(spanErr,"try-error")) return(NULL)
  }

  if(any(span < 1)) {
    span <- floor(length(v)*span)
  }
  if(length(v) <=3) return(NULL)
  
  if(max(span) > length(v)) {
    warning("span is longer then data series. It will be truncated")
    span[span>length(v)] <- length(v)
  }	  
  
  if(is.null(span)) stop("Invalid window span")
  
  v[is.na(v)] <- -999
  
  temp <- .C("mspan_mov_binom_mx",
             as.double(v),
             as.integer(length(v)),
             as.integer(span),
             as.integer(length(span)),
             mat = double(length(v)*length(span)),PACKAGE="reb")
  temp <- matrix(temp$mat,nrow=length(span))
  temp[temp == 0] <- NA
  temp[temp == 999] <- 0
  
  if(!is.null(summarize)) temp <- apply(t(temp),1,summarize)
  
  return(temp)
}


rmAmbigMappings <- function(cL) {
  if(class(cL) != "chromLocation")
    stop("this function acts on a ",sQuote("chromLocation"),"object")
  package <- paste(cL@dataSource,"CHRLOC",sep="")
  ids <- ls(env=get(package))
  locs <- mget(ids,get(package))
  chr <- unlist(lapply(locs,function(x) length(unique(names(x)))))
  ambigN <- names(chr)[chr > 1]
  chrLocs <- cL@chromLocs
  for(i in 1:length(chrLocs)) {
    ix <- which(names(chrLocs[[i]]) %in% ambigN)
    if(length(ix) > 0) {
      chrLocs[[i]] <- chrLocs[[i]][-ix]
    }
  }
  cL@chromLocs <- chrLocs
  return(cL)
}

absMax <- function(x) {
  if(all(is.na(x))) return(NA)
  a <- max(x,na.rm=TRUE)
  b <- min(x,na.rm=TRUE)
  if(abs(a) > abs(b)){
    return(a)
  } else {
    return(b)
  }
}

naMean <- function(x) {
  mean(x,na.rm=T)
}

smoothByRegion <- function (eset, genome, chrom = "ALL", ref = NULL, center = FALSE, 
    aggrfun = absMax, method = c("movbin", "supsmu", "lowess","movt"),...) 
{
    
	if (chrom == "ALL") {
	    chrom <- names(attr(genome, "chromInfo"))
	    if (is.null(chrom) || is.na(chrom)) 
		stop("chromLoc object does not contain any chromInfo\n")
	}else if (chrom == "arms") {
		chrom <- genome@chromLocs$armList
	}else if (chrom == "bands") {
		chrom <- genome@chromLocs$bandList
	}else if (chrom == "mb") {
		chrom <- genome@chromLocs$mbList
	}
    
	if(class(eset) == "exprSet"){
		exprs <- eset@exprs
		.Deprecated(msg=Biobase:::EXPRSET_DEPR_MSG)
	}
	else if(class(eset) == "ExpressionSet") exprs <- assayData(eset)$exprs
	else exprs <- eset
		
	
    if (!is.null(ref)) {
        if (!is.numeric(ref)) 
            stop("column index's required")
    }
    method <- match.arg(method)
    if (!exists(method)) 
        stop(sQuote(method), " is not found")
    if (!is.null(ref)) {
        cat("Creating ratios...", "\n")
        ref_mean <- apply(exprs[, ref], 1, mean, na.rm = TRUE)
        exprs <- sweep(exprs, 1, ref_mean, "-")
    }
    if (center) 
        exprs <- scale(exprs, scale = F)
    #temp.eset <- eset
    temp.exprs <- matrix(ncol = ncol(exprs))
    for (chr in chrom) {
        uCG <- try(.usedChromExprs(exprs, genome, chr, aggrfun))
        if (inherits(uCG, "try-error") || is.null(uCG)) 
            next
        if (ncol(uCG$exprs) != ncol(exprs)) 
            next
        locs <- uCG$locs
        names(locs) <- uCG$simpleIDs
        gx <- uCG$exprs
        dix <- duplicated(locs)
        if (sum(dix) & is.null(aggrfun) & method != "movbin") {
            warning(sum(dix), " duplicate locations found...removing duplicates\n")
            gx <- gx[!dix, ]
            locs <- locs[!dix]
        }
        r.matrix <- matrix(NA, ncol = ncol(gx), nrow(gx))
        colnames(r.matrix) <- colnames(gx)
        rownames(r.matrix) <- names(locs)
        for (i in c(1:ncol(gx))) {
            nas <- is.finite(gx[, i])
            x <- locs[nas]
            y <- gx[nas, i]
            sm <- try(switch(method, movbin = try(movbin(y, ...)), 
                supsmu = try(supsmu(x, y, ...)$y), lowess = try(lowess(x, 
                  y, ...)$y),movt = try(movt(y,...))), silent = FALSE)
            if (!inherits(sm, "try-error")) {
                aa <- try(approx(x, sm, xout = locs), silent = TRUE)
                if (!inherits(aa, "try-error")) {
                  r.matrix[, i] <- aa$y
                }
            }
        }
        temp.exprs <- rbind(temp.exprs, r.matrix)
    }
    temp.exprs <- temp.exprs[-1, ]
    return(temp.exprs)
}

reb <- smoothByRegion

.one2many <- function(locs) {
  ul <- vector()
  n <- names(locs)
  for (i in 1:length(n)) {
    temp <- locs[[i]]
    names(temp) <- rep(n[i],length(temp))
    ul <- c(ul,temp)
  }
  return(ul)
}

.estimateCHRLEN <- function(x) {
  if(is.null(x)) return(0)
  return(max(x,na.rm=TRUE) - min(x,na.rm=TRUE))  # Not so good an estimate
}

.MAP2chromLoc <- function(mapEnv,chrEnv,regions) {
  map <- contents(mapEnv)
  map <- unlist(map)
  r <- list()
  for(i in regions) {
    pattern <- paste("^",i,sep="")
    ix <- grep(pattern,map,perl=TRUE)
    subnames <- names(map[ix])
    locs <- try(mget(subnames,env=chrEnv,ifnotfound=NA))
    if (inherits(locs,"try-error") || length(locs) == 0) {
      r <- c(r,list(NULL))
    } else {
      unlisted.locs <- .one2many(locs)
      r <- c(r,list(unlisted.locs))
    }
  }
  names(r) <- regions
  return(r)
}

buildChromMap <- function(dataPkg,regions) {
  if (!require(dataPkg, character.only=TRUE))
    stop(paste("Package:",dataPkg,"is not available"))
  
  pEnv <- paste("package",dataPkg,sep=":")
  
  chrLocList <- .MAP2chromLoc(get(paste(dataPkg,"MAP",sep=""),pos=pEnv), get(paste(dataPkg,"CHRLOC",sep=""),pos=pEnv),regions)
  
  ## !!! Need to get the version info for dataSource
  newCC <- new("chromLocation",
               organism=get(paste(dataPkg,"ORGANISM",sep=""),pos=pEnv),
               dataSource=dataPkg,
               chromLocs=chrLocList,
               chromInfo=unlist(sapply(chrLocList,.estimateCHRLEN)),
               geneSymbols=get(paste(dataPkg,"SYMBOL",sep=""),pos=pEnv))
  return(newCC)
}



#####   gff and acgh conversion functions 


### to BAND format functions

isAbnormal <- function(x,percent=0.5){
	x[is.na(x)] <- 0
	if(sum(x>0) > percent*length(x))return (1)
	if(sum(x<0) > percent*length(x))return (-1)
	return(0)
}


cset2band <- function(exprs,genome,chr="ALL",organism=NULL,FUN=isAbnormal,...){
  
  if(is.null(organism)) {
    organism <- tolower(substr(genome@organism,1,1))
  }
  if(class(exprs) == "exprSet") {
    exprs <- exprs@exprs
    .Deprecated(msg=Biobase:::EXPRSET_DEPR_MSG)
  }
  cytoEnv <- NULL
  cytoEnv <- switch(organism,
                    "h"=get("Hs.cytoband","package:idiogram"),
                    "r"=get("Rn.cytoband","package:idiogram"),
                    "m"=get("Mm.cytoband","package:idiogram"),
                    NULL)
  if(is.null(cytoEnv))
    stop("Cannot determine organism type, please specify (h)uman, (r)at, or (m)ouse")

  	if (chr == "ALL") {
	    chr <- names(attr(genome, "chromInfo"))
	    if (is.null(chr) || is.na(chr)) 
		stop("chrLoc object does not contain any chrInfo\n")
	}else if (chr == "arms") {
		chr <- genome@chromLocs$armList
	}else if (chr == "bands") {
		chr <- genome@chromLocs$bandList
	}else if (chr == "mb") {
		chr <- genome@chromLocs$mbList
	}

	returnMat <- NULL
	
	for(chr in chr){  
	  bands <- paste(chr,(gsub("\\..*", "",attr(get(chr,cytoEnv),"band"))),sep="")
	  start <- as.numeric(attr(get(chr,cytoEnv),"start"),sep="")
	  end <- as.numeric(attr(get(chr,cytoEnv),"end"),sep="")
	  actualEnd <- end[length(end)] 
	  start <- start[!duplicated(bands)]
	  end <- end[!duplicated(bands)]
	  bands <- bands[!duplicated(bands)]	
	  len <- length(start)
	  
	  end[1:len-1] <- start[2:len]
	  end[len] <- actualEnd
	  abc <- .usedChromExprs(exprs,genome,chr)
	  
	  ids <- abc$geneIDs
	  names <- vector(length=length(ids))
	
	  for(i in 1:length(names)){
		counter <- 1
		cur <- abc$locs[i]
		try(while(cur < start[counter] | cur > end[counter]) counter <- counter + 1,silent=T)
		names[i] <- bands[counter]
	  }
	
	  bands <- names
	  names(bands) <- ids
	
	  aggregated <- aggregate(abc$exprs,list(bands),FUN,...)
	  tempMat <- as.matrix(aggregated[2:length(aggregated)])
	  rownames(tempMat) <- as.character(aggregated[,1])
	  returnMat <- rbind(returnMat,tempMat)
	  }
return(returnMat)	
}

####     gff3     #####

writeGFF3 <- function(cset,genome,chr,file.prefix="temp.gff",organism=NULL){
  
  if(is.null(organism)) {
    organism <- tolower(substr(genome@organism,1,1))
  }
  cytoEnv <- NULL
	cytoEnv <- switch(organism,
                          "h"=get("Hs.cytoband","package:idiogram"),
                          "r"=get("Rn.cytoband","package:idiogram"),
                          "m"=get("Mm.cytoband","package:idiogram"),
                          NULL)
  if(is.null(cytoEnv))
    stop("Cannot determine organism type, please specify (h)uman, (r)at, or (m)ouse")
  
  tempMat <- cset2band(cset,genome,chr)
  
  gffList <- list()
  counter <- 0
  rows <- match(paste(chr,unique(gsub("\\..*", "",attr(get(chr,cytoEnv),"band"))),sep=""),rownames(tempMat))
  locs <- attr(get(chr,cytoEnv),"start")
  
  for(i in 1:ncol(tempMat)){
    sID <- colnames(tempMat)[i]
    ##for(j in c(1:22,"X","Y")){
    temp <- tempMat[rows,i]
    inAb <- FALSE
    for(k in 1:length(temp)){
      ## start aberration
      if(((temp[k]!=0)|is.na(temp[k])) & !inAb) {
        inAb <- TRUE
        start <- locs[k]
      }
      
### in aberration
      if(inAb) {
        
### continue aberration
        cont <- FALSE
        try(if(temp[k] ==  temp[k+1]) cont <- TRUE,silent=TRUE)
        if(cont) next()
        
### end aberration
        a <- try(end <- locs[k+1],silent=TRUE)
        if(inherits(a, "try-error")) end <- locs[k]
        if(is.na(end)) end <- locs[k-1]
        score <- temp[k]
        atributeA <- paste("chr=",chr,sep="")
        atributeB <- paste("id=",sID,sep="")
        attributes <- paste(atributeA,atributeB,"scoreType=FCrelativeToMock; sticky=false; coeff=-0.1; expName=Twist; coeffType=pValue;",sep=";")
        seqID <- paste(sID,counter,sep="_")
        counter <- counter+1
        gffList[[counter]] <- as.vector(c(seqID,".","genomic_DNA",start,end,score,".",".",attributes))
        inAb <- FALSE
      }
    }	
    ##}
  }
  
  toPlot <- gffList
  m <- t(data.frame(toPlot))  ## this converts it to a matrix using transpose
  colnames(m) <- c("SEQ_ID","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE","ATTRIBUTES")
  if(!is.null(file.prefix)) write.table(m,file=file.prefix,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
  else	return(m)
  
  invisible(gffList)
  
}

### CGH strings ###

revish <- function(cset,genome,chr,organism=NULL){
  
  if(is.null(organism)) {
    organism <- tolower(substr(genome@organism,1,1))
  }
  cytoEnv <- NULL
  cytoEnv <- switch(organism,
                    "h"=get("Hs.cytoband","package:idiogram"),
                    "r"=get("Rn.cytoband","package:idiogram"),
                    "m"=get("Mm.cytoband","package:idiogram"),
                    NULL)
  if(is.null(cytoEnv))
    stop("Cannot determine organism type, please specify (h)uman, (r)at, or (m)ouse")
  
  tempMat <- cset2band(cset,genome,chr)
  
  rows <- match(paste(chr,unique(gsub("\\..*", "",attr(get(chr,cytoEnv),"band"))),sep=""),rownames(tempMat))
  rows <- rows[!is.na(rows)]
  
  enhList <- list()
  dimList <- list()
  for(i in 1:ncol(tempMat)){
    
    inEnh <- FALSE
    inDim <- FALSE
    enhString <- ""
    dimString <- ""
    
    ##for(j in c(1:22,"X","Y")){
    temp <- tempMat[rows,i]
    ids <- names(temp)
    for(q in 1:length(ids)){
      a <- nchar(ids[q])
      ids[q] <- substr(ids[q],a-2,a)	
    }		
    seperator <- paste(",",chr,sep="")
    
    for(k in 1:length(temp)){
      if(inEnh | inDim){
        try(if(inEnh & (temp[k+1] >0)) next(),silent=T)
        try(if(inDim & (temp[k+1] <0)) next(),silent=T)
        if(inEnh){
          inEnh <- FALSE
          inDim <- FALSE
          enhString <- paste(enhString,ids[k],sep="")
        }
        if(inDim){
          inDim <- FALSE
          inEnh <- FALSE
          dimString <- paste(dimString,ids[k],sep="")
        }
        next()
      }			
      if(temp[k]==0) {
        inEnh <- FALSE
        inDim <- FALSE
        next()
      }
      if(temp[k]>0) {
        inEnh <- FALSE
        inDim <- FALSE
        try(if(temp[k+1]>0) inEnh <- TRUE,silent=TRUE)
        enhString <- paste(enhString,ids[k],sep=seperator)
        next()
      }
      if(temp[k]<0) {
        inEnh <- FALSE
        inDim <- FALSE
        try(if(temp[k+1]<0) inDim <- TRUE, silent=TRUE)
        dimString <- paste(dimString,ids[k],sep=seperator)
        next()
      }		
    }	
    ##}
    enhList[i] <- substr(enhString,2,nchar(enhString))
    dimList[i] <- substr(dimString,2,nchar(dimString))
    ##cat(".")
  }
  names(enhList) <- colnames(tempMat)
  names(dimList) <- colnames(tempMat)
  return(list(enh=enhList,dim=dimList))	
}



fromRevIsh <- function(enhList,dimList,chr,organism="h"){

  cytoEnv <- switch(organism,
                    "h"=get("Hs.cytoband","package:idiogram"),
                    "r"=get("Rn.cytoband","package:idiogram"),
                    "m"=get("Mm.cytoband","package:idiogram"),
                    NULL)
  if(is.null(cytoEnv))
    stop("Cannot determine organism type, please specify (h)uman, (r)at, or (m)ouse")

  bands <- paste(chr,(gsub("\\..*", "",attr(get(chr,cytoEnv),"band"))),sep="")
  start <- as.numeric(attr(get(chr,cytoEnv),"start"),sep="")

  start <- start[!duplicated(bands)]
  bands <- bands[!duplicated(bands)]	
    names(start) <- bands
  startABC <- start
  
  allgenes <- ls(vai.chr@geneSymbols)
  mat1 <- matrix(data=NA,ncol=1,nrow=length(allgenes))
  rownames(mat1) <- allgenes
    
  genes <- .usedChromExprs(mat1,vai.chr,chr)
  names(genes$locs) <- genes$geneIDs
  genes <- genes$locs

  tempMat2 <- matrix(ncol=length(enhList),nrow=length(start))
  rownames(tempMat2) <- names(start)
  tempMat2[,] <- 0
  for(i in 1:ncol(tempMat2)){
    eList <- unlist(strsplit(unlist(enhList[i]),","))
    dList <- unlist(strsplit(unlist(dimList[i]),","))
    test <- 0
##### enh
    for(j in 1:length(eList)) {
      try(test <- nchar(eList[j]),silent=T)
      if (length(test)>0) {
        if(test %in% c(4,5)){
          tempMat2[eList[j],i] <- 1
        }
        
        if(test %in% c(7,8)){
          temp <- eList[j]
          #chr <- substr(temp,1,nchar(temp)-6)
          start <- substr(temp,nchar(temp)-5,nchar(temp)-3)
          stop <- substr(temp,nchar(temp)-2,nchar(temp))
          
	#temp <- names(mb.chr@chromLocs[[chr]][order(mb.chr@chromLocs[[chr]])])
	#	str(temp)
	#	str(startABC)	
	#	stop()
	temp <- names(startABC)
		
          first <- match(paste(chr,start,sep=""),temp)
          last <- match(paste(chr,stop,sep=""),temp)

          tempMat2[temp[first:last],i] <- 1
          
        }
      }
    }
    test <- 0
##### dim
    for(j in 1:length(dList)) {
      
      try(test <- nchar(dList[j]),silent=T)
      if (length(test)>0) {
        if(test %in% c(4,5)){
          tempMat2[dList[j],i] <- -1
        }
        if(test %in% c(7,8)){
          temp <- dList[j]
          #chr <- substr(temp,1,nchar(temp)-6)
          start <- substr(temp,nchar(temp)-5,nchar(temp)-3)
          stop <- substr(temp,nchar(temp)-2,nchar(temp))
          #temp <- names(mb.chr@chromLocs[[chr]][order(mb.chr@chromLocs[[chr]])])
	  temp <- names(startABC)
          first <- match(paste(chr,start,sep=""),temp)
          last <- match(paste(chr,stop,sep=""),temp)

          tempMat2[temp[first:last],i] <- -1
          
        }
      }
    }
  }
  colnames(tempMat2) <- names(enhList)
  return(tempMat2)
}



buildChromCytoband <- function(organism="h"){
tempLoc <- new("chromLocation")
	if(organism=="h") cytoEnv <- Hs.cytoband
		else if(organism=="r") cytoEnv <- Rn.cytoband
		else if(organism=="m")cytoEnv <- Mm.cytoband
		else stop("organism must be h, m, or r")
	
	chrList <- ls(cytoEnv)
	chrLocs <- list()
		
	chromInfo <- vector(length=length(chrList))
	names(chromInfo) <- chrList
	allNames <- vector()
	for(i in 1:length(chrList))
	{
		names <- paste(chrList[i],(gsub("\\..*", "",attr(get(chrList[i],cytoEnv),"band"))),sep="")
		locs <- as.numeric(attr(get(chrList[i],cytoEnv),"start"),sep="")
		names(locs) <- names
		temp <- locs[!duplicated(names)]
		chrLocs[[i]] <- temp
		chromInfo[i] <- max(as.numeric(attr(get(chrList[i],cytoEnv),"end"),sep=""))
		allNames <- c(allNames,unique(names))
	}
	
	names(chrLocs) <- chrList

	
	tempLoc@organism <- "Human"
	tempLoc@dataSource <- "Hs.cytoband?"
	tempLoc@chromLocs <- chrLocs
	tempLoc@chromInfo <- chromInfo

	
return(tempLoc)
}


