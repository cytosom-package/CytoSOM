heatmap.2ClustMat = function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, clustMat = x,
                              distfun = dist, hclustfun = hclust, dendrogram = c("both",
                                                                                 "row", "column", "none"), reorderfun = function(d, w) reorder(d,
                                                                                                                                               w), symm = FALSE, scale = c("none", "row", "column"),
                              na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks,
                              symbreaks = any(x < 0, na.rm = TRUE) || scale != "none",
                              col = "heat.colors", colsep, rowsep, sepcolor = "white",
                              sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan",
                              na.color = par("bg"), trace = c("column", "row", "both",
                                                              "none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks),
                              linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors,
                              cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
                              labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0,
                                                                                      NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5,
                              colRow = NULL, colCol = NULL, key = TRUE, keysize = 1.5,
                              density.info = c("histogram", "density", "none"), denscol = tracecol,
                              symkey = any(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25,
                              key.title = NULL, key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL,
                              key.ytickfun = NULL, key.par = list(), main = NULL, xlab = NULL,
                              ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL,
                              ...)
{
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && any(duplicated(breaks)))
    stop("breaks may not contian duplicate values")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || any(is.na(Rowv)))
    Rowv <- FALSE
  if (is.null(Colv) || any(is.na(Colv)))
    Colv <- FALSE
  else if (all(Colv == "Rowv"))
    Colv <- Rowv
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) &&
        (dendrogram %in% c("both", "row"))) {
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
      if (dendrogram == "both")
        dendrogram <- "column"
      else dendrogram <- "none"
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) &&
        (dendrogram %in% c("both", "column"))) {
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
      if (dendrogram == "both")
        dendrogram <- "row"
      else dendrogram <- "none"
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
    if (length(rowInd) > nr || any(rowInd < 1 | rowInd >
                                   nr))
      stop("Rowv dendrogram doesn't match size of x")
    if (length(rowInd) < nr)
      nr <- length(rowInd)
  }
  else if (is.integer(Rowv)) {
    distr <- distfun(clustMat)
    hcr <- hclustfun(distr)
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(clustMat, na.rm = na.rm)
    distr <- distfun(clustMat)
    hcr <- hclustfun(distr)
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (!isTRUE(Rowv)) {
    rowInd <- nr:1
    ddr <- as.dendrogram(hclust(dist(diag(nr))))
  }
  else {
    rowInd <- nr:1
    ddr <- as.dendrogram(Rowv)
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
    if (length(colInd) > nc || any(colInd < 1 | colInd >
                                   nc))
      stop("Colv dendrogram doesn't match size of x")
    if (length(colInd) < nc)
      nc <- length(colInd)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    distc <- distfun(if (symm)
      clustMat
      else t(clustMat))
    hcc <- hclustfun(distc)
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(clustMat, na.rm = na.rm)
    distc <- distfun(if (symm)
      clustMat
      else t(clustMat))
    hcc <- hclustfun(distc)
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (!isTRUE(Colv)) {
    colInd <- 1:nc
    ddc <- as.dendrogram(hclust(dist(diag(nc))))
  }
  else {
    colInd <- 1:nc
    ddc <- as.dendrogram(Colv)
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (!is.null(colRow))
    colRow <- colRow[rowInd]
  if (!is.null(colCol))
    colCol <- colCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) <
      1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    if (!missing(ColSideColors)) {
      if (!is.character(ColSideColors) || length(ColSideColors) !=
          nc)
        stop("'ColSideColors' must be a character vector of length ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
                      1)
      lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
      if (!is.character(RowSideColors) || length(RowSideColors) !=
          nr)
        stop("'RowSideColors' must be a character vector of length nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
                                           1), 1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  plot.index <- 1
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    plot.index <- plot.index + 1
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    plot.index <- plot.index + 1
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
          c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
        breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!gtools::invalid(na.color) & any(is.na(x))) {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  if (is.null(srtCol) && is.null(colCol))
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 +
           offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1],
         padj = adjCol[2])
  else {
    if (is.null(srtCol) || is.numeric(srtCol)) {
      if (missing(adjCol) || is.null(adjCol))
        adjCol = c(1, NA)
      if (is.null(srtCol))
        srtCol <- 90
      xpd.orig <- par("xpd")
      par(xpd = NA)
      xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2,
                   tick = 0)
      text(x = xpos, y = par("usr")[3] - (1 + offsetCol) *
             strheight("M"), labels = labCol, adj = adjCol,
           cex = cexCol, srt = srtCol, col = colCol)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtCol ignored.")
  }
  if (is.null(srtRow) && is.null(colRow)) {
    axis(4, iy, labels = labRow, las = 2, line = -0.5 + offsetRow,
         tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2])
  }
  else {
    if (is.null(srtRow) || is.numeric(srtRow)) {
      xpd.orig <- par("xpd")
      par(xpd = NA)
      ypos <- axis(4, iy, labels = rep("", nr), las = 2,
                   line = -0.5, tick = 0)
      text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"),
           y = ypos, labels = labRow, adj = adjRow, cex = cexRow,
           srt = srtRow, col = colRow)
      par(xpd = xpd.orig)
    }
    else warning("Invalid value for srtRow ignored.")
  }
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0,
                              xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) +
                                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
                                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
                                                                                                       1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
                              col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in 1:length(colInd)) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in 1:length(rowInd)) {
      if (!is.null(hline)) {
        abline(h = i - 0.5 + hline.vals, col = linecol,
               lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  plot.index <- plot.index + 1
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    ##flag <- try(plot.dendrogram(ddr, horiz = TRUE, axes = FALSE,
    flag <- try(plot(ddr, horiz = TRUE, axes = FALSE,

                                yaxs = "i", leaflab = "none"))
    if ("try-error" %in% class(flag)) {
      cond <- attr(flag, "condition")
      if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
        stop("Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
    }
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    ##flag <- try(plot.dendrogram(ddc, axes = FALSE, xaxs = "i",
    flag <- try(plot(ddc, axes = FALSE, xaxs = "i",
                                leaflab = "none"))
    if ("try-error" %in% class(flag)) {
      cond <- attr(flag, "condition")
      if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?")
        stop("Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
    }
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    mar <- c(5, 4, 2, 1)
    if (!is.null(key.xlab) && is.na(key.xlab))
      mar[1] <- 2
    if (!is.null(key.ylab) && is.na(key.ylab))
      mar[2] <- 2
    if (!is.null(key.title) && is.na(key.title))
      mar[3] <- 1
    par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
    if (length(key.par) > 0)
      do.call(par, key.par)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min.breaks
      max.raw <- max.breaks
    }
    z <- seq(min.raw, max.raw, by = min(diff(breaks)/100))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    if (is.null(key.xtickfun)) {
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      xargs <- list(at = xv, labels = lv)
    }
    else {
      xargs <- key.xtickfun()
    }
    xargs$side <- 1
    do.call(axis, xargs)
    if (is.null(key.xlab)) {
      if (scale == "row")
        key.xlab <- "Row Z-Score"
      else if (scale == "column")
        key.xlab <- "Column Z-Score"
      else key.xlab <- "Value"
    }
    if (!is.na(key.xlab)) {
      mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5,
            cex = par("cex") * par("cex.lab"))
    }
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE,
                      from = min.scale, to = max.scale)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[!omit]
      dens$y <- dens$y[!omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      if (is.null(key.ytickfun)) {
        yargs <- list(at = pretty(dens$y)/max(dens$y) *
                        0.95, labels = pretty(dens$y))
      }
      else {
        yargs <- key.ytickfun()
      }
      yargs$side <- 2
      do.call(axis, yargs)
      if (is.null(key.title))
        key.title <- "Color Key\nand Density Plot"
      if (!is.na(key.title))
        title(key.title)
      par(cex = 0.5)
      if (is.null(key.ylab))
        key.ylab <- "Density"
      if (!is.na(key.ylab))
        mtext(side = 2, key.ylab, line = par("mgp")[1],
              padj = 0.5, cex = par("cex") * par("cex.lab"))
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      if (is.null(key.ytickfun)) {
        yargs <- list(at = pretty(hy)/max(hy) * 0.95,
                      labels = pretty(hy))
      }
      else {
        yargs <- key.ytickfun()
      }
      yargs$side <- 2
      do.call(axis, yargs)
      if (is.null(key.title))
        key.title <- "Color Key\nand Histogram"
      if (!is.na(key.title))
        title(key.title)
      par(cex = 0.5)
      if (is.null(key.ylab))
        key.ylab <- "Count"
      if (!is.na(key.ylab))
        mtext(side = 2, key.ylab, line = par("mgp")[1],
              padj = 0.5, cex = par("cex") * par("cex.lab"))
    }
    else if (is.null(key.title))
      title("Color Key")
    if (trace %in% c("both", "column")) {
      vline.vals <- scale01(vline, min.raw, max.raw)
      if (!is.null(vline)) {
        abline(v = vline.vals, col = linecol, lty = 2)
      }
    }
    if (trace %in% c("both", "row")) {
      hline.vals <- scale01(hline, min.raw, max.raw)
      if (!is.null(hline)) {
        abline(v = hline.vals, col = linecol, lty = 2)
      }
    }
  }
  else {
    par(mar = c(0, 0, 0, 0))
    plot.new()
  }
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  retval$layout <- list(lmat = lmat, lhei = lhei, lwid = lwid)
  if (!is.null(extrafun))
    extrafun()
  invisible(retval)
}

PlotStarsBigLeg <- function(fsom,
                      markers=fsom$map$colsUsed,
                      view="MST", #"grid","tSNE"
                      colorPalette=grDevices::colorRampPalette(
                        c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                          "yellow", "#FF7F00", "red", "#7F0000")),
                      starBg = "white",
                      backgroundValues = NULL,
                      backgroundColor = function(n){
                        grDevices::rainbow(n, alpha=0.3)},
                      backgroundLim = NULL,
                      backgroundBreaks = NULL,
                      backgroundSize = NULL,
                      thresholds=NULL,
                      legend=TRUE,
                      query=NULL,
                      main="",
                      smallTree=F){
   # Add star chart option to iGraph
    igraph::add.vertex.shape("star", clip=igraph::igraph.shape.noclip, plot=mystarBL,
                    parameters=list(vertex.data=NULL,vertex.cP = colorPalette,
                                    vertex.scale=TRUE, vertex.bg = starBg))

    if(is.null(thresholds)){
        # Use MFIs
        data <- fsom$map$medianValues[, markers,drop=FALSE]
        scale <- TRUE
    } else {
        # scale thresholds same as data
        if(fsom$transform){
            warning("Thresholds should be given in the transformed space")
        }
        if(!is.null(fsom$scaled.center)){
          thresholds <- scale(t(thresholds),
                              center = fsom$scaled.center[markers],
                              scale = fsom$scaled.scale[markers])
        }
        # Use pctgs of cells above threshold as star plot values
        data <-
            t(sapply(seq_len(fsom$map$nNodes), function(i) {
                res = NULL
                for(m in seq_along(markers)){
                    res = c(res,
                            sum(subset(fsom$data,
                               fsom$map$mapping[,1] == i)[,
                                                  markers[m]] > thresholds[m])/
                            sum(fsom$map$mapping[,1] == i))

                }
                res
            }))
        scale <- FALSE
    }

    # Choose layout type
    switch(view,
        MST  = { layout <- fsom$MST$l
                    lty <- 1},
        grid = { layout <- as.matrix(fsom$map$grid)
                    lty <- 0},
        tSNE = { layout <- fsom$MST$l2
                    lty <- 0},
        stop("The view should be MST, grid or tSNE. tSNE will only work
                   if you specified this when building the MST.")
    )

    # Choose background colour
    if (!is.null(backgroundValues)) {
        background <- computeBackgroundColorBL(backgroundValues,backgroundColor,
                                             backgroundLim, backgroundBreaks)
        if (is.null(backgroundSize)) {
          backgroundSize <- fsom$MST$size
          backgroundSize[backgroundSize == 0] <- 3
        }
    } else {
        background <- NULL
    }

    # Save plot window settings and minimize margin
    oldpar <- graphics::par(no.readonly = TRUE)
    graphics::par(mar=c(1,1,1,1))

    # Add legend
    if(legend){
        if(!is.null(backgroundValues)){
            # Make plot with folowing layout
            # 1 3
            # 2 3
          if (smallTree) {graphics::layout(matrix(c(1,1,3,3,2,2), 3, 2, byrow = TRUE),
                                           widths=c(1), heights=c(1,1,3))}
          else {
            graphics::layout(matrix(c(1,1,3,3,2,2), 3, 2, byrow = TRUE),
                    widths=c(1), heights=c(1,3,1))}
        } else {
            graphics::layout(matrix(c(1,2), 1, 2, byrow = TRUE),
                   widths=c(1,2), heights=c(1))
        }

       if(is.null(query)){
            plotStarLegendBL(fsom$prettyColnames[markers],
                            colorPalette(ncol(data)))
        } else {
            plotStarQuery(fsom$prettyColnames[markers],
                            values=query == "high",
                            colorPalette(ncol(data)))
        }

        if(!is.null(backgroundValues)){
            if (smallTree){PlotBackgroundLegendBL(backgroundValues,background,cexLegend=2)}
          else
            {PlotBackgroundLegendBL(backgroundValues,background)}
        }
    }

    # Plot the actual graph
    igraph::plot.igraph(fsom$MST$g,
                        vertex.shape = "star",
                        vertex.label = NA,
                        vertex.size = fsom$MST$size,
                        vertex.data = data,
                        vertex.cP = colorPalette(ncol(data)),
                        vertex.scale = scale,
                        layout = layout,
                        edge.lty = lty,
                        mark.groups = background$groups,
                        mark.col = background$col[background$values],
                        mark.border = background$col[background$values],
                        mark.expand	= backgroundSize,
                        main=main,
                        margin=c(0,0,.3,0)
    )
    # Reset plot window
    graphics::par(oldpar)
    graphics::layout(1)
}
## Internal tool, for BigLegendPlot
computeBackgroundColorBL <- function(backgroundValues,
                                    backgroundColor,
                                    backgroundLim = NULL,
                                    backgroundBreaks = NULL){
    # Choose background colour
    backgroundList <- list()
    backgroundColors <- NULL
    if(!is.null(backgroundValues)){
        if(is.numeric(backgroundValues)){
            backgroundList <- as.list(seq_along(backgroundValues))

            if(class(backgroundColor)=="function" &
               !is.null(backgroundBreaks) &
               length(backgroundBreaks)>1)
            {
                backgroundColor <- backgroundColor(length(backgroundBreaks))
            } else if (class(backgroundColor)=="function"){
                backgroundColor <- backgroundColor(100)
                backgroundBreaks <- length(backgroundColor)
            } else if (is.null(backgroundBreaks)){
                backgroundBreaks <- length(backgroundColor)
            }

            if(length(backgroundLim) > 0){
                ids <- cut(c(backgroundLim,backgroundValues),
                           backgroundBreaks
                )[-c(seq_along(backgroundLim))]
            } else {
                ids <- cut(backgroundValues,
                           backgroundBreaks)
            }
            backgroundValues <- ids
#             backgroundColors <- backgroundColor[ids]
        } else {
            if(! is.factor(backgroundValues)){
                backgroundValues <- as.factor(backgroundValues)
            }

            backgroundList <- as.list(seq_along(backgroundValues))

            if(class(backgroundColor)=="function"){
                backgroundColor <- backgroundColor(
                    length(levels(backgroundValues)))
            }

            if(length(backgroundColor) < length(levels(backgroundValues))){
                stop("You specified less backgroundcolours than groups.")
            }

        }
    }
    backgroundColors <- backgroundColor[backgroundValues]

    list(values=backgroundValues,
         col=backgroundColor,
         groups=backgroundList)
}
## Internal tool, for BL
plotStarLegendBL <- function(labels, colors=grDevices::rainbow(length(labels)),
                            main=""){
    graphics::plot(1, type="n", xlab="", ylab="",
        xlim=c(-10, 10), ylim=c(-3, 3),asp=1,
        bty="n",xaxt="n",yaxt="n",main=main)

    graphics::stars(matrix(c(1:(2*length(labels))),nrow=2),col.segments=colors,
        locations=c(0,0),draw.segments = TRUE,add=TRUE,
        inches=FALSE)
    n <- length(labels)
    angle <- 2*pi / n
    angles <- seq(angle/2,2*pi,by=angle)

    left <- (angles > (pi/2) & angles < (3*pi/2))
    x <- c(2,-2)[left+1]
    y_tmp <- c(seq(-2,2,by= 4/(sum(!left)+1))[-c(1,sum(!left)+2)],
                seq(2,-2,by=-4/(sum(left)+1))[-c(1,sum(left)+2)])
    y <- shiftFunctionBL(y_tmp,max((cummax(y_tmp)<0)*seq_along(y_tmp)))

    for(i in seq_along(labels)){
        graphics::text(x= x[i],
            y= y[i],
            labels=labels[i],
            adj = c(as.numeric(left)[i],0.5),
            cex = 0.5)

        graphics::lines(x=c(x[i]+c(-0.2,0.2)[left[i]+1],
                c(1.5,-1.5)[left[i]+1],
                cos(angles[i])),
            y=c(y[i],
                y[i],
                sin(angles[i])),
            col=colors[i],
            lwd=2)
    }
}

##Internal tool, for Big Legend
shiftFunctionBL <- function(x,n){
    c(x[(n+1):length(x)],x[1:n])
}
##Internal tool, for Big Legend
PlotBackgroundLegendBL <- function(backgroundValues, background,
                                 main="Background",cexLegend=1){
    graphics::plot.new()
    if(is.numeric(backgroundValues)) {
        legendContinuous(background$col,
                         as.numeric(gsub(".*,","",
                                         gsub("].*","",
                                              levels(background$values)))))
    } else {
        relCex=exp(-(length(levels(background$values)))/96)*exp(-max(nchar(levels(background$values)))/92)
        orderIndex = stringr::str_order(levels(background$values),numeric=T)
        graphics::legend("center", legend=levels(background$values)[orderIndex],
               fill=background$col[orderIndex],
               cex=relCex*cexLegend,
               ncol =  ceiling(length(levels(background$values)) / (10*cexLegend)),
               bty="n",
               title=main)
    }
}
##Internal tool, for Big Legend
mystarBL <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
    }
    vertex.size    <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
    }
    data <- params("vertex", "data")
    cP <- params("vertex","cP")
    scale <- params("vertex","scale")
    bg <- params("vertex","bg")
    graphics::symbols(coords[, 1], coords[, 2], circles = vertex.size,
                      inches = FALSE, bg = bg, bty='n', add=TRUE)
    graphics::stars(data, locations = coords, labels = NULL,scale=scale,
            len = vertex.size, col.segments = cP,
            draw.segments = TRUE, mar = c(0, 0, 0, 0), add=TRUE,
            inches=FALSE)

}

## Internal tool: gating subset from FlowSOMworshop

gating_subset_toolBox<- function(flowjo_res, gate){

  if(!is.null(flowjo_res$flowSet)){
    res <- lapply(seq_len(length(flowjo_res$flowSet)),
                  function(i){
                    flowjo_res$flowSet[[i]][flowjo_res$gates[[i]][,gate],]
                  })
    names(res) <- flowCore::sampleNames(flowjo_res$flowSet)
    return(list(
      flowSet = flowCore::flowSet(res),
      gates = lapply(flowjo_res$gates, function(x)x[x[,gate], ])
    ))

  } else {
    return(list(flowFrame = flowjo_res$flowFrame[flowjo_res$gates[,gate], ],
                gates = flowjo_res$gates[flowjo_res$gates[,gate], ]))
  }
}

## Internal tool: new parse_flowjo for CytoML 1.12.0

parse_flowjo_CytoML_v12 <- function (files, wsp_file, group = "All Samples")
{
  wsp <- CytoML::open_flowjo_xml(wsp_file)
  o <- capture.output(gates <- suppressMessages(CytoML::flowjo_to_gatingset(wsp,
                                                                       group)))
  files_in_wsp <- gates@data@origSampleVector 
  counts <- as.numeric(gsub(".*_([0-9]*)$", "\\1", files_in_wsp))
  result <- list()
  for (file in files) {
    print(paste0("Processing ", file))
    file_id <- grep(gsub(".*/", "", file), files_in_wsp,fixed=T)
    if (length(file_id) == 0) {
      stop("File not found. Files available: ", gsub("_[0-9]*$",
                                                   "\n", files_in_wsp))
    }
    gate_names <- flowWorkspace::gs_get_pop_paths(gates[[file_id]],path = "auto")
    gatingMatrix <- matrix(FALSE, nrow = counts[file_id],
      ncol = length(gate_names), dimnames = list(NULL,
      gate_names))
    for (gate in gate_names) {
      gatingMatrix[, gate] <- flowWorkspace::gh_pop_get_indices(gates[[file_id]],gate)
    }
if ((unlist(packageVersion("flowWorkspace"))[1] == 3) && (unlist(packageVersion("flowWorkspace"))[2] <= (32)))
{ff <- flowWorkspace::getData(gates[[file_id]], "root")}
    else {ff <- flowWorkspace::gh_pop_get_data(gates[[file_id]], "root")}

    ff@exprs[, "Time"] <- ff@exprs[, "Time"] * 100
    result[[file]] <- list(flowFrame = ff, gates = gatingMatrix)
  }
  if (length(files) == 1) {
    result <- result[[1]]
  }
  else {
    result <- list(flowSet = flowCore::flowSet(lapply(result,
                function(x) x$flowFrame)), gates = lapply(result,
                     function(x) x$gates))
  }
  return(result)
}
## Internal tool: corrected parse_flowjo

    parse_flowjo_CytoML <- function (files, wsp_file, group = "All Samples", plot = FALSE)
    {
       ##wsp <- flowWorkspace::openWorkspace(wsp_file)
      wsp <- CytoML::openWorkspace(wsp_file)
        o <- capture.output(gates <- suppressMessages(CytoML::parseWorkspace(wsp,
                                                                             group)))
        files_in_wsp <- gates@data@origSampleVector 
        counts <- as.numeric(gsub(".*_([0-9]*)$", "\\1", files_in_wsp))
        result <- list()
        for (file in files) {
            print(paste0("Processing ", file))
            file_id <- grep(gsub(".*/", "", file), files_in_wsp,fixed=T)
            if (length(file_id) == 0) {
                stop("File not found. Files available: ", gsub("_[0-9]*$",
                    "\n", files_in_wsp))
            }
            gate_names <- flowWorkspace::getNodes(gates[[file_id]],
                path = "auto")
            gatingMatrix <- matrix(FALSE, nrow = counts[file_id],
                ncol = length(gate_names), dimnames = list(NULL,
                    gate_names))
            for (gate in gate_names) {
                gatingMatrix[, gate] <- flowWorkspace::getIndiceMat(gates[[file_id]],
                    gate)
            }
            ff <- flowWorkspace::getData(gates[[file_id]], "root")
            if (match("Time",colnames(ff@exprs),nomatch = F)) {
                print("rescale time")
               ff@exprs[, "Time"] <- ff@exprs[, "Time"] * 100
            }
            else {print("Time not found for rescaling")}
            result[[file]] <- list(flowFrame = ff, gates = gatingMatrix)
            if (plot) {
                flowWorkspace::plot(gates[[file_id]])
            }
        }
        if (length(files) == 1) {
            result <- result[[1]]
        }
        else {
            result <- list(flowSet = flowCore::flowSet(lapply(result,
                function(x) x$flowFrame)), gates = lapply(result,
                                                          function(x) x$gates))
        }
        return(result)
    }

## Internal tool: GetClusters ?
##if(!exists("GetClusters",mode="function")) {
    GetClusters <- function(fsom) {
      if (class(fsom) == "list" & !is.null(fsom$FlowSOM)) {
        fsom <- fsom$FlowSOM
      }
      if (class(fsom) != "FlowSOM") {
        stop("fsom should be a FlowSOM object.")
      }
      return(fsom$map$mapping[,1])
    }
##}

## Internal tool: GetMetClusters ?
##if(!exists("GetMetaclusters",mode="function")) {
    GetMetaclusters <- function(fsom, meta = NULL){

      if (class(fsom) == "list" & !is.null(fsom$FlowSOM)) {
        if (is.null(meta) & !is.null(fsom$metaclustering)) {
          meta <- fsom$metaclustering
        }
        fsom <- fsom$FlowSOM
      }
      if (class(fsom) != "FlowSOM"){
        stop("fsom should be a FlowSOM object.")
      }
      if(is.null(meta)){
        stop("No metaclustering found.")
      }

      return(meta[fsom$map$mapping[,1]])
    }
##}

    ## Internal tool: Seems that PlotLabels diseapear...
##if(!exists("PlotLabels",mode="function")) {

    PlotLabels <- function(fsom,
                       labels,
                       view="MST",
                       main=NULL,
                       nodeSize=fsom$MST$size,
                       fontSize = 1,
                       backgroundValues = NULL,
                       backgroundColor = function(n){
                         grDevices::rainbow(n,alpha=0.3)},
                       backgroundLim = NULL,
                       backgroundBreaks = NULL){
        switch(view,
               MST  = { layout <- fsom$MST$l
                   lty <- 1},
               grid = { layout <- as.matrix(fsom$map$grid)
                   lty <- 0},
               tSNE = { layout <- fsom$MST$l2
                   lty <- 0},
               stop("The view should be MST, grid or tSNE. tSNE will only work
                   if you specified this when building the MST.")
               )

                                        # Choose background colour
        if(!is.null(backgroundValues)){
            background <- computeBackgroundColor(backgroundValues,backgroundColor,
                                                 backgroundLim, backgroundBreaks)
        } else {
            background <- NULL
        }

        igraph::plot.igraph(fsom$MST$graph,
                            layout=layout,
                            vertex.size=nodeSize,
                            vertex.label=labels,
                            vertex.label.cex = fontSize,
                            edge.lty=lty,
                            mark.groups=background$groups,
                            mark.col=background$col[background$values],
                            mark.border=background$col[background$values],
                            main=main)

    }
##}


## Internal tool: extract meta-clusters count ratio in percent
get_pctgsMT <- function(fSOM,metacl, meta_names = NULL){
  `%>%` <- magrittr::`%>%`
  cell_ids <- fSOM$metaData
  files <- sapply(seq_len(length(cell_ids)),
                  function(i){
                    rep(gsub(".*/", "", names(cell_ids)[i]),
                        cell_ids[[i]][2] - cell_ids[[i]][1] + 1)
                  }) %>%
    unlist()
  pctgs <- table(files, GetClusters(fSOM)) %>%
    as.matrix() %>%
    apply(1, function(x){x/sum(x)}) %>%
    t()
  pctgs_meta <- table(files, GetMetaclusters(fSOM,meta = metacl)) %>%
    as.matrix() %>%
    apply(1, function(x){x/sum(x)}) %>%
    t()
  if(!is.null(meta_names)) colnames(pctgs_meta) <- meta_names
  return(list("pctgs" = as.matrix(pctgs),
              "pctgs_meta" = as.matrix(pctgs_meta)))
}

## Internal tool: extract absolute count of meta-clusters
get_abstgsMT <- function(fSOM,metacl, meta_names = NULL){
  `%>%` <- magrittr::`%>%`
  cell_ids <- fSOM$metaData
  files <- sapply(seq_len(length(cell_ids)),
                  function(i){
                    rep(gsub(".*/", "", names(cell_ids)[i]),
                        cell_ids[[i]][2] - cell_ids[[i]][1] + 1)
                  }) %>%
      unlist()
  pctgs <- table(files, GetClusters(fSOM)) %>%
      as.matrix() %>% apply(1,as.numeric) %>% t()
  pctgs_meta <- table(files, GetMetaclusters(fSOM,meta = metacl)) %>%
      as.matrix()
  ## %>% apply(1,as.numeric) %>% t()
  if(!is.null(meta_names)) colnames(pctgs_meta) <- meta_names
  return(list("abstgs" = pctgs,
              "abstgs_meta" = pctgs_meta))
}

##Internal tool: return p-value of Tukey test, given metacluster names
TukeyTestSarah = function(fSOMTable, metaClust){TukeyHSD(aov(as.formula(paste(metaClust,"~ Treatment")),data=fSOMTable))$Treatment[,4]}

##Internal tool every kind of boxplots-heatmaps
BoxPlotMetaClustFull <- function(TreeMetaCl,Title,treatmentTable,ControlTreatment,BottomMargin,yLab,Norm=FALSE,Marker="",Robust,ClustHeat,ExportData)
{
    ## Search for the marker
    treatmentTable$Treatment=gsub(" ","_",treatmentTable$Treatment,fixed=T)
    if (length(unique(treatmentTable$Treatment)) < 2) {stop("Single treatment")}
    ControlTreatment = gsub(" ","_",ControlTreatment,fixed=T)
    MarkerIndex = which(TreeMetaCl$fSOMTree$prettyColnames == Marker)
    if(length(MarkerIndex) == 1) {
        print(paste("User marker:",Marker,":",names(TreeMetaCl$fSOMTree$prettyColnames)[MarkerIndex]))
        fSOMnbrs = sapply(unique(TreeMetaCl$metaCl),function(metaClust){
            clusterList=which(TreeMetaCl$metaCl == metaClust)
            metaClustIndices=unlist(sapply(clusterList,function(cluster){which(TreeMetaCl$fSOMTree$map$mapping[,1] == cluster)}))
            sapply(TreeMetaCl$fSOMTree$metaData,function(StartEnd){
                indices = intersect((StartEnd[1]:StartEnd[2]),metaClustIndices)
                return(median(TreeMetaCl$fSOMTree$data[indices,MarkerIndex]))
            })
        })
        colnames(fSOMnbrs)=unique(TreeMetaCl$metaCl)
        row.names(fSOMnbrs)=gsub(".*/","",row.names(fSOMnbrs))
        PlotLab=Marker
    }
    else {
    ## constuct fSOMnbrs, according to Norm (false: percentage, true: normalized)
    if (Norm) {
        abstgs=get_abstgsMT(TreeMetaCl$fSOMTree,TreeMetaCl$metaCl)
        fSOMnbrs<-abstgs$abstgs_meta
        if (is.null(treatmentTable$NormalizationFactor)) {stop("No column NormalizationFactor in annotation table")}
        NormFactors = sapply(row.names(fSOMnbrs),function(fileFCS){treatmentTable$NormalizationFactor[which(treatmentTable$files == fileFCS)]})
        for (index in 1:length(NormFactors)) {fSOMnbrs[index,] = fSOMnbrs[index,]/NormFactors[index] }
        PlotLab=paste("size of ",yLab,sep="")
    }
    else {
        pctgs<-get_pctgsMT(TreeMetaCl$fSOMTree,TreeMetaCl$metaCl)
        fSOMnbrs<-pctgs$pctgs_meta
        fSOMnbrs<-fSOMnbrs*100
        PlotLab=paste("% of ",yLab,sep="")
    }
        ##colnames(fSOMnbrs)=unique(TreeMetaCl$metaCl)
        ##row.names(fSOMnbrs)=gsub(".*/","",row.names(fSOMnbrs))
        }
    treatmentsFSOM=sapply(row.names(fSOMnbrs),function(fileFCS){treatmentTable$Treatment[which(treatmentTable$files == fileFCS)]})
    Treatments=unique(treatmentTable$Treatment)
    if(length(which(Treatments == ControlTreatment)) == 0) {stop(paste("No",ControlTreatment,"in annotation table"))}
    treatmentsFSOM=factor(treatmentsFSOM,levels=c(ControlTreatment,setdiff(Treatments,ControlTreatment))) # set control treatment at first
    if(length(MarkerIndex) == 1) {
      pdf(file=NoSpCharForFile(gsub("/","_",paste(Title,"_BoxPlot",Marker,"Metacl.pdf",sep=""),fixed=T)))
      } else {
    if (Norm) {pdf(file=NoSpCharForFile(paste(Title,"_BoxPlotNormMetacl.pdf",sep="")))}
    else {pdf(file=NoSpCharForFile(paste(Title,"_BoxPlotPercentMetacl.pdf",sep="")))}}
    metaclNumber=length(fSOMnbrs[1,])

    par(mfrow=c(6,6),las=2,mar=c(BottomMargin,3,1+max(sapply(colnames(fSOMnbrs),nchar))%/%24,.5),mgp=c(1.8,.8,0)) ## page have 6x6 boxplots
    fSOMnbrs=fSOMnbrs[,unique(unique(TreeMetaCl$metaCl))]
    mClustNames4Plot = unlist(lapply(as.character(colnames(fSOMnbrs)),function(cName){
      if (nchar(cName) < 27) {return(as.character(cName))}
      else {return(as.character(paste(substring(
        cName,(0:((nchar(cName)%/%24))*24),c(1:((nchar(cName)%/%24))*24,nchar(cName))),
        collapse = "\n")))}
    }))
    cex4Title=exp(-min(27,max(sapply(colnames(fSOMnbrs),nchar)))/40)
    for (metaCl in (1:metaclNumber)){ ## boxplots with no annotations
        plotDf=data.frame(PP=fSOMnbrs[,metaCl],TreatmentFSOM=treatmentsFSOM) ## dataframe for box plot
        boxplot(PP ~ TreatmentFSOM,data=plotDf,main=mClustNames4Plot[metaCl],xlab="",ylab=PlotLab,cex.axis=.5,cex.main=cex4Title,cex.lab=.5)
        beeswarm::beeswarm(PP ~ TreatmentFSOM,data=plotDf,main=paste("mtcl",colnames(fSOMnbrs)[metaCl],sep="_"),add=T,cex=.5,col="red")
    }
    par(mfrow=c(6,6),las=2,mar=c(BottomMargin,3,1+max(sapply(colnames(fSOMnbrs),nchar))%/%24,.5),mgp=c(1.8,.8,0))
    PvalPairwiseTable = sapply((1:metaclNumber),function(metaCl) ## construct pval table of tukey pairwise comparison test, boxplots with p-values annotation
    {
        plotDf=data.frame(PP=fSOMnbrs[,metaCl],TreatmentFSOM=treatmentsFSOM,noData=rep(1,length(treatmentsFSOM)))
        if (Robust) {
          invisible(capture.output(tmpNoData <- dunn.test::dunn.test(plotDf$noData,plotDf$TreatmentFSOM,table=F)))
          pairwisePval=tmpNoData$P
          names(pairwisePval) = gsub(" ","",tmpNoData$comparisons,fixed=T) ## create template
          if (length(unique(plotDf$TreatmentFSOM[which(is.finite(plotDf$PP))])) > 1 ){
          invisible(capture.output(tmp <- dunn.test::dunn.test(plotDf$PP,plotDf$TreatmentFSOM)
                                       ##tryCatch({dunn.test::dunn.test(plotDf$PP,plotDf$TreatmentFSOM)},
                  ##error = function(e){
                  ##  print(str(unique(treatmentsFSOM)))
                  ##  print("haha")
                  ##  data.frame(P=rep(1,length(plotDf$PP)*(length(plotDf$PP)-1)/2),comparisons = combn(unique(treatmentsFSOM),2,function(x){paste(x[1],x[2],sep=" - ")}))
                  ##  })
                  ))
            pairwiseSignPval=tmp$P
            names(pairwiseSignPval) = gsub(" ","",tmp$comparisons,fixed=T)
            for(name in names(pairwiseSignPval)){
              pairwisePval[[name]] = pairwiseSignPval[[name]]
              }}
        } else {
          pairwisePval=TukeyHSD(aov(noData ~ TreatmentFSOM,data=plotDf))$TreatmentFSOM[,4]
          pairwisePval[]=1
          names(pairwisePval)=row.names(TukeyHSD(aov(noData ~ TreatmentFSOM,data=plotDf))$TreatmentFSOM) ##create template
          if (length(unique(plotDf$TreatmentFSOM[which(is.finite(plotDf$PP))])) > 1 ){
          pairwiseSignPval=TukeyHSD(aov(PP ~ TreatmentFSOM,data=plotDf))$TreatmentFSOM[,4]
          names(pairwiseSignPval)=row.names(TukeyHSD(aov(PP ~ TreatmentFSOM,data=plotDf))$TreatmentFSOM)
          for(name in names(pairwiseSignPval)){
            pairwisePval[[name]] = pairwiseSignPval[[name]]
          }
          }

        }
        ListSignif=(sapply(1:length(pairwisePval),function(index){
            if (!is.finite(pairwisePval[index])){return(c())}
            if(pairwisePval[index] < 0.0001){return(c("****",strsplit(names(pairwisePval)[index],split = "-")[[1]]))}
            else if(pairwisePval[index] < 0.001){return(c("***",strsplit(names(pairwisePval)[index],split = "-")[[1]]))}
            else if(pairwisePval[index] < 0.01){return(c("**",strsplit(names(pairwisePval)[index],split = "-")[[1]]))}
            else if(pairwisePval[index] < 0.05){return(c("*",strsplit(names(pairwisePval)[index],split = "-")[[1]]))}
        }))
        ListSignif = ListSignif[which(sapply(ListSignif,length) > 0)]
        ListSignifPosIndex = lapply(ListSignif,function(hit){
            return(c(which(levels(plotDf$TreatmentFSOM) == hit[2]),which(levels(plotDf$TreatmentFSOM) == hit[3])))})
        minTr=min(plotDf$PP,na.rm=T)
        maxTr=max(plotDf$PP,na.rm=T)
        boxplot(PP ~ TreatmentFSOM,
                data=plotDf,main=mClustNames4Plot[metaCl],
                xlab="",
                ylab=PlotLab,
                cex.axis=.5,
                cex.main=cex4Title,
                cex.lab=.5,
                ylim=c(minTr,length(ListSignif)*abs(maxTr-minTr)*.2+maxTr)
                )
        if (length(ListSignif) > 0) {
            if (length(pairwisePval) > 1) ## more than one pair of comparison
                {
                    for (signifIndex in (1:length(ListSignif))) {
                        ##  print(maxTr+(signifIndex-.4)*abs(maxTr-minTr)*.2)
                        segments(y0=maxTr+(signifIndex-.4)*abs(maxTr-minTr)*.2,
                                 x0= ListSignifPosIndex[[signifIndex]][1],x1=ListSignifPosIndex[[signifIndex]][2])
                        text(x=(ListSignifPosIndex[[signifIndex]][1]+ListSignifPosIndex[[signifIndex]][2])/2,y=maxTr+(signifIndex-.1)*abs(maxTr-minTr)*.2,
                             labels=ListSignif[[signifIndex]][1])
                    }
                } else {
                    segments(y0=maxTr+(1-.4)*abs(maxTr-minTr)*.2,
                             x0= 1,x1=2)
                    text(x=1+1/2,y=maxTr+(1-.1)*abs(maxTr-minTr)*.2,
                             labels=ListSignif[1])
                    }
       }
        beeswarm::beeswarm(PP ~ TreatmentFSOM,data=plotDf,add=T,cex=.5,col="red")
       return(pairwisePval)
    })
    ## finish the construction of PvalTable, write csv files
    if(is.matrix(PvalPairwiseTable)) {
        PvalPairwiseTable=as.data.frame(PvalPairwiseTable)
    }
    else if (is.list(PvalPairwiseTable)) {PvalPairwiseTable = as.data.frame(do.call(cbind,PvalPairwiseTable))}
    else {tmpName=names(PvalPairwiseTable)[1]
    PvalPairwiseTable=as.data.frame(t(PvalPairwiseTable))
    row.names(PvalPairwiseTable) = c(tmpName)} ## case of only two treatments
    names(PvalPairwiseTable)=paste("mtcl",colnames(fSOMnbrs)[1:metaclNumber],sep="_")
    par(mfrow=c(1,1),mar=c(3,2,3,1),cex=.5)

    if (ExportData) {
         if(length(MarkerIndex) == 1) {
           write.table(PvalPairwiseTable,NoSpCharForFile(gsub("/","_",paste(Title,"_PairwisePval",Marker,"Metacl.csv",sep=""),fixed=T)),sep=";",col.names = NA)} else {
          if (Norm) {write.table(PvalPairwiseTable,NoSpCharForFile(paste(Title,"_PairwisePvalNormMetacl.csv",sep="")),sep=";",col.names = NA)}
        else {write.table(PvalPairwiseTable,NoSpCharForFile(paste(Title,"_PairwisePvalPercentMetacl.csv",sep="")),sep=";",col.names = NA)}}
    }
    
    DF4lm = data.frame(y=c(fSOMnbrs),metaCl = c(sapply(colnames(fSOMnbrs),function(name){rep(name,length(fSOMnbrs[,1]))})),treat = rep(c(sapply(row.names(fSOMnbrs),function(name){as.character(treatmentTable$Treatment[which(treatmentTable$files == name)])})),length(fSOMnbrs[1,])))
    DF4lm$treat = factor(DF4lm$treat,levels=c(ControlTreatment,setdiff(unique(DF4lm$treat),ControlTreatment)))
    if (Robust) {
        if (length(PvalPairwiseTable[,1]) == 1) {pvalLmMatrix=PvalPairwiseTable}
        else {
          pvalLmMatrix=as.matrix(PvalPairwiseTable)[which(sapply(row.names(PvalPairwiseTable),function(name){
            spName = strsplit(name,split = "-")[[1]];((spName[1] == ControlTreatment) |(spName[2] == ControlTreatment) )})),]
        }
        row.names(pvalLmMatrix)=gsub(paste("^",ControlTreatment,"-",sep=""),"",row.names(pvalLmMatrix))
        row.names(pvalLmMatrix)=gsub(paste("-",ControlTreatment,"$",sep=""),"",row.names(pvalLmMatrix))
    } else {
        pvalLmMatrix = t(do.call(rbind,lapply(colnames(fSOMnbrs),function(metaCl){
          SubDF4Lm = DF4lm[which(DF4lm$metaCl == metaCl),]
          SubDF4Lm$noData = rep(1,length(SubDF4Lm$y))
          retPval = summary(lm(noData ~ treat,data = SubDF4Lm))$coefficient[-1,4] ##create template
          if (length(unique(SubDF4Lm$treat[which(is.finite(SubDF4Lm$y))])) > 2) {
            SignRetPval = summary(lm(y ~ treat,data = DF4lm[which(DF4lm$metaCl == metaCl),]))$coefficient[-1,4]
            for(name in names(SignRetPval)){retPval[[name]] = SignRetPval[[name]]}
          }
          return(retPval)
        })))
        row.names(pvalLmMatrix) = setdiff(unique(DF4lm$treat),ControlTreatment)
    }
    colnames(pvalLmMatrix) = paste("mtcl",colnames(fSOMnbrs),sep= "_")
    pvalLmMatrix = rbind(rep(1,length(fSOMnbrs[1,])),pvalLmMatrix)
    row.names(pvalLmMatrix)[1] = ControlTreatment
    DF4lm$metaCl = factor(DF4lm$metaCl,level=unique(DF4lm$metaCl)) ## to get the right ordering after "by" function
    if (Robust) {
        meanMatrix  = t(by(DF4lm$y,list(DF4lm$metaCl,DF4lm$treat),function(x){median(x,na.rm=T)}))
        pvalLmMatrix=pvalLmMatrix[row.names(meanMatrix),]
    } else {
        meanMatrix  = t(by(DF4lm$y,list(DF4lm$metaCl,DF4lm$treat),function(x){mean(x,na.rm=T)})) }
    attr(meanMatrix,"class") = NULL
    attr(meanMatrix,"call") = NULL
    colnames(meanMatrix) = paste("mtcl",colnames(meanMatrix),sep= "_")
    pvalAnnotationMatrix = apply(pvalLmMatrix,c(1,2),function(x){
    if (!is.finite(x)){return("")}
            else if (x < 0.0001){return("****")}
            else if(x < 0.001){return("***")}
            else if(x < 0.01){return("**")}
            else if(x < 0.05){return("*")} else {return("")}})
    if (length(MarkerIndex) == 1) {
        colMarginSize=20-18*exp(-max(sapply(colnames(meanMatrix),nchar))/10)
        rowMarginSize=20-18*exp(-max(sapply(row.names(meanMatrix),nchar))/10)
        colCex4Plot=exp(-max(sapply(colnames(meanMatrix),nchar))/70)*exp(-length(colnames(meanMatrix))/50)
        rowCex4Plot=exp(-max(sapply(row.names(meanMatrix),nchar))/70)*exp(-length(row.names(meanMatrix))/50)

        if (Robust) {heatTitle = paste("Median MFI of ",PlotLab,sep="")} else {heatTitle = paste("Mean MFI of ",PlotLab,sep="")}
        par(cex.main=exp(-nchar(heatTitle)/70))
        if (ClustHeat) {
            gplots::heatmap.2(meanMatrix,Rowv=F,Colv=T,dendrogram = "column",scale="none",col = heat.colors(100),cellnote = pvalAnnotationMatrix,
                      notecol = "black",trace = "none",cexRow = rowCex4Plot,cexCol=colCex4Plot,density.info="none",main=heatTitle,
                      notecex=.5,margins=c(colMarginSize,rowMarginSize),key.xlab = "",key.title="")
        }
        gplots::heatmap.2(meanMatrix,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = heat.colors(100),cellnote = pvalAnnotationMatrix,
                  notecol = "black",trace = "none",cexRow = rowCex4Plot,cexCol=colCex4Plot,density.info="none",main=heatTitle,
                  notecex=.5,margins=c(colMarginSize,rowMarginSize),key.xlab = "",key.title="")
    } else {
        if (Robust) { heatTitle = paste("Median ",PlotLab,sep="")}
        else {heatTitle = paste("Mean ",PlotLab,sep="")}
        par(cex.main=exp(-max(c(nchar(heatTitle),(nchar(ControlTreatment)+18)))/70))
        heatTitle=paste(heatTitle,"\n(rel. to ",ControlTreatment,", scaled)",sep="")
        meanMatrix=apply(meanMatrix,2,function(x){
          varCol = sd(x,na.rm=T)
          if (varCol == 0) {varCol = 1} ## contant column, no scale
          return((x-x[1])/varCol)
          })
        meanMatrix=meanMatrix[,paste("mtcl_",unique(TreeMetaCl$metaCl),sep="")] ## get the correct ordering
        if (length(meanMatrix[,1]) > 2) {
          meanMatrix = meanMatrix[-1,,drop=F]
          pvalAnnotationMatrix = pvalAnnotationMatrix[-1,,drop=F]}
        colMarginSize=20-18*exp(-max(sapply(colnames(meanMatrix),nchar))/10)
        rowMarginSize=20-18*exp(-max(sapply(row.names(meanMatrix),nchar))/10)
        colCex4Plot=exp(-max(sapply(colnames(meanMatrix),nchar))/70)*exp(-length(colnames(meanMatrix))/50)
        rowCex4Plot=exp(-max(sapply(row.names(meanMatrix),nchar))/70)*exp(-length(row.names(meanMatrix))/50)
        if (ClustHeat) {
          par(cex.main=.5)
            gplots::heatmap.2(meanMatrix,Rowv=F,Colv=T,dendrogram = "column",scale="none",col = gplots::bluered(100),cellnote = pvalAnnotationMatrix,
                      notecol = "black",trace = "none",cexRow = rowCex4Plot,cexCol=colCex4Plot,density.info="none",main=heatTitle,
                      notecex=.5,margins=c(colMarginSize,rowMarginSize),key.xlab = "",key.title="")
        }
            gplots::heatmap.2(meanMatrix,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = gplots::bluered(100),cellnote = pvalAnnotationMatrix,
                      notecol = "black",trace = "none",cexRow = rowCex4Plot,cexCol=colCex4Plot,density.info="none",main=heatTitle,
                      distfun=function(x){dist(t(apply(meanMatrix,2,function(y){scale(y)})))},notecex=.5,margins=c(colMarginSize,rowMarginSize),key.xlab = "",key.title="")
    }
    par(cex.main=1)

    if (length(as.matrix(PvalPairwiseTable)[,1])>1) {
      matrixPval4Heat=apply(as.matrix(PvalPairwiseTable)[,paste("mtcl_",unique(TreeMetaCl$metaCl),sep="")],c(1,2),function(x){
      if (!is.finite(x)) {return(0)}
      else if (x < 0.0001) {return(-log10(0.0001))}
      else {return(-log10(x))}
    })}
    else
    {
      matrixPval4Heat=apply(t(as.matrix(as.matrix(PvalPairwiseTable)[,paste("mtcl_",unique(TreeMetaCl$metaCl),sep="")])),c(1,2),function(x){
        if (!is.finite(x)) {return(0)}
        else if (x < 0.0001) {return(-log10(0.0001))}
        else {return(-log10(x))}
      })
      matrixPval4Heat=rbind(rep(0,length(matrixPval4Heat[1,])),matrixPval4Heat)
      row.names(matrixPval4Heat)=c("",row.names(PvalPairwiseTable)[1])
    }
    if (length(as.matrix(PvalPairwiseTable)[,1])>1) {
      matrixAnnot4Heat=apply(as.matrix(PvalPairwiseTable)[,paste("mtcl_",unique(TreeMetaCl$metaCl),sep="")],c(1,2),function(x){
      if (!is.finite(x)){return("")}
      else if (x < 0.0001){return("****")}
      else if(x < 0.001){return("***")}
      else if(x < 0.01){return("**")}
      else if(x < 0.05){return("*")}
      else if (x < 0.1){return(".")} else {return("")}})
    }
    else {
      matrixAnnot4Heat=apply(t(as.matrix(as.matrix(PvalPairwiseTable)[,paste("mtcl_",unique(TreeMetaCl$metaCl),sep="")])),c(1,2),function(x){
        if (!is.finite(x)){return("")}
        else if (x < 0.0001){return("****")}
        else if(x < 0.001){return("***")}
        else if(x < 0.01){return("**")}
        else if(x < 0.05){return("*")}
        else if (x < 0.1){return(".")} else {return("")}})
      matrixAnnot4Heat=rbind(rep("",length(matrixAnnot4Heat[1,])),matrixAnnot4Heat)
      row.names(matrixAnnot4Heat)=c("",row.names(PvalPairwiseTable)[1])
    }


    maxLogPval = max(unlist(matrixPval4Heat))
    colMarginSize=20-18*exp(-max(sapply(colnames(matrixPval4Heat),nchar))/10)
    rowMarginSize=20-18*exp(-max(sapply(row.names(matrixPval4Heat),nchar))/10)
    colCex4Plot=exp(-max(sapply(colnames(matrixPval4Heat),nchar))/70)*exp(-length(colnames(matrixPval4Heat))/50)
    rowCex4Plot=exp(-max(sapply(row.names(matrixPval4Heat),nchar))/70)*exp(-length(row.names(matrixPval4Heat))/50)
    if (Robust) {
        if (ClustHeat) {
            gplots::heatmap.2(matrixPval4Heat,Rowv=T,Colv=T,dendrogram = "both",scale="none",col = gray(1-((0:100)/100*maxLogPval/(-log10(0.0001)))),
                      trace="none",main="-log10(Dunn p-values)",cexRow = rowCex4Plot,cexCol=colCex4Plot,margins=c(colMarginSize,rowMarginSize),density.info="none",
                      cellnote = matrixAnnot4Heat,notecol = "blue",notecex = .5,key.xlab = "",key.title="")
        }
            gplots::heatmap.2(matrixPval4Heat,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = gray(1-((0:100)/100*maxLogPval/(-log10(0.0001)))),
                      trace="none",cexRow = rowCex4Plot,cexCol=colCex4Plot,main="-log10(Dunn p-values)",margins=c(colMarginSize,rowMarginSize),density.info="none",
                      cellnote = matrixAnnot4Heat,notecol = "blue",notecex = .5,key.xlab = "",key.title="")
    } else {
             if (ClustHeat) {
                 gplots::heatmap.2(matrixPval4Heat,Rowv=T,Colv=T,dendrogram = "both",scale="none",col = gray(1-((0:100)/100*maxLogPval/(-log10(0.0001)))),
                           trace="none",cexRow = rowCex4Plot,cexCol=colCex4Plot,main="-log10(Tukey p-values)",margins=c(colMarginSize,rowMarginSize),density.info="none",
                           cellnote = matrixAnnot4Heat,notecol = "blue",notecex = .5,key.xlab = "",key.title="")
             }
             gplots::heatmap.2(matrixPval4Heat,Rowv=F,Colv=F,dendrogram = "none",scale="none",col = gray(1-((0:100)/100*maxLogPval/(-log10(0.0001)))),
                       trace="none",cexRow = rowCex4Plot,cexCol=colCex4Plot,main="-log10(Tukey p-values)",margins=c(colMarginSize,rowMarginSize),density.info="none",
                       cellnote = matrixAnnot4Heat,notecol = "blue",notecex = .5,key.xlab = "",key.title="") }
    dev.off()

    if (is.table(fSOMnbrs)) {
      DFSizes = as.data.frame(as.matrix.data.frame(fSOMnbrs))
      row.names(DFSizes) = row.names(fSOMnbrs)
      colnames(DFSizes) = colnames(fSOMnbrs)
    } else {DFSizes = as.data.frame(fSOMnbrs)}

    retData=list(DFSizes,as.data.frame(meanMatrix),as.data.frame(PvalPairwiseTable),as.data.frame(pvalLmMatrix))

    if (ExportData) {
      if(length(MarkerIndex) == 1) {
        write.table(DFSizes,NoSpCharForFile(gsub("/","_",paste(Title,"_MFI",Marker,"Metacl.csv",sep=""),fixed=T)),sep=";",col.names = NA)
        write.table(as.data.frame(meanMatrix),NoSpCharForFile(gsub("/","_",paste(Title,"_MedianMFI",Marker,"Metacl.csv",sep=""),fixed=T)),sep=";",col.names = NA)
        write.table(pvalLmMatrix,NoSpCharForFile(gsub("/","_",paste(Title,"_LmPval",Marker,"Metacl.csv",sep=""),fixed=T)),sep=";",col.names = NA)} 
      else {
      if (Norm) {
        write.table(DFSizes,NoSpCharForFile(paste(Title,"_NormSizesMetacl.csv",sep="")),sep=";",col.names = NA)
        write.table(as.data.frame(meanMatrix),NoSpCharForFile(paste(Title,"_MedianNormSizesMetacl.csv",sep="")),sep=";",col.names = NA)
        write.table(pvalLmMatrix,NoSpCharForFile(paste(Title,"_LmPvalNormMetacl.csv",sep="")),sep=";",col.names = NA)}
     else {
       write.table(DFSizes,NoSpCharForFile(paste(Title,"_PercentMetacl.csv",sep="")),sep=";",col.names = NA)
       write.table(as.data.frame(meanMatrix),NoSpCharForFile(paste(Title,"_MedianPercentMetacl.csv",sep="")),sep=";",col.names = NA)
       write.table(pvalLmMatrix,NoSpCharForFile(paste(Title,"_LmPvalPercentMetacl.csv",sep="")),sep=";",col.names = NA)}}
}
    names(retData)=c("Sizes","Median","PvalPairwise","PvalLm")
    return(retData)
}
## Internal tool: replace special characters by "_"
NoSpCharForFile = function(char){return(gsub("[*()#$><%!&|{}\\[\\]?/:@]","_",char,perl=T))}
