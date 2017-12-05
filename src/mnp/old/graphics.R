
A4PortraitWidth <- function(height)
{
	height/sqrt(2)
}

A4PortraitHeight <- function(width)
{
	sqrt(2)*width
}


build.mat <- function(x, y, g)
{
	ux <- sort(unique(x))
	ug <- sort(unique(g))
	mat <- matrix(nrow=length(ux), ncol=length(ug), )
	
	for (ix in 1:length(ux))
	{
		for (ig in 1:length(ug))
		{
			i <- which(x==ux[ix] & g==ug[ig])
			if (0==length(i)) next
			if (1<length(i))
			{
				stop("Cannot process replicate combinations: ",ux[ix],",",ug[ig],"\n")
			}
			mat[ix,ig] <- y[i]
		}
	}
	rownames(mat) <- ux
	colnames(mat) <- ug
	mat
}

plotmany <- function(formula, data=NULL, ...)
{
	form.string <- deparse(formula)
	y.form <- sub("~.*","", form.string)
	x.form <- sub("\\|.*", "", sub(".*~","", form.string))
	g.form <- sub(".*~.*\\|", "", form.string)
	
	y <- apply.transform(y.form, data)
	x <- apply.transform(x.form, data)
	g <- apply.transform(g.form, data)
	
	mat <- build.mat(x,y,g)
	thing <- matplot(as.numeric(rownames(mat)), mat, ...)
	invisible(mat)
}

my.truehist <- function (data, nbins = "Scott", h, x0 = -h/1000, breaks, prob = TRUE,
		xlim = range(breaks), ymax = max(est), col = 5, xlab = deparse(substitute(data)),
		bty = "n", add=FALSE, ...)
# truehist() from library MASS but with an "add" argument
{
	plot.truehist <- function(breaks, est, xlim, ymax, bty, xlab,
			ylab = "", density = NULL, angle = 45, col = NULL, border = NULL,
			lty = NULL, lwd = par("lwd"), ...) {
		n <- length(breaks)
		
		if (!add)
		{
			plot(xlim, c(0, ymax), type = "n", xlab = xlab, ylab = ylab,
					bty = bty, ...)
		}
		rect(breaks[-n], 0, breaks[-1], est, density = density,
				angle = angle, col = col, border = border, lty = lty,
				lwd = lwd)
	}
	xlab
	data <- data[is.finite(data)]
	if (missing(breaks)) {
		if (missing(h)) {
			if (is.character(nbins))
				nbins <- switch(casefold(nbins), scott = nclass.scott(data),
						"freedman-diaconis" = , fd = nclass.FD(data))
			h <- diff(pretty(data, nbins))[1]
		}
		first <- floor((min(data) - x0)/h)
		last <- ceiling((max(data) - x0)/h)
		breaks <- x0 + h * c(first:last)
	}
	if (any(diff(breaks) <= 0))
		stop("breaks must be strictly increasing")
	if (min(data) < min(breaks) || max(data) > max(breaks))
		stop("breaks do not cover the data")
	db <- diff(breaks)
	if (!prob && sqrt(var(db)) > mean(db)/1000)
		warning("Uneven breaks with prob = FALSE will give a misleading plot")
	bin <- cut(data, breaks, include.lowest = TRUE)
	est <- tabulate(bin, length(levels(bin)))
	if (prob)
		est <- est/(diff(breaks) * length(data))
	plot.truehist(breaks, est, xlim, ymax, bty = bty, xlab = xlab,
			col = col, ...)
	invisible()
}

plot.binary.interact <- function(form, data, lty=1:10, col="black",
		ylim=c(0,1),
		pch=20, main=NULL, ylab=NULL, xlab=NULL,
		legend.xy  = NULL,
		legend.txt = NULL,
		uex=0.1,
		ci=NULL, ...)
# f( case ~ fact1 * fact2 , data)
{
	spl        <- split.formula(form)
	response   <- spl$response
	y          <- data[,response]
	predictors <- string.trim(unlist(strsplit(spl$predictors, "\\*")))
	x1    <- as.factor(data[,predictors[1]])
	x2    <- as.factor(data[,predictors[2]])
	
	n1 <- length(levels(x1))
	n2 <- length(levels(x2))
	if (is.null(main)) main <- paste(response, "~", predictors[1], "x", predictors[2])
	if (is.null(xlab)) xlab <- predictors[2]
	if (is.null(ylab)) ylab <- response
	
	pad <- 0.25
	plot(c(1-pad,n2+pad), c(0,1), type="n", axes=FALSE,
			main=main, ylab=ylab, xlab=xlab, ylim=ylim, ...)
	axis(1, labels=levels(x2), at=1:n2)
	axis(2, las=1)
	
	if (length(lty)<n1) lty <- rep(lty, length.out=n1)
	if (length(col)<n1) col <- rep(col, length.out=n1)
	if (length(pch)<n1) pch <- rep(pch, length.out=n1)
	
	se.factor <- ifelse(is.null(ci), 1, qnorm(1-(1-ci)/2))
	cat(se.factor,"\n")
	for (i1 in 1:length(levels(x1)))
	{
		a1 <- levels(x1)[i1]
		y.p  <- NULL
		y.se <- NULL
		for (a2 in levels(x2))
		{
			i <- x1==a1 & x2==a2
			i <- i & !is.na(x1) & !is.na(x2)
			p  <- mean(y[i], na.rm=TRUE)
			se <- p*(1-p)/sqrt(sum(i))
			push.back(y.p, p)
			push.back(y.se, se*se.factor)
		}
		points.errbar(1:n2, y.p, y.plus=y.se, y.minus=y.se,
				pch=pch[i1], col=col[i1], uex=uex)
		lines(1:n2, y.p, lty=lty[i1], col=col[i1])
	}
	if (is.null(legend.xy)) legend.xy <- c(1, ylim[2])
	if (is.null(legend.txt)) legend.txt <- levels(x1)
	
	legend(legend.xy[1], legend.xy[2],
			legend=legend.txt,
			pch=pch, lty=lty, col=col)
}


plot.table <- function(data, main="", indent=0.1, box=TRUE, line=TRUE, header=colnames(data), ...)
{
	plot(c(0,ncol(data)), c(0,nrow(data)), type="n", axes=FALSE, ylab="", xlab="", main=main)
	if (line) abline(h=nrow(data)-0.5)
	if (box) box()
	for (ic in 1:ncol(data))
	{
		text(ic-1+ncol(data)*indent, nrow(data):0, c(header[ic], data[,ic]),...)
	}
}

plot.text <- function(text, add=FALSE, mar=c(2,2,2,2), ...)
{
	lines <- unlist(strsplit(text, "\n"))
	n <- length(lines)
	oldmar <- par(mar=mar); on.exit(par(mar=oldmar))
	if (!add) plot(0:1,c(0,n), type="n", axes=FALSE, xlab="", ylab="")
	text(x=rep(0,n), y=(n:1)-0.5, adj=c(0,1), labels=lines, ...)
} 

plot.ci <- function(midvals, narrow.intervals, wide.intervals, 
		names=1:length(midvals),
		add=FALSE,
		axes=sides(top=TRUE, bottom=TRUE, left=yaxis, right=FALSE),
		xlab="Estimate",
		xlab.line=2.5,
		xlim=NULL,
		ylab="",
		yaxis=TRUE,
		ylim=c(0, length(midvals)),
		name.line=4,
		pch.midvals=19,
		col="black",
		col.midvals=col,
		cex.labels=1,
		type="p",
		name.margin=6.1,
		title.margin=4.1,
		bottom.margin=5.1,
		right.margin=2.1,
		mar=sides(left=name.margin, bottom=bottom.margin, top=title.margin, right=right.margin),
		mar.update=sides(),
		before.data=function(){},
		...)
# Example: plot.ci( c(0,10), narrow.intervals=rbind(c(-1,1), c(8,12)), wide.intervals=rbind(c(-3,4), c(5,15)), names=c("Fred", "Barney"))
{
	nvals <- length(midvals)
	col.midvals <- rep(col.midvals, length.out=nvals)
	y.pos <- (1:nvals)-0.5
	if (!add)
	{
		if (is.null(xlim))
		{
			xlim <- range(c(wide.intervals,narrow.intervals,midvals), na.rm=TRUE)
			xlim <- c(-1,1) * diff(xlim)*0.1 + xlim
		}
		
		mar <- c(bottom.margin, name.margin, title.margin, right.margin)+0.1
		mar=update.sides(mar, mar.update)
		oldmar <- par(mar=mar); on.exit(par(mar=oldmar))
		
		plot(xlim, ylim, type="n", axes=FALSE, ylab=ylab, ylim=ylim, xlab="", ...)
		title(xlab=xlab, line=xlab.line)
		if (axes["top"])
		{
			axis(3, line=-1)
		}
		if (axes["bottom"])
		{
			axis(1)
		}
		if (axes["left"])
		{ 
			axis(2, at=y.pos, labels=rev(names), las=1, lty=0, hadj=0, line=name.line, cex.axis=cex.labels)
		}
	}
	before.data()
	if ("p"==type)
	{
		for (i in 1:nvals)
		{
			pos <- nvals-i + 0.5
			lines(wide.intervals[i,], rep(pos,2))
			lines(narrow.intervals[i,], rep(pos,2), lwd=3)
			points(midvals[i], pos, pch=pch.midvals, col=col.midvals[i])
		}
	}
	invisible(rev(y.pos))
}


plot.corr <- function(mat,
		range=c(-1,1),
		border.lwd=0,
		num.colors=200,
		main="correlations (/=+ve, \\=-ve)",
		line.col="white",
		name.margin=6.1,
		name.line=4,
		mar=sides(bottom=5.1, left=name.margin, top=4.1, right=2.1),
		mar.update=sides(),
		names=colnames(mat),		
		darkness.extent=1,
		positive.increasing=TRUE,
		show.range=TRUE,
		show.sign=TRUE,
		xlab="correlation",
		xlab.line=2)
{
	num.colors <- 2*ceiling(num.colors/2)
	if (TRUE)
	{
		darkness <- c( ((num.colors/2):1)/(num.colors/2), 0, (1:(num.colors/2))/(num.colors/2))
		darkness <- darkness*darkness.extent
		colors <- gray(1-darkness)
	}
	else
	{
		colors <- rainbow(num.colors+1, start=4/6, end=0)
	}
	limits <- 0:num.colors * (diff(range)/num.colors) + min(range)
	n <- nrow(mat)
	
	# margins
	mar=update.sides(mar, mar.update)
	oldmar <- par(mar=mar); on.exit(par(mar=oldmar))
	
	y.offset <- 0.5
	
	plot(c(0,n), c(0,n+y.offset), axes=FALSE, type="n", ylab="", xlab="")
	title(xlab=xlab, line=xlab.line)
	title(main=main, line=xlab.line)
	image(0:n, c(0:n)+y.offset, mat[,n:1], zlim=range, col=colors, add=TRUE)
	
	if (show.sign)
	{
		x0=matrix(nrow=n, ncol=n)
		y0=matrix(nrow=n, ncol=n)
		x1=matrix(nrow=n, ncol=n)
		y1=matrix(nrow=n, ncol=n)
		for (j in 1:n)
		{
			for (k in 1:n)
			{
				x <- mat[j,k]
				x.start <- j-1
				x.end   <- j
				y.start <- y.offset + n-k
				y.end   <- y.offset + n-k+1
				
				x0[j,k]=x.start
				x1[j,k]=x.end
				if ((x>0 & positive.increasing) | (x<0 & !positive.increasing))
				{
					y0[j,k]=y.start
					y1[j,k]=y.end
				}
				else if ((x<0 & positive.increasing)|(x>0 & !positive.increasing))
				{
					y0[j,k]=y.end
					y1[j,k]=y.start
				}
			}
		}
		segments(c(x0), c(y0), c(x1), c(y1), col=line.col)
	}
	if (0<border.lwd)
	{
		abline(h=(0:n)+y.offset, lwd=border.lwd, col="white")
		abline(v=0:n, lwd=border.lwd, col="white")
	}
	axis(3, at=(1:(n))-.5, labels=names, las=3, lty=0, line=-1)
	axis(2, at=(1:(n))-0.5+y.offset, labels=rev(names), las=1, lty=0, hadj=0, line=name.line)
	if (show.range)
	{
		points((0:num.colors)/num.colors*n, rep(0, num.colors+1), col=colors, pch=19)
		axis(1, at=c(0,n/2, n), labels=as.character(signif(c(min(range), mean(range), max(range))),digits=3))
	}
}


# plot.ci <- function(midvals, narrow.intervals, wide.intervals, names=1:length(midvals),
#     add=FALSE,
#     xlab="Estimate",
#     ylab="",
#     yaxis=TRUE,
#     name.margin=6,
#     name.line=4,
#     pch.midvals=19,
#     col="black",
#     col.midvals=rep(col,length(midvals)),
#     type="p",
#     ...)
# {
#   nvals <- length(midvals)
#   y.pos <- (1:nvals)-0.5
#   if (!add)
#   {
#     lim <- range(c(wide.intervals,narrow.intervals,midvals), na.rm=TRUE)
#     lim <- c(-1,1) * diff(lim)*0.1 + lim
#   
#     mar <- c(5, name.margin, 4, 2)+0.1
#     oldmar <- par(mar=mar); on.exit(par(mar=oldmar))
#     plot(lim, c(0,nvals+0.5), type="n", axes=FALSE, xlab=xlab, ylab=ylab, ...)
#     axis(1)
#     axis(3, line=-1)
#     if (yaxis)
#     {
#       y <- axis(2, at=y.pos, labels=rev(names), las=1, lty=0, hadj=0, line=name.line)
#     }
#   }
#   if ("p"==type)
#   {
#     for (i in 1:nvals)
#     {
#       pos <- nvals-i + 0.5
#       lines(wide.intervals[i,], rep(pos,2))
#       lines(narrow.intervals[i,], rep(pos,2), lwd=3)
#       points(midvals[i], pos, pch=pch.midvals, col=col.midvals[i])
#     }
#   }
#   invisible(rev(y))
# }


panel.errbar <- function(x,
		y=NULL,
		y.plus=0,
		y.minus=0,
		y.lower=y-y.plus,
		y.upper=y+y.minus,
		uex = 0.2,
		umbrella = uex*(range(x)%*%c(1,-1))/length(x),
		col="black",
		lwd=1,
		...)
{
	if (!is.null(y)) panel.points(x,y,col=col,...)
	for (i in 1:length(x))
	{
		panel.lines(rep(x[i],2), c(y.lower[i], y.upper[i]), col=col, lwd=lwd)
		panel.lines(x[i] + umbrella*c(-1,1), rep(y.lower[i],2), col=col, lwd=lwd)
		panel.lines(x[i] + umbrella*c(-1,1), rep(y.upper[i],2), col=col, lwd=lwd)
	}
}


points.errbar <- function(x,
		y=NULL,
		y.plus=0,
		y.minus=0,
		y.lower=y-y.minus,
		y.upper=y+y.plus,
		uex = 0.2,
		umbrella = c(uex*(range(x)%*%c(1,-1))/length(x)),
		col="black",
		lwd=1,
		...)
{
	col <- rep(col, length.out=length(y))
	segments(x, y.lower, x, y.upper, col=col, lwd=lwd)
	
	if (!is.na(umbrella))
	{
		segments(x - umbrella, y.lower, x + umbrella, y.lower, col=col, lwd=lwd)
		segments(x - umbrella, y.upper, x + umbrella, y.upper, col=col, lwd=lwd)
	}
	points(x,y,col=col,...)
}

rescaled.axis <- function(side,
		m=1,
		c=0,
		fun=function(x){c+m*x},
		inv=function(y){y/m-c},
		lim=NULL,
		unit="",
		at = NULL,
		labels=NULL,
		...)
{
	if (is.null(at) | is.null(labels))
	{
		if (is.null(lim))
		{
			if (1==side | 3==side)
			{
				lim <- par("usr")[c(1,2)]
			}
			else
			{
				lim <- par("usr")[c(3,4)]
			}
		}
		labels <- pretty(sort(fun(lim)))
		at <- inv(labels)
	}
	axis(side, at=at, labels=paste(labels,unit,sep=""), ...)
}

sides=function(default=NA, bottom=default, left=default, top=default, right=default)
{
	x=c(bottom, left, top, right)
	names(x)=c("bottom", "left", "top", "right")
	x
}

update.sides=function(old=par("mar"), new=rep(NA, 4))
{
	old[!is.na(new)]=new[!is.na(new)]
	old
}

