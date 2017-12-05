


central.portion <- function(x, p=0.5, value=TRUE)
{
    r <- length(x)*0.5*(1 + p*c(-1,1))
    i <- ceiling(r[1]):floor(r[2])
    if (value) return ( x[ ceiling(r[1]):floor(r[2])] )
    else return (i)
}


first <- function(x)
{
    return(x[1])
}

indep.factor.cov <- function(fact)
{
    X <- incidence.matrix(fact)
    X %*% t(X)
}

last <- function(x)
{
    return(x[length(x)])
}

map <- function(x, lookup=NULL, from=NULL, to=NULL)
{
    stop("use map.eq()\n")
}

normalize <- function(x, ...)
{
    qnorm(rank(x)/(length(x)+1), ...)
}


one <- function()
{
	par(mfrow=c(1,1))
	invisible()
}

pop.back <- function(x)
{
    r <- x[length(x)]
    eval.parent(substitute(x<-x[-length(x)]))
    return (r)
}

pop.front <- function(x)
{
    r <- x[1]
    eval.parent(substitute(x<-x[-1]))
    return (r)
}

rev.rank <- function(x, ...)
{
    length(x) - rank(x, ...) + 1
}

squarelist2matrix <- function(x)
{
	warning("Deprecated, use sapply(, function(x){x})\n")
	ncol <- length(x)
	nrow <- unique(lapply(x, length))
	if (1!=length(nrow))
	{
		stop("Cannot turn ragged list into a square matrix. ",
				"Got column lengths: {", paste(collapse=",", nrow),
				"}\n")
	}
	mat <- matrix(ncol=ncol, nrow=nrow)
	for (i in 1:ncol)
	{
		mat[,i] <- x[[i]]
	}
	mat
}


SS <- function(x, mean.correct=T)
{
    if (mean.correct)
    {
        x <- x - mean(x)
    }
    x %*% x
}

strcat <- function(vec, sep="")
{
	vec <- as.character(vec)
	paste(vec, collapse=sep)
}

substr.pattern <- function(prefix, wanted, suffix, text)
{
	out <- rep(NA, length(text))
	is.match <- igrep(paste(prefix,wanted,suffix,sep=""), text, perl=TRUE)
	for (i in which(is.match))
	{
		s <- sub(prefix, "", text[i])
		r <- regexpr(wanted, s, perl=TRUE)
		out[i] <- substr(s, r, attr(r, "match.length")+r-1)
	}
	out
}



unique.factor.pairs <- function(data)
{
	facts <- colnames(data)
	n     <- length(facts) 
	cmat <- matrix(NA, nrow=n, ncol=n)
		
	for (i in 1:n)
	{
		cmat[i,i] <- length(unique(data[,facts[i]]))
	}
	for (i in 1:(n-1))
	{
		for (j in (i+1):n)
		{
			nu <- nrow(unique(data[,c(facts[i], facts[j])]))
			cmat[i,j] <- nu
			cmat[j,i] <- nu - max(cmat[i,i], cmat[j,j])
		}
	}
	
	colnames(cmat) <- facts
	rownames(cmat) <- facts
	cmat
}


# random stuff

shannon <- function(x)
{
    p <- x/sum(x)
    p <- p[0!=p]
    - sum( p*log(p,2))
}

sic <- function(x1, x0=NULL, rescale=NULL)
# selective information content
# I = Shannon(before) - Shannon(after)
{
    s0 <- NULL
    if (is.null(x0))
    {
        s0 <- log(length(x1),2)
    }
    else
    {
        s0 <- shannon(x0)
    }
    ic <-  s0 - shannon(x1)

    if (!is.null(rescale) & 2==length(rescale))
    {
        ic <- diff(rescale)*ic/log(length(x1),2) + rescale[1]
    }
    ic
}


#--------------------------------------------------------------------------
# Longer functions (to be moved at some point)

calc.high.region.single <- function(x, y, center.index, threshold)
{
    if (y[center.index] < threshold)
    {
        return (data.frame(
                lower.index=NA,
                lower.bound=NA,
                upper.index=NA,
                upper.bound=NA,
                length=0,
                center.index=center.index,
                center=x[center.index],
                threshold=threshold))
    }

    # find lower bound
    for (i in center.index:1)
    {
        if (is.na(x[i]) | is.na(y[i])) next
        if (y[i] >= threshold) lower.index <- i
        else break
    }

    # find upper bound
    for (i in center.index:length(x))
    {
        if (is.na(x[i]) | is.na(y[i])) next
        if (y[i] >= threshold) upper.index <- i
        else break
    }

    return (data.frame(
            lower.index=lower.index,
            lower.bound=x[lower.index],
            upper.index=upper.index,
            upper.bound=x[upper.index],
            length= x[upper.index] - x[lower.index],
            center.index=center.index,
            center=x[center.index],
            threshold=threshold))
}

calc.high.region <- function(x, y, center.index, thresholds,
        row.names=1:length(thresholds))
{
    regions.df <- NULL
    for (i in 1:length(thresholds))
    {
        regions.df <- rbind(regions.df,
                calc.high.region.single(x, y, center.index, thresholds[i]))
    }
    rownames(regions.df) <- row.names
    return (regions.df)
}
