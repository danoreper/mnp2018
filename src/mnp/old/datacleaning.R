trim.factor.levels <- function(x)
{
	if (is.factor(x))
	{
		x <- as.factor(string.trim(as.character(x)))
	}
	else if (is.data.frame(x))
	{
		for (i in 1:ncol(x))
		{
			xo <- x[,i]
			xn <- trim.factor.levels(xo) 			
			if (any(as.character(xn)!=as.character(xo), na.rm=TRUE))
			{
				x[,i] <- xn
				cat(colnames(x)[i], " lost '", sep="", paste(
						collapse="'", setdiff(levels(xo), levels(xn))),
						"'\n")
			}
		}
	}
	return (x); 
}

merge.factor.levels <- function(x, lookup)
{
	as.factor(map.eq(as.character(x), lookup))
}

view.factors <- function(d, max.levels=20)
{
	for (i in 1:ncol(d))
	{
		name <- colnames(d)[i]
		x    <- d[,i]
		if (is.factor(x))
		{
			cat(name, ":")
			if (nlevels(x)>max.levels)
			{
				cat("Too many levels to display:", nlevels(x), "\n")
				next
			}
			print(table(x))
		}
	}
}