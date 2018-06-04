library(data.table)

##Add num shuffles with seed
generateSynchPerm <- function(df, f1name, f2name)
{
##    set.seed(seed)
    I = length(unique(df[[f1name]]))
    J = length(unique(df[[f2name]])) 

    df.perm = data.table(df)
    df.perm$neworder = 1:nrow(df.perm)

    df.perm[,srt:=sample(1:.N, size = .N, replace = F),
            by = paste0(eval(parse(text=f1name)), eval(parse(text=f2name)))]
    
    f1.pairs = combn(unique(df.perm[[f1name]]), m=2)


    N = min(table(df[[f1name]], df[[f2name]]))
    ##nu = get.nu(N, I, J)
    nu = floor(N/2)

    for(p in 1:ncol(f1.pairs))
    {
        f11 = f1.pairs[1,p]
        f12 = f1.pairs[2,p]

        ## N = min(sum(df.sub[[f1name]]==f11), sum(df.sub[[f1name]]==f12))
        for(f2l in unique(df.perm[[f2name]]))
        {
            df.sub = df.perm[eval(parse(text=f2name))==f2l]
            
        
            ## nu = get.nu(N, I, J)
            if(nu>0)
            {
                for(swp in 1:nu)
                {
                    
                    elem1 = df.perm[eval(parse(text=f2name))==f2l & eval(parse(text=f1name)) ==f11 & srt ==swp,][["neworder"]]
                    elem2 = df.perm[eval(parse(text=f2name))==f2l & eval(parse(text=f1name)) ==f12 & srt ==swp,][["neworder"]]
                        
                    df.perm[eval(parse(text=f2name))==f2l & eval(parse(text=f1name)) ==f11 & srt ==swp,][["neworder"]] = elem2
                    df.perm[eval(parse(text=f2name))==f2l & eval(parse(text=f1name)) ==f12 & srt ==swp,][["neworder"]] = elem1
                }
            }
        }
    }
    return(df.perm$neworder)
}

get.nu <- function(N,I,J, sz=1)
{
    prob = nu.prob(nu=0:N, N, I, J)
##    browser()
    out = sample(x = 0:N, size=sz, replace = T, prob=prob)
    return(out)
}


nu.prob <- function(nu, n, I, J)
{
    if(length(n)>1)
    {
        stop()
    }
    
    aterm = function(nu, n, I, J)
    {
        out = (choose(n, nu))^((J*I)*(I-1))
        return(out)
    }

    if(n %%2 ==1)
    {
        upper = as.integer(round((n-1)/2))
        denom = sum(unlist(lapply((0:upper), aterm, n=n, I=I, J=J)))
        return( (nu<=upper) * aterm(nu, n, I, J)/denom)
        
    } else {

        upper = as.integer(round(n/2))
        
        aterm2 = function(nu,n,I,J)
        {
            if(nu==as.integer(round(n/2)))
            {
                x = (choose(n, n/2)^(2*J))/2
                x = x^(choose(I,2))
                return(x)
            } else {
                x = (((nu<=upper)*aterm(nu,n,I,J)))
                return(x)
            }
        }
        
        denom = sum(unlist(lapply((0:upper), aterm2, n=n, I=I, J=J)))
        out= (unlist(lapply(nu, aterm2, n=n, I=I, J=J))/denom)
        return(out)
    }
}





## n.levels = c(4,2)
## n.samples = 3


## df = list()
## for(i in 1:length(n.levels))
## {
##     factorname = paste0("f", i)
##     df[[i]] = paste0(factorname, 1:n.levels[i])
## }
## df = do.call(expand.grid, df)
## colnames(df) = paste0("f", 1:length(n.levels))

## df = df[rep(1:nrow(df), n.samples),]
## rownames(df) = NULL

## df$ID = 1:nrow(df)
## ##df$ID = sample(1:nrow(df), replace = F)
## df$Y = 1:nrow(df)

## numperm = 1
## f1name      = "f1"
## f2name      = "f2"
## outcomename = "Y"


## for(i in 1:numperm)
## {
##     perm1a = generateSynchPerm(df,f1name,f2name)
##     perm1b = generateSynchPerm(df,f2name,f1name)
## }



