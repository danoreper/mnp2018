library(data.table)
nu.prob <- function(nu, n, I, J)
{
    aterm = function(nu, n, I, J)
    {
        out = (choose(n, nu))^((J*I)*(I-1))
        return(out)
    }

    if(n %%2 ==1)
    {

        upper = as.integer(round((n-1)/2))
        if(nu>upper)
        {
            return(0)
        }
        return( aterm(nu, n, I, J)/sum(aterm(0:upper, n, I, J)))
        
    } else {
        
        upper = as.integer(round(n/2))
        
        aterm2 = function(nu,n,I,J)
        {
            if(nu>upper)
            {
                return(0)
            }
            if(nu==as.integer(round(n/2)))
            {
                x = (choose(n, n/2)^(2*J))/2
                x = x^(choose(I,2))
                return(x)
            } else {
                return(aterm(nu,n,I,J))
            }
        }
        ##debug(aterm2)
        
        denom = sum(unlist(lapply((0:upper), aterm2, n=n, I=I, J=J)))
        out= (unlist(lapply(nu, aterm2, n=n, I=I, J=J))/denom)
        return(out)
    }
}

get.nu <- function(N,I,J, sz=1)
{
    prob = nu.prob(nu=0:N, N, I, J)
    out = sample(x = 0:N, size=sz, replace = T, prob=prob)
    return(out)
}


generateSynchPerm <- function(df, f1name, f2name, outcomename)
{
    
    I = length(unique(df[[f1name]]))
    J = length(unique(df[[f2name]])) 

    df.perm = data.table(df)
    ##    df.perm[[outcomename]] = sample(df.perm[[outcomename]], size = nrow(df.perm), replace = F)
    df.perm[,srt:=sample(1:.N, size = .N, replace = F),
            by = paste0(eval(parse(text=f1name)), eval(parse(text=f2name)))]
    
    f1.pairs = combn(unique(df.perm[[f1name]]), m=2)
    for(p in 1:ncol(f1.pairs))
    {
        f11 = f1.pairs[1,p]
        f12 = f1.pairs[2,p]
        
        for(f2l in unique(df.perm[[f2name]]))
        {
            df.sub = df.perm[eval(parse(text=f2name))==f2l]
            
            N = min(sum(df.sub[[f1name]]==f11), sum(df.sub[[f1name]]==f12))
            nu = get.nu(N, I, J)
            if(nu>0)
            {
                for(swp in 1:nu)
                {
                    
                    elem1 = df.perm[eval(parse(text=f2name))==f2l & eval(parse(text=f1name)) ==f11 & srt ==swp,][[outcomename]]
                    elem2 = df.perm[eval(parse(text=f2name))==f2l & eval(parse(text=f1name)) ==f12 & srt ==swp,][[outcomename]]
                        
                    df.perm[eval(parse(text=f2name))==f2l & eval(parse(text=f1name)) ==f11 & srt ==swp,][[outcomename]] = elem2
                    df.perm[eval(parse(text=f2name))==f2l & eval(parse(text=f1name)) ==f12 & srt ==swp,][[outcomename]] = elem1
                }
            }
        }
    }
    return(df.perm)
}



n.levels = c(4,2)
n.samples = 4


df = list()
for(i in 1:length(n.levels))
{
    factorname = paste0("f", i)
    df[[i]] = paste0(factorname, 1:n.levels[i])
}
df = do.call(expand.grid, df)
colnames(df) = paste0("f", 1:length(n.levels))

df = df[rep(1:nrow(df), n.samples),]
rownames(df) = NULL

df$ID = 1:nrow(df)
##df$ID = sample(1:nrow(df), replace = F)
df$Y = 1:nrow(df)

numperm = 100
f1name      = "f1"
f2name      = "f2"
outcomename = "Y"


for(i in 1:numperm)
{
    perm1a = generateSynchPerm(df,f1name,f2name,outcomename)
    perm1b = generateSynchPerm(df,f2name,f1name,outcomename)
}



