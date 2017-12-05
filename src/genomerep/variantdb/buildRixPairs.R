# TODO: Add comment
# 
# Author: doreper
###############################################################################
library("data.table")
library(reshape2)

get.matnut.RIX.strains <- function()
{
    x = read.table("../data/offspring_behavior/AllMice_Behavior_GeneExpression5_11_15.csv", sep=",", header=T)
    crosses = unique(paste0(pmin(x$Sire.Line, x$Dam.Line), ".", pmax(x$Sire.Line, x$Dam.Line)))
    source("./idConversion2.R")
    idConversion$init(probsDir)
    parents = idConversion$toParents(crosses)
    strain1 = unname(idConversion$strainNames[parents$p1])
    strain1 = c(strain1, "NOD_ShiLtJ")
    strain2 = unname(idConversion$strainNames[parents$p2])
    strain2 = c(strain2, "C57BL6J")
    return(list(strain1=strain1, strain2=strain2))
}



create.RIX.table <- function(strain1, strain2)
{
    dbSendQuery(con, "drop table if exists temp_strain;")
    q = paste0("create table temp_strain(strain_name_1 varchar(30), strain_name_2 varchar(30)) engine = ", engine, ";");
    dbSendQuery(con, q)
    values  = paste(paste0("(","'",strain1,"'", ", '", strain2,"')"),collapse=", ") 
    dbSendQuery(con, paste0("INSERT IGNORE INTO temp_strain (strain_name_1, strain_name_2) VALUES ", values,";"))
    dbCommit(con)
    
    dbSendQuery(con, "drop table if exists RIX;")
    query = "create table RIX(prob float(3) unsigned not null) select 
		s1.strain_id as strain_id_1, 
		s2.strain_id as strain_id_2, 
		ilh1.variant_id as variant_id, 
		ilh1.allele_index as allele_index_1,  
		ilh2.allele_index as allele_index_2, 
		ilh1.prob*ilh2.prob as prob 

		from inbred_line_allele_hap ilh1 
		inner join strain s1 on s1.strain_id = ilh1.strain_id 
		
		
		inner join inbred_line_allele_hap ilh2 on ilh2.variant_id=ilh1.variant_id 
		inner join strain s2 on s2.strain_id=ilh2.strain_id 
		
		inner join temp_strain on s1.strain_name=temp_strain.strain_name_1 and s2.strain_name=temp_strain.strain_name_2;"


    pracma::tic()
    dbSendQuery(con, query)
    
    ##filters out all variants that are not segregating between any of the parents of any RIX
    createIndex(con, "RIX", c("allele_index_1", "allele_index_2"))
    dbSendQuery(con, "create table vs select distinct variant_id from RIX where allele_index_1!=allele_index_2;")
    createIndex(con, "vs", c("variant_id"))
    
    dbSendQuery(con, "create table RIX2 select RIX.* from RIX, vs where RIX.variant_id=vs.variant_id;")
    dbSendQuery(con, "drop table vs;")
    dbSendQuery(con, "drop table RIX;")
    dbSendQuery(con, "rename table RIX2 to RIX;")
    
    ##this doesnt actually seem to happen, so not going to worry about it
    ##dbSendQuery(con, "create table vs select variant_id, count(distinct concat(ifnull(allele_index_1,'n'), ifnull(allele_index_2,'n'))) as counts from RIX2  having counts=1;")

    createIndexShort(con, "RIX", c("strain_id_1","strain_id_2","variant_id","allele_index_1", "allele_index_2","prob"), "ix1")

    print("built RIXes")
    pracma::toc()
}


writeRIX_info_tofile <- function(whereclause, con, fname=NULL) 
{
    query = paste0(
        "select s1.strain_name as strain_name_1, 
					s2.strain_name as strain_name_2, 
					variant.variant_id, 
					chrom, 
					pos, 
					gene.gene_name, 
					transcript.transcript_name, 
					c1.consequence_name as consequence_1,
					c2.consequence_name as consequence_2, 
					a1.allele           as allele_1,
					a2.allele           as allele_2,
					f1.fixed            as fixed_1,
					f2.fixed            as fixed_2,
					prob		    as prob
					from RIX 
					
					inner join strain s1 on s1.strain_id = RIX.strain_id_1 
					inner join strain s2 on s2.strain_id = RIX.strain_id_2 
					inner join variant on RIX.variant_id = variant.variant_id
					
					inner join inbred_line_fixed as f1 on RIX.strain_id_1 = f1.strain_id and RIX.variant_id = f1.variant_id
					inner join inbred_line_fixed as f2 on RIX.strain_id_2 = f2.strain_id and RIX.variant_id = f2.variant_id
					
					inner join consequence_info as varinfo_1 on RIX.variant_id=varinfo_1.variant_id and RIX.allele_index_1 = varinfo_1.allele_index 
					inner join consequence c1 on c1.consequence_id = consequence_info_1.consequence_id
					inner join allele a1 on a1.variant_id = RIX.variant_id and a1.allele_index = RIX.allele_index_1 
					
					inner join consequence_info as consequence_info_2 on RIX.variant_id=consequence_info_2.variant_id and consequence_info_2.transcript_id = consequence_info_1.transcript_id and RIX.allele_index_2 = consequence_info_2.allele_index 
					inner join consequence c2 on c2.consequence_id = consequence_info_2.consequence_id
					inner join allele a2 on a2.variant_id = RIX.variant_id and a2.allele_index = RIX.allele_index_2 
					
					inner join gene on gene.gene_id = consequence_info_1.gene_id 
					inner join transcript on transcript.transcript_id =consequence_info_1.transcript_id ",
			whereclause, ";")

    ##TODO should we be able to get the transcript and gene wihout having an allele (eg NA case)? need to refactor consequence_info building if that is so, make a separate table for variant->(gene,transcript)
    pracma::tic()
    x = dbGetQuery(con, query)
    print(length(unique(x$variant_id)))
    pracma::toc()
             
    pracma::tic()
    if(!is.null(fname))
    {
        write.table(x, fname, row.names = F, sep="\t")
    }
    pracma::toc()
    return(x)
}




strains = get.matnut.RIX.strains()
create.RIX.table(strains$strain1, strain$strain2)

##write out various subset of the RIX data in long table form, depending on the where clause
sharefolder = "../doreper_matnut_tempshare/"
dir.create(sharefolder) 
hetfolder = file.path(sharefolder,"het")
dir.create(hetfolder)
allfolder = file.path(sharefolder,"all")
dir.create(allfolder)

whereclause = ""
fname = file.path(allfolder, "alldb.csv")	
fname =NULL
z = writeRIX_info_tofile(whereclause = whereclause, con = con, fname = fname)

whereclause = "where (RIX.allele_index_1!=RIX.allele_index_2 or (RIX.allele_index_1 is null or RIX.allele_index_2 is null))"
fname = file.path(hetfolder, "hetdb.csv")	
z = writeRIX_info_tofile(whereclause = whereclause, con = con, fname = fname)



whereclause = "where ((RIX.allele_index_1!=RIX.allele_index_2 or 
		(RIX.allele_index_1 is null or RIX.allele_index_2 is null)) and 
		(c1.consequence_name like '%stop%' or c2.consequence_name like '%stop%')) "
fname = file.path(hetfolder, "hetstop_db.csv")	
z = writeRIX_info_tofile(whereclause = whereclause, con = con, fname = fname)


chroms = dbGetQuery(con, "Select distinct chrom from variant;")$chrom
for (chrom in chroms)
{
    print(chrom)
    whereclause = paste0("where chrom = '",chrom,"'")
    fname = file.path(allfolder, paste0("chrom_",chrom,".csv"))
    writeRIX_info_tofile(whereclause, con, fname)
}

for (chrom in chroms)
{
    print(chrom)
    whereclause = "where (RIX.allele_index_1!=RIX.allele_index_2 or (RIX.allele_index_1 is null or RIX.allele_index_2 is null))"
    whereclause = paste0(whereclause, " and chrom = '",chrom,"'")
    fname = file.path(hetfolder, paste0("chrom_",chrom,".csv"))
    writeRIX_info_tofile(whereclause, con, fname)
}








df = fread(file.path(allfolder, "alldb.csv"), header=T, sep="\t")
df$gene_name = as.factor(df$gene_name)
df$transcript_name = as.factor(df$transcript_name)
df$strain_name_1 = as.factor(df$strain_name_1)
df$strain_name_2 = as.factor(df$strain_name_2)
df$consequence_1 = as.factor(df$consequence_1)
df$consequence_2 = as.factor(df$consequence_2)
df$allele_1      = as.factor(df$allele_1)
df$allele_2      = as.factor(df$allele_2)
setkey(df, "transcript_name")

#TODO this block wont be necessary once we preprocess the founders table properly to 
#only inlcude the right subset of transcripts that are on the imprinted strand
imprintedExonFeatures = rtracklayer::import("../output/imprintedExonFeatures.gtf")
trans =unique(imprintedExonFeatures$transcript_id)
trans = intersect(df$transcript_name, trans)
df = df[J(trans)]

#setkey(df, "strain_name_1", "strain_name_2", "variant_id", "transcript_name", "allele_1", "allele_2")
#df[,prob := sum(prob), by=c("strain_name_1", "strain_name_2", "variant_id", "transcript_name","allele_1","allele_2")]



setkey(df, "strain_name_1", "strain_name_2", "variant_id", "transcript_name")
df[,ord := order(prob,decreasing=T), by=c("strain_name_1", "strain_name_2", "variant_id", "transcript_name")]



topranked = df[df$ord==1]
topranked$prob=round(100*topranked$prob)
topranked$diplotype       = as.factor(ifelse((topranked$allele_1==topranked$allele_2),           ".", paste0(topranked$allele_1,"/", topranked$allele_2))) 
topranked$consequences    = as.factor(ifelse((topranked$consequence_1==topranked$consequence_2), ".", paste0(topranked$consequence_1,"/", topranked$consequence_2)))
topranked$certainty       = as.factor(paste0(topranked$prob, ",", paste0(topranked$fixed_1,"/", topranked$fixed_2)))

dips = data.table(dcast(topranked, "variant_id + chrom + pos +gene_name + transcript_name  ~strain_name_1+strain_name_2", value.var=c("diplotype")))
colnamesbool = !(colnames(dips) %in% c("variant_id" ,"chrom", "pos", "gene_name", "transcript_name"))
newnames = colnames(dips)
newnames[colnamesbool] = paste0(newnames[colnamesbool], "_diplo")
setnames(dips,colnames(dips), newnames)
setkey(dips, "variant_id","transcript_name")

csqs = data.table(dcast(topranked, "variant_id + transcript_name  ~strain_name_1+strain_name_2", value.var=c("consequences")))
colnamesbool = !(colnames(csqs) %in% c("variant_id" ,"transcript_name"))
newnames = colnames(csqs)
newnames[colnamesbool] = paste0(newnames[colnamesbool], "_csq")
setnames(csqs,colnames(csqs), newnames)
setkey(csqs, "variant_id","transcript_name")

cert = data.table(dcast(topranked, "variant_id + transcript_name  ~strain_name_1+strain_name_2", value.var=c("certainty")))
colnamesbool = !(colnames(cert) %in% c("variant_id" ,"transcript_name"))
newnames = colnames(cert)
newnames[colnamesbool] = paste0(newnames[colnamesbool], "_cert")
setnames(cert,colnames(cert), newnames)
setkey(cert, "variant_id","transcript_name")


merged = dips[csqs]
merged = merged[cert]
setkey(merged, "variant_id")

merged$NA_NA_diplo = NULL
merged$NA_NA_csq   = NULL
merged$NA_NA_cert  = NULL

write.table(merged,file.path(allfolder,"colformat.csv"),quote=F,row.names=F, sep="\t")
z = merged[,colnames(merged)[grepl("csq", colnames(merged))] ,with=F]
stoprows = c()
for (i in 1:ncol(z))
{
	stoprows =c(stoprows, which(grepl("stop", unname(unlist(z[,i, with=F])))))
}
stoprows = unique(stoprows)
kk = merged[stoprows]
write.table(kk,file.path(allfolder, "stops.csv"), row.names=F, sep="\t")

small = df
small$consequence_1 = NULL
small$consequence_2 = NULL
small$transcript_name = NULL
small = unique(small)

hets = small[ ,list(isSegregating = max(prob*(as.character(allele_1)!=as.character(allele_2)))),by="variant_id"]
