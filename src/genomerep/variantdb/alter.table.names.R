library("RMySQL")
user     = "root"
password = "FrozenMixedVegetables"
db       = "rachel_x"
host     = "127.0.0.1"
##an unused port. Someday, write code to check which ports are available (i.e. rotation student)
port     = 3308

##either get system call to work without hanging, or do this manually on the CL to open an SSH tunnel
##from port 3308 on this machine to port 3306 on the valdardb machine. Then when rmysql connects to port 3308 locally, it sees valdardb
sshcommand = paste0("ssh -fnNT -L ",port,":localhost:3306 doreper@valdardb.its.unc.edu")
print(sshcommand)
system(sshcommand)

con = dbConnect(MySQL(), user=user, port = port, password =password, dbname=db, host=host)
dbSendQuery(con, 'set @@max_heap_table_size=4294967296;')

## renameMap = list(founder_allele = "founder_genotype",
##                  inbred_line_allele = "genotype",
##                  inbred_line_allele_hap = "allele_sampling",
##                  inbred_line_fixed      = "fixed",
##                  inbred_line_founder    = "diplotype",
##                  inbred_line_founder_hap= "haplo_sampling")

renameMap = list(
                 varinfo= "consequence_info")

for(nam in names(renameMap))
{
    query = paste0("RENAME table ", nam, " TO ", renameMap[[nam]])
    dbSendQuery(con, query)
}
dbCommit(con)

out = create.F1.table(con, strain1, strain2)
write.table(out, file.path(mnp$data, mnp$NOD.B6.variantsFile), sep="\t", row.names=F)
