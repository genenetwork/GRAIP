start = 0  #starting number of permutation
end = 1000 #ending number of permuatation

library(qtl) #loads rQTL

  
for (i in start:end) {
    myDir=""  #replace with directory path appropriate to OS

    #genotype and phenotype file names.  Note that the genotype file name is variable
    myGenoFile= paste(i,".txt",sep="")
    myPhenoFile= "pheno.dat"

    #generating maps
    permdata <- read.cross("gary", dir=myDir, genfile=myGenoFile,phefile=myPhenoFile, chridfile="chrid.dat", mnamesfile="mnames.txt", pnamesfile="pnames.txt")
    permdata <- sim.geno(permdata, n.draws=16, step=0)
    summary.cross(permdata)
    out.em <- scanone(permdata, method=c("imp"), pheno.col=3)

 
    #the rest of the code captures output files

    #output file name   
    outPutFile = paste(myDir,"out",i,".txt",sep="")
 
    #capture output file
    capture.output(print(out.em), file=outPutFile,append=FALSE)

    closeAllConnections()
    }
}






