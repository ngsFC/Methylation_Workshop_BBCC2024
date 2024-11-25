#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ChAMP")

library(ChAMP)

setwd("YOUR_DIRECTORY")

type <- "EPIC"
norm <- "BMIQ"
testDir <- "YOUR_DIRECTORY"

myLoad <- champ.load(
  testDir,
  method="ChAMP",
  methValue="B",
  autoimpute=TRUE,
  filterDetP=TRUE,
  ProbeCutoff=0,
  SampleCutoff=0.1,
  detPcut=0.01,
  filterBeads=TRUE,
  beadCutoff=0.05,
  filterNoCG=TRUE,
  filterSNPs=TRUE,
  population=NULL,
  filterMultiHit=TRUE,
  filterXY=TRUE,
  force=FALSE,
  arraytype=type
)

CpG.GUI(CpG=rownames(myLoad$beta),
        arraytype=type)

QC.GUI(beta=myLoad$beta,
       arraytype=type)

#cores <- detectCores()

myNorm <- champ.norm(
  beta=myLoad$beta,
  arraytype=type, 
  method=norm,
  plotBMIQ = TRUE,
  #cores = cores-1
  cores = 1
)

QC.GUI(beta=myNorm,
       arraytype=type)

champ.SVD(beta=as.data.frame(myNorm), 
          pd=myLoad$pd)

myCombat <- champ.runCombat(beta=myNorm,
                            pd=myLoad$pd,
                            batchname=c("Slide"))

QC.GUI(beta=myCombat,
       arraytype=type)

myDMP <- champ.DMP(beta = myCombat,
                   pheno=myLoad$pd$Sample_Group,
                   arraytype = type)

DMP.GUI(DMP=myDMP[[1]],
        beta=myCombat,
        pheno=myLoad$pd$Sample_Group)

myDMR <- champ.DMR(beta=myCombat,
                   pheno=myLoad$pd$Sample_Group,
                   method="Bumphunter",
                   arraytype = type,
                   cores = cores-1,
                   B=50)

DMR.GUI(DMR=myDMR,
        beta=myCombat,
        pheno = myLoad$pd$Sample_Group,
        runDMP = T,
        arraytype = type)

myGSEA <- champ.GSEA(beta=myCombat,
                     DMP=myDMP[[1]], 
                     DMR=myDMR, 
                     arraytype=type,
                     adjPval=0.05, 
                     method="fisher",
                     pheno=myLoad$pd$Sample_Group,
                     cores=cores-1)

View(myGSEA$DMP)
