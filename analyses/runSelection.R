library(dplyr)
library(moiR)
library(data.table)
source("~/safedata/natvar/analyses/phenoselection_multi_FUNCTIONS-copy.R")
## load datasets
pheno <- read.table(file = '../data/atlas1001_phenotypes_matrix.csv', sep=",", header = T)
pheno[1:5,1:5]
load("../data/phenotypenames.rda")

load("../data/d4.rda")
idsfield<-unique(d4$id)
whichfield<-which(pheno$V1 %in% idsfield)
idex515<-whichfield

## removed duplicated phenotype from Martinez-berdeja and Kalladan
pheno <- pheno[,c(1:1706, 1726:1883)]
pheno<- as.data.frame(pheno)
phenotypenames <- phenotypenames[c(1:1705,1725:1882 ), ]
p<- pheno[idex515,-1]
colnames(p) <- paste(phenotypenames[,2])
p<- cbind(id=pheno[idex515,1], p)

## add in variance in Flowering time
FTvariance_df <- readRDS(file="../data/FTvariance_df.rda")
p[,ncol(p)+1] <- FTvariance_df$FT_variance
dim(p)

## fix -9 to NAs
id <- grep("FT_", colnames(p))[1:8]
p[,id][p[, id] == -9]  <- NA
p[,id]


## seperate fitness data from Exposito-Alonso field experiment and other phenotype data
fitness <- colnames(p)[grep("Fitness", colnames(p))][1:8]
survival <- colnames(p)[grep("Survival", colnames(p))][1:8]
seeds <- colnames(p)[grep("rSeeds", colnames(p))][1:8]
all_fitness <- p[,c(fitness, survival, seeds)]
justphenos <- p[,!colnames(p) %in% c(fitness, survival, seeds)]
justphenos <- justphenos[,-1]

all_fitness <- p[,c(fitness, survival, seeds)]
for (k in colnames(all_fitness)){
  all_fitness[,k][all_fitness[, k] == -9] <- NA
}

sum(all_fitness==-9, na.rm=T)
#dim(justphenos)
justphenos <- apply(justphenos, 2, fn)
justphenos <- apply(justphenos, 2, scale)
justphenos[1:5,1:5]
print("Got here")

lsel<-c()
for(i in c("rFitness","rSurvival_fruit","rSeeds")){
  print(i)
  for(j in c("mlp","mli","mhi","mhp", "tlp", "tli","thi","thp")){
    trait<-paste0(i,"_",j)
    print(trait)
    w=all_fitness[,trait]/mean(all_fitness[,trait], na.rm=T)
    z=justphenos
    d1<-preparedata(w,z)
    res<-PHENOSELECTION(Variables=z, Fitness=w, replicates = 50)
    # tmp1<-do.call(cbind,parseformatted(res$gradient_linear)) %>% data.frame
    # tmp1$param<-"beta"
    tmp2<-do.call(cbind,parseformatted(res$coefficient_linear)) %>% data.frame
    tmp2$param<-"s"
    print(tmp2)
    sel<-tmp2
    sel$mean<-fn(sel$mean)
    sel$se<-fn(sel$se)
    sel$lower=sel$mean - sel$se *1.96
    sel$upper=sel$mean + sel$se *1.96
    sel$trait<- colnames(justphenos)
    sel$fitness<-i
    sel$env<-j
    lsel<- rbind(lsel, sel)
  }
}
saveRDS(lsel, file="../data/SelectionRepResults_AllPhenos_083022.rda")
