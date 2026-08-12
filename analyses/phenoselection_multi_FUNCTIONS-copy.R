#### Helper function ####

parseformatted<-function(x){
  tmp<-lapply(x,function(i) strsplit(i,split = " (",fixed=T)[[1]][1])
  mymean<-as.numeric(unlist(tmp))
  tmp2<-lapply(x,function(i) strsplit(i,split = " (",fixed=T)[[1]][2])
  tmp3<-lapply(tmp2,function(i) strsplit(i,split = ")",fixed=T)[[1]][1])
  myse<-as.numeric(unlist(tmp3))
  mystars<-lapply(tmp2,function(i) strsplit(i,split = ")",fixed=T)[[1]][2])
  mystars<-as.character(unlist(mystars))
  mystars[is.na(mystars)]<-"ns"
  return(list(mean=mymean,se=myse,signi=mystars))
}

preparedata<-function(Fitness,Variables){
  # prepare data function
  w=as.numeric(Fitness)
  w=w/mean(w,na.rm=T)

  phenos<-apply(Variables,2,function(x) scale(as.numeric(x)))
  colnames(phenos)<-paste0("Var",1:ncol(Variables))

  df<-data.frame(w,phenos)
  return(df)
}

makeselectiontable<-function(obj){
  tmp1<-rbind(obj$gradient_linear,
        obj$coefficient_linear)
  # tmp3<-cbind(tmp1,tmp2)
  row.names(tmp1)<-c("B","s")
  colnames(tmp1)<-paste0("Var",1:ncol(tmp1))

  if(length(obj$gradient_quadratic)==42){
    tmp2<-rbind(c(
                obj$gradient_quadratic[1],
                obj$gradient_quadratic[4],
                obj$gradient_quadratic[3]) ,
                c(
                  obj$coefficient_quadratic[1],
                  obj$coefficient_quadratic[4],
                  obj$coefficient_quadratic[3])
                )

    row.names(tmp2)<-c("g","c")
    colnames(tmp2)<-c("Var1","Var2","Var1:2")

  }else{
    tmp4<-matrix(obj$coefficient_quadratic,
                 ncol=length(obj$gradient_linear),
                 nrow=length(obj$gradient_linear),byrow = T)
    tmp5<-matrix(obj$gradient_quadratic,
                 ncol=length(obj$gradient_linear),
                 nrow=length(obj$gradient_linear),byrow = T)
    tmp2<-tmp4
    tmp2[lower.tri(tmp2)]<-tmp5[lower.tri(tmp5)]
    colnames(tmp2)<-paste0("cVar",1:length(obj$gradient_linear))
    rownames(tmp2)<-paste0("gVar",1:length(obj$gradient_linear))
  }
  return(list(linear=tmp1,quadratic=tmp2))
}


#### Selection functions ####

gradientlinear<- function(data, indices) {
  d1 <- data[indices,]
  P<-cov(d1[,-1], use="pairwise.complete.obs")
  Pinv<-solve(P)
  s<- cov(d1, use="pairwise.complete.obs")[-1,1]
  B<-Pinv%*%s
  return(B)
}

coeflinear<- function(data, indices) {
  d1 <- data[indices,]
  #P<-cov(d1[,-1])
  #Pinv<-solve(P)
  s<- cov(d1, use="pairwise.complete.obs")[-1,1]
  #B<-Pinv%*%s
  return(s)
}

gradientquadratic<- function(data, indices) {
  d1 <- data[indices,]
  w<-d1[,1]
  pheno<-d1[,-1]
  P<-cov(pheno)
  Pinv<-solve(P)
  n=ncol(pheno)
  c<-matrix(ncol=n,nrow=n)
  for(i in 1:n){
    for(j in 1:n){
      new<-pheno[,i]*pheno[,j]
      c[i,j]<-cov(w, new)
    }
  }
  gamma = Pinv%*% c %*%Pinv
  return(gamma)
}

coefquadratic<- function(data, indices) {
  d1 <- data[indices,]
  w<-d1[,1]
  pheno<-d1[,-1]
  P<-cov(pheno)
  Pinv<-solve(P)
  n=ncol(pheno)
  c<-matrix(ncol=n,nrow=n)
  for(i in 1:n){
    for(j in 1:n){
      new<-pheno[,i]*pheno[,j]
      c[i,j]<-cov(w, new)
    }
  }
  # new<-d1[,-1]^2
  # z12<-d1[,2]*d1[,3] # problem, need to figure out
  # new<-data.frame(new,z12)
  # c<-matrix(ncol=2,nrow=2)
  # rawcov<-cov(d1$w,new)
  # c<-diag(cov(d1$w,new)[1:2])
  # c[1,2]<-rawcov[3]
  # c[2,1]<-rawcov[3]
  gamma = Pinv%*% c %*%Pinv
  return(gamma)
}

responselinear<- function(data, indices) {
  d1 <- data[indices,]
  P<-cov(d1[,-1])
  Pinv<-solve(cov(d1[,-1]))
  s<- cov(d1)[-1,1]
  B<-Pinv%*%s
  deltaZ<-Gmatrix%*%B
  deltaZ
  newZ<-deltaZ
  meanchange<-newZ
  return(meanchange)
}

responsequadratic<- function(data, indices) {
  d1 <- data[indices,]
  P<-cov(d1[,-1])
  Pinv<-solve(cov(d1[,-1]))

  d1 <- data[indices,]
  P<-cov(d1[,-1])
  Pinv<-solve(cov(d1[,-1]))
  s<- cov(d1)[-1,1]
  B<-Pinv%*%s


  new<-d1[,2:3]^2
  z12<-d1[,2]*d1[,3]
  new<-data.frame(new,z12)
  c<-matrix(ncol=2,nrow=2)
  rawcov<-cov(d1$w,new)
  c<-diag(cov(d1$w,new)[1:2])
  c[1,2]<-rawcov[3]
  c[2,1]<-rawcov[3]
  gamma = Pinv%*% c %*%Pinv


  deltaP<- (P*gamma*P) - s%*%t(s)


  deltaG<-Gmatrix%*%(gamma-(B%*%t(B)))%*%Gmatrix
  deltaG+Gmatrix

  return(deltaP)

}

responsequadratic_gmatrix<- function(data, indices) {

  d1 <- data[indices,]
  P<-cov(d1[,2:3])
  Pinv<-solve(cov(d1[,2:3]))
  s<- cov(d1)[2:3,1]
  B<-Pinv%*%s


  new<-d1[,2:3]^2
  z12<-d1[,2]*d1[,3]
  new<-data.frame(new,z12)
  c<-matrix(ncol=2,nrow=2)
  rawcov<-cov(d1$w,new)
  c<-diag(cov(d1$w,new)[1:2])
  c[1,2]<-rawcov[3]
  c[2,1]<-rawcov[3]
  gamma = Pinv%*% c %*%Pinv
  return(gamma)


  deltaP<- (P*gamma*P) - s%*%t(s)  # why variance is not 1??
  deltaP
  newP<-deltaP+P
  covariancechange<-deltaP

  deltaG<-Gmatrix%*%(gamma-(B%*%t(B)))%*%Gmatrix
  deltaG+Gmatrix # the change in heritability and genetic correlation!

  return(deltaG)

}

#### Bootstrap function ####
#x= abs(rnorm(n=100, mean = .5, sd = .5))

treatbootstrap<-function(x){
    five<-quantile(x,p=c(0.05,0.95) ,na.rm=T)
    one<-quantile(x,p=c(0.01,0.99) ,na.rm=T)
    zeroone<-quantile(x,p=c(0.001,0.999) ,na.rm=T)
    tempsign<-(five[1]/five[2]) / abs(five[1]/five[2])
    tempsign<-c(tempsign,(one[1]/one[2]) / abs(one[1]/one[2]))
    tempsign<-c(tempsign,(zeroone[1]/zeroone[2]) / abs(zeroone[1]/zeroone[2]))
    tempsign[tempsign==1]<-"*"
    tempsign[tempsign==-1]<-""
    tempsign[is.na(tempsign)]<-""
    tempsign<-paste(as.character(tempsign)[1],as.character(tempsign)[2],as.character(tempsign)[3],sep="")
    as.character(tempsign)
    se<-round(sqrt(var(x)),digits = 3)
    media<-round(mean(x),digits = 3)
    pasted<-paste(media," (",se, ")",tempsign,sep="")
    return(pasted)
  }

extractbootstrap<-function(bootstrapresults){
  extracted<-apply(bootstrapresults$t,2,treatbootstrap)
  return(extracted)
}

extractbootstrap_numeric<-function(bootstrapresults){
  treatbootstrap<-function(x){
    media<-round(mean(x),digits = 3)
    return(media)
  }
  extracted<-apply(bootstrapresults$t,2,treatbootstrap)
  return(extracted)
}

#### Master function ####

PHENOSELECTION<-function(Variables,Fitness,Gmatrix=NULL,replicates=100){
  require(boot)
  #@ start LITTLE BIT THAT ACTUALLY DO ANALYSES @#
  d1<-preparedata(Fitness,Variables)
  if (is.null(Gmatrix)==T){
    # print ("Heritabilities not provided, only selection analyses reported" )

    result_gradient_linear<- boot(data=d1, statistic=gradientlinear, R=replicates)
    #result_gradient_quadratic<- boot(data=d1, statistic=gradientquadratic, R=replicates)
    result_coefficient_linear<- boot(data=d1, statistic=coeflinear, R=replicates)
    #result_coefficient_quadratic<- boot(data=d1, statistic=coefquadratic, R=replicates)
    resa<-extractbootstrap(result_gradient_linear)
    #resb<-extractbootstrap(result_gradient_quadratic)
    resc<-extractbootstrap(result_coefficient_linear)
    #resd<-extractbootstrap(result_coefficient_quadratic)
    # analysislist<-list(gradient_linear=resa,gradient_quadratic=resb,
    #                    coefficient_linear=resc,coefficient_quadratic=resd)
    analysislist<-list(gradient_linear=resa,coefficient_linear=resc)

  }
  if (is.null(Gmatrix)==F){
    # print (" Heritabilities provided, selection gradients and response to selection analyses reported" )

    result_gradient_linear<- boot(data=d1, statistic=gradientlinear, R=replicates)
    #result_gradient_quadratic<- boot(data=d1, statistic=gradientquadratic, R=replicates)
    result_coefficient_linear<- boot(data=d1, statistic=coeflinear, R=replicates)
    #result_coefficient_quadratic<- boot(data=d1, statistic=coefquadratic, R=replicates)

    resa<-extractbootstrap(result_gradient_linear)
    #resb<-extractbootstrap(result_gradient_quadratic)
    resc<-extractbootstrap(result_coefficient_linear)
    #resd<-extractbootstrap(result_coefficient_quadratic)

    result_response_linear<- boot(data=d1, statistic=responselinear, R=replicates)
    #result_response_quadratic<- boot(data=d1, statistic=responsequadratic, R=replicates)
    #result_response_quadratic_gmatrix<- boot(data=d1, statistic=responsequadratic_gmatrix, R=replicates)

     rese<-extractbootstrap(result_response_linear)
    # resf<-extractbootstrap(result_response_quadratic)
    # resg<-extractbootstrap(result_response_quadratic_gmatrix)

    # analysislist<-list(gradient_linear=resa,gradient_quadratic=resb,
    #                    coefficient_linear=resc,coefficient_quadratic=resd,
    #                    response_linear=rese,
    #                    response_quadratic_Vpheno=resf,
    #                    response_quadratic_Vaddit=resg)
    analysislist<-list(gradient_linear=resa,coefficient_linear=resc,response_linear=rese)

  }
  #@ end LITTLE BIT THAT ACTUALLY DO ANALYSES @#
  return(analysislist)
}


#### C function ####

# library(Rcpp)
# sourceCpp("src/phenotypeselection.cpp")

PHENOSELECTIONcwrap<-function( w=rnorm(100),
                               z=cbind(rnorm(100),rnorm(100)),
                               g= diag(ncol(z)),
                               replicates=10,
                               debug = FALSE){
  # res<-PHENOSELECTIONc(w,z,g,replicates)
  # metricsnames<-names(res[[1]])
  # pars<- lapply(1:length(metricsnames),function(x){
  #               tmp<-purrr::map(1:replicates,function(i) unlist(res[[i]][x])) %>% do.call(rbind,.)
  #              })
  # names(pars) = metricsnames
  pars<-PHENOSELECTIONcboot(w,z,g,replicates,debug)
  return(pars)
}


meanse_<-function(x){
  m<-apply(x,2,mean)
  lower<-apply(x,2,function(x)quantile(x,p=c(0.025) ))
  upper<-apply(x,2,function(x)quantile(x,p=c(0.975) ))
  pseudop<-apply(x,2,function(x){
    if(mean(x)>0){
      pval=sum(x<0) / length(x)
    }else{
      pval=sum(x>0) / length(x)
    }
    return(pval)
  })
  se<-apply(x,2,function(i) sd(i) ) # in bootstrap se=sd not sd/sqrt n
  lowerwse<-m - se * 1.96
  upperwse<-m + se * 1.96
  pse<-pnorm((m/se), mean = 0, sd = 1, lower.tail = F)
  toreport<-rbind(m,lower,upper,pseudop,
                  se,lowerwse,upperwse,pse)
  row.names(toreport)<-c("estimate","lower","upper","pseudop",
                         "se","lowerwse","uppwerwse","pse")
  return(toreport)
}
meanse<-function(x){
  if(is.list(x)){
    tmp<-lapply(1:length(x),function(i){
      meanse_(x[[i]])
      })
    names(tmp)=names(pars)
  }else if(is.matrix(x)){
    tmp<-meanse_(x)
  }else{
    stop("Not a list nor a matrix")
  }
    return(tmp)
}
