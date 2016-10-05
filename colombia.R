set.seed(10123010)
setwd("/home/domingue/colombia/anchor")
fun<-function(x,nm) {
    x[,1]->gr
    x[,-1]->resp
    J<-nrow(resp)
    K<-ncol(resp)
    y<-list()
    jj<-list()
    kk<-list()
    for (i in 1:ncol(resp)) {
        resp[,i]->y[[i]]
        1:nrow(resp)->jj[[i]]
        rep(i,nrow(resp))->kk[[i]]
    }
    do.call("c",y)->y
    do.call("c",jj)->jj
    do.call("c",kk)->kk
    #
    L<-list(J=J,K=K,N=length(y),jj=jj,kk=kk,y=y,g=gr)
                                        #
    library(rstan)
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
                                        #
    fit <- stan(file='irt-ml.stan',data=L,iter=4000,chains=5)
    save(fit,file=paste("colfit-",nm,".Rdata",sep=""))
}

for (nm in c("QR_String","CR_String")) {
    load(paste("/home/domingue/colombia/community_colleges/",nm,".Rdata",sep=""))
    rowSums(is.na(tmp))==0->test
    tmp[test,]->tmp
    sample(1:nrow(tmp),10000)->index
    tmp[index,]->x
    fun(x,nm)
}



anch<-function(fit) {
    extract(fit)$beta->B
    colMeans(B[,,1])->b1
    colMeans(B[,,2])->b2
    f<-function(x) {
        apply(x,2,quantile,c(.025,.975))->qu
    }
    f(B[,,1])->ci1
    f(B[,,2])->ci2
    overlap<-function(ci1,ci2) {
        ci1[1] >= ci2[1] & ci1[1]<=ci2[2] -> test1
        ci1[2] >= ci2[1] & ci1[2]<=ci2[2] -> test2
        ci2[1] >= ci1[1] & ci2[1]<=ci1[2] -> test3
        ci2[2] >= ci1[1] & ci2[2]<=ci1[2] -> test4
        (test1 | test2) | (test3 | test4)
    }
    anchor<-logical()
    for (i in 1:ncol(ci1)) overlap(ci1[,i],ci2[,i])->anchor[i]
    anchor
}

library(rstan)
anch.out<-list()
for (nm in c("QR_String","CR_String")) {
    load(paste("colfit-",nm,".Rdata",sep=""))
    anch(fit)->anch.out[[nm]]
}
