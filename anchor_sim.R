set.seed(1234)

sim_data<-function(N,g2,ni=30,offsets) {
    ##g1 mean and var are 0/1
    diff<-rnorm(ni)
    n<-length(offsets)
    ##simulate responses for first group
    rnorm(N,0,1)->th1
    th.mat<-matrix(th1,N,n,byrow=FALSE)
    diff.mat<-matrix(diff,N,n,byrow=TRUE)
    kern<-th.mat-diff.mat
    pv<-exp(kern)/(1+exp(kern))
    test<-matrix(runif(N*n),N,n)
    ifelse(pv>test,1,0)->resp1
    ##now group 2
    rnorm(N,g2[1],g2[2])->th2
    th.mat<-matrix(th2,N,n,byrow=FALSE)
    diff.mat<-matrix(diff,N,n,byrow=TRUE)
    for (i in 1:ncol(diff.mat)) diff.mat[,i]<-diff.mat[,i]+offsets[i]
    kern<-th.mat-diff.mat
    pv<-exp(kern)/(1+exp(kern))
    test<-matrix(runif(N*n),N,n)
    ifelse(pv>test,1,0)->resp2
    ##
    rbind(resp1,resp2)->resp
    paste("i",1:n,sep="")->colnames(resp)
    c(rep("g1",N),rep("g2",N))->gr
#####################################################
##stan
    J<-nrow(resp)
    K<-ncol(resp)
    N <- J*K;
                                        #
    y <- rep(-1,N);
    jj <- rep(-1,N);
    kk <- rep(-1,N);
    n <- 1;
    for (j in 1:J) {
        for (k in 1:K) {
            y[n] <- resp[j,k];
            jj[n] = j;
            kk[n] = k;
            n <- n + 1;
        }
    }
    L<-list(J=J,K=K,N=N,jj=jj,kk=kk,y=y,g=ifelse(gr=="g1",1,2))
    list(stan=L,true=list(th=c(th1,th2),g2=g2,diff=diff,offsets=offsets))
}

#######################################################################
##simulating data
library(rstan)

ni<-30




L<-list()
for (N in c(1000,2500)) {
    for (mu in c(0,-.25,-1)) {
        for (M in c(.1,.5,1)) {
            offsets<-c(0,runif(ni-1,-1*M,M)) #always set one to 0 just for fun.
            g2<-c(mu,.6)
            sim_data(N=N,g2=g2,offsets=offsets,ni=ni)->L[[paste(N,mu,M)]]
        }
    }
}

f<-function(L) {
    library(rstan)
    #fit <- stan(file='~/colombia/anchor/irt-ml.stan',data=L$stan,iter=2000,chains=3)
    fit <- stan(file='~/colombia/anchor/irt-ml.stan',data=L$stan,iter=2000,chains=2)
    list(L=L,fit=fit)
}

length(L)->np
min(30,np)->np
library(parallel)
makeCluster(np)->cl
clusterApply(cl,L,f)->out
stopCluster(cl)
names(L)->names(out)
save(out,file="sim-stan.Rdata")



load("sim-stan.Rdata")

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

fun<-function(nm,out) {
    out[[nm]]->ll
    library(rstan)
    tr<-list()
                                        #
    ll$fit->fit
    ll$L->L
    extract(fit)$mu_alpha->foo
    mean(foo[,1])->m1
    mean(foo[,2])->m2
    c(m1,m2,m2-m1,L$true$g2[1])->tr$means #estimated group 2 mean versus true group 2 mean
                                        #
    extract(fit)$alpha->M
    cor(L$true$th,colMeans(M))->tr$coor #should be high, this is row scores and theta correlatin
                                        #
    anch(fit)->anchor
    extract(fit)$beta->B
    colMeans(B[,,1])->b1
    colMeans(B[,,2])->b2
    c(sum(b1), sum(b2))->tr$sums #should be around 0, identification constraint
                                        #
    cbind(b1,b2)->tmp
    col<-ifelse(anchor,"red","black")
                                        #col<-c("red",rep("black",100))
    #layout(matrix(c(1,2,3,3),2,2,byrow=TRUE))
    #par(mgp=c(2,1,0),mar=c(3,3,1,1))
    #plot(L$true$diff,b1,pch=19); abline(0,1)
    #plot(L$true$diff+L$true$offsets,b1,pch=19); abline(0,1)
    plot(L$true$offsets,b2-b1,pch=19,xlim=c(-1,1),ylim=c(-1,1),col=col); abline(0,1,lty=2)
    mtext(side=3,nm,line=0)
    abline(v=0)
    abline(v=mean(L$true$offsets),col="red")
    abline(h=L$true$g2[1],col="blue")
    unlist(tr)
}
pdf("/tmp/pic.pdf")
lapply(names(out),fun,out)->zz
dev.off()

do.call("rbind",zz)->tab
strsplit(names(out)," ")->txt
par(mgp=c(2,1,0))
plot(tab[,4],tab[,3],type="n",xlab="True diff in means",ylab="Est diff in means")
abline(0,1)
text(tab[,4],tab[,3],sapply(txt,"[",3))


fun<-function(nm,out) {
    out[[nm]]->ll
    library(rstan)
    tr<-list()
                                        #
    strsplit(nm," ")->tr$base
    ll$fit->fit
    ll$L->L
    extract(fit)$mu_alpha->foo
    mean(foo[,1])->m1
    mean(foo[,2])->m2
    #m2-m1 -> tr$delta
    extract(fit)$beta->B
    anch(fit)->anchor
                                        #mean(b2-b1) -> tr$b.diff
                                        #mean(L$true$offsets[anchor])->tr$anchor
    sum(anchor)->tr$SA
    mean(abs(L$true$offsets[anchor]))->tr$anch
                                        #max(abs(L$true$offsets[anchor]))->tr$anch.max
    quantile(abs(L$true$offsets[anchor]),.9)->tr$anch.max
    mean(abs(L$true$offsets[!anchor]))->tr$not.anch
                                        #min(abs(L$true$offsets[!anchor]))->foo
    quantile(abs(L$true$offsets[!anchor]),.1)->foo
    ifelse(is.finite(as.numeric(foo)),foo,NA)->tr$not.min
    unlist(tr)
}
lapply(names(out),fun,out)->zz
do.call("rbind",zz)->tab

par(mgp=c(2,1,0),mar=c(3.2,10,2,1))
plot(tab[,5],1:nrow(tab),type="b",pch=19,xlim=c(0,.7),xlab="mean abs offset",ylab="",yaxt="n")
points(tab[,6],1:nrow(tab),type="b",pch=1,col="black",lty=2)
points(tab[,7],1:nrow(tab),type="b",pch=19,col="red")
points(tab[,8],1:nrow(tab),type="b",pch=1,col="red",lty=2)
mtext(adj=.5,side=2,line=8,tab[,1],at=1:nrow(tab),cex=.7,las=1)
mtext(adj=.5,side=2,line=6,tab[,2],at=1:nrow(tab),cex=.7,las=1)
mtext(adj=.5,side=2,line=4,tab[,3],at=1:nrow(tab),cex=.7,las=1)
mtext(adj=.5,side=2,line=2,round(as.numeric(tab[,4])/30,digits=2),at=1:nrow(tab),cex=.7,las=1)
mtext(adj=.5,side=2,line=8,at=0.5,"N/group",las=1,cex=.7)
mtext(adj=.5,side=2,line=6,at=nrow(tab)+.5,"mean diff",las=1,cex=.7)
mtext(adj=.5,side=2,line=4,at=0.5,"offsets",las=1,cex=.7)
mtext(adj=.5,side=2,line=2,at=nrow(tab)+.5,"% anchor",las=1,cex=.7)
for (i in 1:nrow(tab)) abline(h=i,lty=1,lwd=.4,cex=.7)
legend(x=-.3,y=0,c("anchor","variant","mean","extrema"),pch=c(19,19,NA,NA),col=c("black","red","black","black"),lty=c(NA,NA,1,2),xpd=NA,ncol=2,bty="n",cex=.8)

#######################################################################
##mirt

## ##base model, no adjustment for dtf
## library(mirt)
## models <- paste('F1 = 1-',ncol(resp)-1,sep="")
## mod.base <- multipleGroup(resp, models,itemtype="Rasch", group = gr, SE=TRUE,invariance=c("free_means","free_var","intercepts"))
## ##item params are constant, for example:
## coef(mod.base)[[1]][28]
## coef(mod.base)[[2]][28]
## ##and we can see that the group mean for g2 is off
## coef(mod.base)[[1]][n+1]
## coef(mod.base)[[2]][n+1] #should only be -.3

## ##separate model, can't compare group means
## library(mirt)
## models <- paste('F1 = 1-',ncol(resp),sep="")
## mod.sep <- multipleGroup(resp, models,itemtype="Rasch", group = gr, SE=TRUE)

## coef(mod.sep)->co
## pars<-list()
## for (i in 1:2) { ##get coefficients
##     co[[i]]->tmp
##     length(tmp)->j
##     print(tmp[[j]]) #these are the mean/variance components, notice that g1/g2 are 0/1 here.
##     tmp[-j]->tmp
##     f<-function(x) x[1,]
##     lapply(tmp,f)->tmp
##     do.call("rbind",tmp)[,2]->pars[[i]]
## }
## do.call("cbind",pars)->pars
## plot(pars); abline(0,1) #items are just a little easier for group 1 than group 2. 
## mean(pars[,2]-pars[,1]) #should just be -.2 given way offsets are defined, but is now -.2 + -.3 (mean offset plus latent mean diff)

## ## est<-list()
## ## for (j in 1:ncol(resp)) {
## ##     print(j)
## ##     est[[j]]<-multipleGroup(resp,group=gr,method="EM",itemtype="Rasch",models,invariance=c("free_means","free_var",colnames(resp)[j]),verbose=FALSE)
## ## }
## ## anovas <- lapply(est, anova, object2=mod.sep, verbose=FALSE)!
## ## sapply(anovas,function(x) x$p[2])->p
## ## p

## DIF(mod.sep, which.par = c( 'd'), Wald=FALSE, scheme="drop_sequential",p.adjust = 'fdr')->dif

## $adj_pvals->pv.r
## plot(offsets,log10(pv.r),pch=19,col=c("red",rep("black",100)))
## abline(h=log10(0.05))

## plot(offsets,pv.r,pch=19,col=c("red",rep("black",100)))












## ##for fun let's look at j=2
## coef(est[[5]][[1]][1:3])
## coef(est[[5]][[2]][1:3]) ##note that the second param are held constant





## a <- matrix(rep(15,1), ncol=1)
## d <- matrix(rnorm(15,0,.7),ncol=1)
## itemtype <- rep('dich', nrow(a))
## N <- 500
## dataset1 <- simdata(a, d, N, itemtype)
## dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
## dat <- rbind(dataset1, dataset2)
## group <- c(rep('D1', N), rep('D2', N))
## models <- 'F1 = 1-15'

## ##separate model, can't compare group means
## library(mirt)
## models <- paste('F1 = 1-',ncol(dat),sep="")
## mod.sep <- multipleGroup(dat, models,itemtype="Rasch", group = gr, SE=TRUE)

## coef(mod.sep)->co
## pars<-list()
## for (i in 1:2) { ##get coefficients
##     co[[i]]->tmp
##     length(tmp)->j
##     print(tmp[[j]]) #these are the mean/variance components, notice that g1/g2 are 0/1 here.
##     tmp[-j]->tmp
##     f<-function(x) x[1,]
##     lapply(tmp,f)->tmp
##     do.call("rbind",tmp)[,2]->pars[[i]]
## }
## do.call("cbind",pars)->pars
## plot(pars); abline(0,1) #items are just a little easier for group 1 than group 2. 
## mean(pars[,2]-pars[,1]) #should just be -.2 given way offsets are defined, but is now -.2 + -.3 (mean offset plus latent mean diff)
