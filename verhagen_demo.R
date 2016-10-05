read.table("baseLPG.txt")->x

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

L<-list(J=J,K=K,N=length(y),jj=jj,kk=kk,y=y,g=gr)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fit <- stan(file='irt-ml.stan',data=L,iter=2000,chains=3)

save(fit,file="verhagen-fit.Rdata")





library(rstan)

                                        #means
extract(fit)$mu_alpha->foo
mean(foo[,1])->m1
mean(foo[,2])->m2


extract(fit)$beta->B
colMeans(B[,,1])->b1
colMeans(B[,,2])->b2
cbind(b1,b2)->tmp


f<-function(x) {
    apply(x,2,quantile,c(.025,.975))->qu
}
f(B[,,1])->ci1
f(B[,,2])->ci2
rbind(ci1,ci2)->ends
t(ends)->ends

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
anch(fit)->anchor

range(ends)->xl
par(mgp=c(2,1,0),mar=c(3.2,3,2,2.5))
plot(NULL,xlim=xl,ylim=c(0,1+nrow(ends)),xlab="parameter estimates",ylab="",yaxt="n")
mtext(side=2,1:nrow(ends),at=1:nrow(ends),las=1,line=.2)
for (i in 1:nrow(ends)) {
    ifelse(anchor[i],1,3)->lty
    segments(ends[i,1],i-.1,ends[i,2],i-.1,col="blue",lwd=2,lty=lty)
    segments(ends[i,3],i+.1,ends[i,4],i+.1,col="black",lty=lty,lwd=2)
    text(b1[i],i-.4,round(b1[i],digits=2),cex=.7)
    text(b2[i],i+.4,round(b2[i],digits=2),cex=.7)
}
##
compare<-function(b1,b2) {
    mean(b1)->m1
    mean(b2)->m2
    sd(b1)->s1
    sd(b2)->s2
    abs(m1-m2)/sqrt(s1^2+s2^2)
}
comp<-numeric()
for (i in 1:11) compare(B[,i,1],B[,i,2])->comp[i]
mtext(side=4,line=.2,las=1,round(comp,digits=2),at=1:nrow(ends))
mtext(side=3,line=.2,"Replication of results from Verhagent et al. Table 2")
legend("topright",c("males","females","anchor","variant"),lty=c(NA,NA,1,3),pch=c(19,19,NA,NA),lwd=2,col=c("blue","black","black","black"),bty="n")
