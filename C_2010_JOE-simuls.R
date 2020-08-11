################################################################################
## R CODE FOR CATTANEO (2009)
## Simulation Study for "Efficient Semiparametric Estimation of Multi-valued Treatment Effects under Exogeneity"
## DATE: 30-Jun-2009
################################################################################
# cd /mainfiles/cattaneo/research/EfficientMultitreat/Simulations

rm(list=ls(all=TRUE))

library(nnet)

################################################################################
## Simulation Setup
################################################################################
S=5000
n=1000
J=3
do.IPW = TRUE; do.EIF = TRUE
x.gen=function(x){return(cbind(rep(1,nrow(x)),x,((x[,1])^2-1)/sqrt(2),x[,1]*x[,2]))}
B0=cbind(rep(0,5),c(1,1,1,1,1),c(2,2,2,2,2))

## TABLES
Table.colnames=c("T0-mean","T0-median","T1-mean","T1-median","T2-mean","T2-median")
Table.model1.IPW = matrix(NA,nrow=S,ncol=length(Table.colnames),dimnames=list(NULL,Table.colnames))
Table.model1.EIF = matrix(NA,nrow=S,ncol=length(Table.colnames),dimnames=list(NULL,Table.colnames))

Table.colnames=c("T0-loc","T1-loc","T2-loc")
Table.model2.IPW = matrix(NA,nrow=S,ncol=length(Table.colnames),dimnames=list(NULL,Table.colnames))
Table.model2.EIF = matrix(NA,nrow=S,ncol=length(Table.colnames),dimnames=list(NULL,Table.colnames))

Table.colnames=c("T0-mean","T0-median","ChangeT-mean","ChangeT-median")
Table.model3.IPW = matrix(NA,nrow=S,ncol=length(Table.colnames),dimnames=list(NULL,Table.colnames))
Table.model3.EIF = matrix(NA,nrow=S,ncol=length(Table.colnames),dimnames=list(NULL,Table.colnames))

Table.colnames=c("T0-loc","ChangeT-loc")
Table.model4.IPW = matrix(NA,nrow=S,ncol=length(Table.colnames),dimnames=list(NULL,Table.colnames))
Table.model4.EIF = matrix(NA,nrow=S,ncol=length(Table.colnames),dimnames=list(NULL,Table.colnames))


################################################################################
## Simulation Study
################################################################################
# s=1;
showwhen = 1; showevery=10
for (s in 1:S) {
    if (s==showwhen) {cat(paste("\nSimulations Completed:",s-1,"of",S,"-", Sys.time())); showwhen=showwhen+showevery}
    #x=x.gen(matrix(rnorm(2*n,0,1),n,2))
    x=x.gen(matrix(runif(2*n,-.5,.5),n,2))
    #y.choice=x%*%B0+cbind(rep(0,n),matrix((-log(-log(runif(2*n)))-EM.const)*sqrt(6/pi^2)*(1/2)^(1/4),n,2))
    y.choice=x%*%B0+matrix((-log(-log(runif(3*n)))),n,3)
    T=apply(y.choice,1,which.max)
    
    ## Generalized Propensity Score (p)
    sink("/dev/null")
    gps.fit = multinom(T ~ x-1, maxit=5000, simplified=T)
    sink()
    #summary(gps.fit)
    gps.hat = gps.fit$fitted.values
    T_gps.hat = cbind((T==1)/gps.hat[,1],(T==2)/gps.hat[,2],(T==3)/gps.hat[,3])
    
    xb = x%*%c(0,1/2,1/2,0,0); x=x[,c(2,3)]
    #u = rnorm(3*n,0,1)
    u = matrix(runif(3*n),n,3); u = -sign(u-.5)*log(1-2*abs(u-.5))
    y0 = xb + u[,1]; y1 = 1 + xb + u[,2]; y2 = 2 + xb + u[,3]
    y = (T==1)*y0 + (T==2)*y1 + (T==3)*y2
    
    ################################################################################
    ## Estimation
    ################################################################################
    ## MODEL 1
    e.IPW.m1 = matrix(NA, n, 6); e.EIF.m1 = matrix(NA, n, 6)
    score.IPW.m1 = matrix(NA,nrow=n,ncol=6); score.EIF.m1 = matrix(NA,nrow=n,ncol=6)

    for (t in 1:3) {
        pos=2*t-1
        ## IPW Estimator - Mean
        if (do.IPW) {
        tmp = function(b) return(abs(sum((T_gps.hat[,t])*(y-b))))
        IPW.m1 = optimize(f=tmp,interval=c(-4,6))
        Table.model1.IPW[s,pos] = IPW.m1$min
        e.IPW.m1[,pos] = x%*%lm((y[T==t]-IPW.m1$min)~x[T==t,]-1)$coeff
        score.IPW.m1[,pos] = (T_gps.hat[,t])*(y-IPW.m1$min)-e.IPW.m1[,pos]*(T_gps.hat[,t]-1)
        }
        ## EIF Estimator - Mean
        if (do.EIF) {
        tmp = function(b) return(abs(sum((T_gps.hat[,t])*(y-b)-x%*%lm((y[T==t]-b)~x[T==t,]-1)$coeff*(T_gps.hat[,t]-1))))
        EIF.m1 = optimize(f=tmp,interval=c(-4,6))
        Table.model1.EIF[s,pos] = EIF.m1$min
        e.EIF.m1[,pos] = x%*%lm((y[T==t]-EIF.m1$min)~x[T==t,]-1)$coeff
        score.EIF.m1[,pos] = (T_gps.hat[,t])*(y-EIF.m1$min)-e.EIF.m1[,pos]*(T_gps.hat[,t]-1)
        }

        pos=2*t
        ## IPW Estimator - Median
        if (do.IPW) {
        tmp = function(b) return(abs(sum((T_gps.hat[,t])*((y<b)-.5))))
        IPW.m1 = optimize(f=tmp,interval=c(-4,6))
        Table.model1.IPW[s,pos] = IPW.m1$min
        e.IPW.m1[,pos] = x%*%lm(((y[T==t]<IPW.m1$min)-.5)~x[T==t,]-1)$coeff
        score.IPW.m1[,pos] = (T_gps.hat[,t])*((y<IPW.m1)-.5)-e.IPW.m1[,pos]*(T_gps.hat[,t]-1)
        }
        ## EIF Estimator - Median
        if (do.EIF) {
        tmp = function(b) return(abs(sum((T_gps.hat[,t])*((y<b)-.5)-x%*%lm(((y[T==t]<b)-.5)~x[T==t,]-1)$coeff*(T_gps.hat[,t]-1))))
        EIF.m1 = optimize(f=tmp,interval=c(-4,6))
        Table.model1.EIF[s,pos] = EIF.m1$min
        e.EIF.m1[,pos] = x%*%lm(((y[T==t]<EIF.m1$min)-.5)~x[T==t,]-1)$coeff
        score.EIF.m1[,pos] = (T_gps.hat[,t])*((y<IPW.m1)-.5)-e.EIF.m1[,pos]*(T_gps.hat[,t]-1)
        }
    }
    if (do.IPW)V.IPW.m1 = crossprod(score.IPW.m1)/n
    if (do.EIF) V.EIF.m1 = crossprod(score.EIF.m1)/n

    ## MODEL 2
    for (t in 1:3) {
        pos=2*t-1
        ## IPW Estimator - Overidentification
        if (do.IPW) {
        W = solve(V.IPW.m1[c(pos,pos+1),c(pos,pos+1)])
        tmp = function(b) {
            moms = colSums(cbind((T_gps.hat[,t])*(y-b),(T_gps.hat[,t])*((y<b)-.5)))
            return(crossprod(moms,W)%*%moms)
        }
        IPW.m2 = optimize(f=tmp,interval=c(-4,6))
        Table.model2.IPW[s,t] = IPW.m2$min
        }

        if (do.EIF) {
        ## EIF Estimator - Overidentification
        W = solve(V.EIF.m1[c(pos,pos+1),c(pos,pos+1)])
        tmp = function(b) {
            moms = colSums(cbind((T_gps.hat[,t])*(y-b)-x%*%lm((y[T==t]-b)~x[T==t,]-1)$coeff*(T_gps.hat[,t]-1),
                                 (T_gps.hat[,t])*((y<b)-.5)-x%*%lm(((y[T==t]<b)-.5)~x[T==t,]-1)$coeff*(T_gps.hat[,t]-1)))
            return(crossprod(moms,W)%*%moms)
        }
        EIF.m2 = optimize(f=tmp,interval=c(-4,6))
        Table.model2.EIF[s,t] = EIF.m2$min
        }
    }

    ## MODEL 3
    ## IPW Estimator - Cross-equation Restriction
    if (do.IPW) {
    W = solve(V.IPW.m1[c(1,3,5),c(1,3,5)])
    tmp = function(b) {
        moms = colSums(cbind(T_gps.hat[,1]*(y-b[1]),
                             T_gps.hat[,2]*(y-(b[1]+b[2])),
                             T_gps.hat[,3]*(y-(b[1]+2*b[2]))))
        return(crossprod(moms,W)%*%moms)
    }
    IPW.m3 = optim(par=c(.5,.5),fn=tmp,method="Nelder-Mead")
    Table.model3.IPW[s,c(1,3)] = IPW.m3$par

    W = solve(V.IPW.m1[c(2,4,6),c(2,4,6)])
    tmp = function(b) {
        moms = colSums(cbind(T_gps.hat[,1]*((y<b[1])-.5),
                             T_gps.hat[,2]*((y<(b[1]+b[2]))-.5),
                             T_gps.hat[,3]*((y<(b[1]+2*b[2]))-.5)))
        return(crossprod(moms,W)%*%moms)
    }
    IPW.m3 = optim(par=c(.5,.5),fn=tmp,method="Nelder-Mead")
    Table.model3.IPW[s,c(2,4)] = IPW.m3$par
    }

    ## EIF Estimator - Cross-equation Restriction
    if (do.EIF) {
    W = solve(V.EIF.m1[c(1,3,5),c(1,3,5)])
    tmp = function(b) {
        moms = colSums(cbind(T_gps.hat[,1]*(y-b[1])-x%*%lm((y[T==1]-b[1])~x[T==1,]-1)$coeff*(T_gps.hat[,1]-1),
                             T_gps.hat[,2]*(y-(b[1]+b[2]))-x%*%lm((y[T==2]-(b[1]+b[2]))~x[T==2,]-1)$coeff*(T_gps.hat[,2]-1),
                             T_gps.hat[,3]*(y-(b[1]+2*b[2]))-x%*%lm((y[T==3]-(b[1]+2*b[2]))~x[T==3,]-1)$coeff*(T_gps.hat[,3]-1)))
        return(crossprod(moms,W)%*%moms)
    }
    EIF.m3 = optim(par=c(.5,.5),fn=tmp,method="Nelder-Mead")
    Table.model3.EIF[s,c(1,3)] = EIF.m3$par

    W = solve(V.EIF.m1[c(2,4,6),c(2,4,6)])
    tmp = function(b) {
        moms = colSums(cbind(T_gps.hat[,1]*((y<b[1])-.5)-x%*%lm(((y[T==1]<b[1])-.5)~x[T==1,]-1)$coeff*(T_gps.hat[,1]-1),
                             T_gps.hat[,2]*((y<(b[1]+b[2]))-.5)-x%*%lm(((y[T==2]<(b[1]+b[2]))-.5)~x[T==2,]-1)$coeff*(T_gps.hat[,2]-1),
                             T_gps.hat[,3]*((y<(b[1]+2*b[2]))-.5)-x%*%lm(((y[T==3]<(b[1]+2*b[2]))-.5)~x[T==3,]-1)$coeff*(T_gps.hat[,3]-1)))
        return(crossprod(moms,W)%*%moms)
    }
    EIF.m3 = optim(par=c(.5,.5),fn=tmp,method="Nelder-Mead")
    Table.model3.EIF[s,c(2,4)] = EIF.m3$par
    }

    ## MODEL 4
    ## IPW Estimator - Overidentification & Cross-equation Restriction
    if (do.IPW) {
    W = solve(V.IPW.m1)
    tmp = function(b) {
        moms = colSums(cbind(T_gps.hat[,1]*(y-b[1]),
                             T_gps.hat[,1]*((y<b[1])-.5),
                             T_gps.hat[,2]*(y-(b[1]+b[2])),
                             T_gps.hat[,2]*((y<(b[1]+b[2]))-.5),
                             T_gps.hat[,3]*(y-(b[1]+2*b[2])),
                             T_gps.hat[,3]*((y<(b[1]+2*b[2]))-.5)))
        return(crossprod(moms,W)%*%moms)
    }
    IPW.m4 = optim(par=c(.5,.5),fn=tmp,method="Nelder-Mead")
    Table.model4.IPW[s,] = IPW.m4$par
    }

    ## EIF Estimator - Overidentification & Cross-equation Restriction
    if (do.EIF) {
    W = solve(V.EIF.m1)
    tmp = function(b) {
        moms = colSums(cbind(T_gps.hat[,1]*(y-b[1])-x%*%lm((y[T==1]-b[1])~x[T==1,]-1)$coeff*(T_gps.hat[,1]-1),
                             T_gps.hat[,1]*((y<b[1])-.5)-x%*%lm(((y[T==1]<b[1])-.5)~x[T==1,]-1)$coeff*(T_gps.hat[,1]-1),
                             T_gps.hat[,2]*(y-(b[1]+b[2]))-x%*%lm((y[T==2]-(b[1]+b[2]))~x[T==2,]-1)$coeff*(T_gps.hat[,2]-1),
                             T_gps.hat[,2]*((y<(b[1]+b[2]))-.5)-x%*%lm(((y[T==2]<(b[1]+b[2]))-.5)~x[T==2,]-1)$coeff*(T_gps.hat[,2]-1),
                             T_gps.hat[,3]*(y-(b[1]+2*b[2]))-x%*%lm((y[T==3]-(b[1]+2*b[2]))~x[T==3,]-1)$coeff*(T_gps.hat[,3]-1),
                             T_gps.hat[,3]*((y<(b[1]+2*b[2]))-.5)-x%*%lm(((y[T==3]<(b[1]+2*b[2]))-.5)~x[T==3,]-1)$coeff*(T_gps.hat[,3]-1)))
        return(crossprod(moms,W)%*%moms)
    }
    EIF.m4 = optim(par=c(.5,.5),fn=tmp,method="Nelder-Mead")
    Table.model4.EIF[s,] = EIF.m4$par
    }
}
cat("\n")


################################################################################
## Figures
################################################################################
colores = c("blue", "red", "black", "darkgreen")
ltyset  = c(1,2,3,4)
lwdset  = c(1,1,1,1)
pchset  = c(21,22,23,24)

## FIGURE 1: IPW (T2-T1)
if (do.IPW) {
postscript(file="./Figure1-IPW-T2-T1.ps", horizontal=F)
to.plot = list(density(Table.model1.IPW[,5]-Table.model1.IPW[,3], adjust=1.50),
               density(Table.model2.IPW[,3]-Table.model2.IPW[,2]),
               density(Table.model3.IPW[,3]),
               density(Table.model4.IPW[,2]))
matplot(cbind(to.plot[[1]]$x,to.plot[[2]]$x,to.plot[[3]]$x,to.plot[[4]]$x),
        cbind(to.plot[[1]]$y,to.plot[[2]]$y,to.plot[[3]]$y,to.plot[[4]]$y),
        type="l", lty=ltyset, lwd=lwdset, col=colores, ylab="", xlab="", xlim=c(0.5,1.5))
abline(v=1, col="grey")
legend("topright",c("Model 1","Model 2","Model 3","Model 4"),lty=ltyset, lwd=lwdset, col=colores)
dev.off()
}

## FIGURE 2: EIF (T2-T1)
if (do.EIF) {
postscript(file="./Figure2-EIF-T2-T1.ps", horizontal=F)
to.plot = list(density(Table.model1.EIF[,5]-Table.model1.EIF[,3], adjust=1.50),
               density(Table.model2.EIF[,3]-Table.model2.EIF[,2]),
               density(Table.model3.EIF[,3]),
               density(Table.model4.EIF[,2]))
matplot(cbind(to.plot[[1]]$x,to.plot[[2]]$x,to.plot[[3]]$x,to.plot[[4]]$x),
        cbind(to.plot[[1]]$y,to.plot[[2]]$y,to.plot[[3]]$y,to.plot[[4]]$y),
        type="l", lty=ltyset, lwd=lwdset, col=colores, ylab="", xlab="", xlim=c(0.5,1.5))
abline(v=1, col="grey")
legend("topright",c("Model 1","Model 2","Model 3","Model 4"),lty=ltyset, lwd=lwdset, col=colores)
dev.off()
}

## FIGURE 3: BOTH (T2-T1)
postscript(file="./Figure3-IPW-T2-T1.ps", horizontal=T)
par(mfcol=c(1,2))
to.plot = list(density(Table.model1.IPW[,5]-Table.model1.IPW[,3], adjust=1.50),
               density(Table.model2.IPW[,3]-Table.model2.IPW[,2]),
               density(Table.model3.IPW[,3]),
               density(Table.model4.IPW[,2]))
matplot(cbind(to.plot[[1]]$x,to.plot[[2]]$x,to.plot[[3]]$x,to.plot[[4]]$x),
        cbind(to.plot[[1]]$y,to.plot[[2]]$y,to.plot[[3]]$y,to.plot[[4]]$y),
        type="l", lty=ltyset, lwd=lwdset, col=colores, ylab="", xlab="", xlim=c(0.5,1.5),
        main="IPW Estimator")
abline(v=1, col="grey")
legend("topright",c("Model 1","Model 2","Model 3","Model 4"),lty=ltyset, lwd=lwdset, col=colores)

to.plot = list(density(Table.model1.EIF[,5]-Table.model1.EIF[,3], adjust=1.50),
               density(Table.model2.EIF[,3]-Table.model2.EIF[,2]),
               density(Table.model3.EIF[,3]),
               density(Table.model4.EIF[,2]))
matplot(cbind(to.plot[[1]]$x,to.plot[[2]]$x,to.plot[[3]]$x,to.plot[[4]]$x),
        cbind(to.plot[[1]]$y,to.plot[[2]]$y,to.plot[[3]]$y,to.plot[[4]]$y),
        type="l", lty=ltyset, lwd=lwdset, col=colores, ylab="", xlab="", xlim=c(0.5,1.5),
        main="EIF Estimator")
abline(v=1, col="grey")
legend("topright",c("Model 1","Model 2","Model 3","Model 4"),lty=ltyset, lwd=lwdset, col=colores)
dev.off()

