################################################################################
## R CODE FOR CATTANEO (2007)
## DATE: 15-Jun-2007
################################################################################
rm(list=ls(all=TRUE))

## Server:
#install.packages("rgenoud",destdir="~/Rlibloc",lib="~/Rlibloc",dependencies=T)
library(SparseM,lib.loc="~/Rlibloc"); library(quantreg,lib.loc="~/Rlibloc"); library(rgenoud,lib.loc="~/Rlibloc")

library(foreign);
library(nnet); library(splines); library(quantreg); library(Hmisc)

################################################################################
## General Functions
################################################################################
## Just-identified Case
m.IPW = function(m,bt,pt,t) return(((T==T.levels[t])/pt)*m(Y,bt))

m.EIF = function(m,bt,pt,t,e) return(((T==T.levels[t])/pt)*m(Y,bt)-e(m,bt,t)*((T==T.levels[t])/pt-1))

V.EIF.JID = function(m,b,p,e) {
        out = matrix(NA,nrow=n,ncol=length(m)*d.T)
        for (t in 1:d.T) for (i in 1:length(m)) out[,i+(t-1)*d.bt] = m.EIF(m=m[[i]],b[i+(t-1)*d.bt],p[,t],t,e)
        #return(tcrossprod(solve(qr.R(qr(out))))*n)
        return(crossprod(out)/n)
}

## Over-identified Case
M.IPW = function(m,b,p) {
        out = matrix(NA,nrow=1,ncol=d.mt*d.T)
        for (t in 1:d.T) out[(1+(t-1)*d.mt):(t*d.mt)] = colSums(m.IPW(m,b[(1+(t-1)*d.bt):(t*d.bt)],p[,t],t))
        return(out)
}

M.EIF = function(m,b,p,e) {
        out = rep.int(NA,d.m*d.T)
        for (t in 1:d.T) out[(1+(t-1)*d.mt):(t*d.mt)] = colSums(m.EIF(m,b[(1+(t-1)*d.bt):(t*d.bt)],p[,t],t,e))
        return(out)
}

Vinv.EIF.OID = function(m,b,p,e) {
        out = matrix(NA,nrow=n,ncol=d.mt*d.T)
        for (t in 1:d.T) out[,(1+(t-1)*d.mt):(t*d.mt)] = m.EIF(m,b[(1+(t-1)*d.bt):(t*d.bt)],p[,t],t,e)
        #return(out)
        return(tcrossprod(solve(qr.R(qr(out))))*n)
        #return(crossprod(out)/n)
}

## Numerical Derivative
Gamma = function(m,d.m,b,d.bt,e,h){
        out = matrix(0,nrow=d.bt*d.T,ncol=d.m*d.T)
        for (t in 1:d.T) for (k in 1: d.bt){
                Dh = rep.int(0,d.bt); Dh[k] = h
                out[((t-1)*d.bt+k),(1+(t-1)*d.m):(t*d.m)] = 
                        colMeans((e(m,b[(1+(t-1)*d.bt):(t*d.bt)]+Dh,t)-e(m,b[(1+(t-1)*d.bt):(t*d.bt)]-Dh,t))/(2*h))
        }
        return(out)

}

## Put starts to coeff
putstars = function(x) {
        t = abs(x[1]/x[2]); stars="";
        if (qnorm(.925)<=t & t<qnorm(.975)) stars=""
        if (qnorm(.975)<=t & t<qnorm(.995)) stars="*"
        if (qnorm(.995)<=t                ) stars="*"
        return(paste(round(x[1]),stars,sep=""))
}

################################################################################
## Load Data - 5-cigarette bins
################################################################################
data = read.dta("data-5cigar.dta")
Y = as.vector(data[,1]); T = as.vector(data[,2])
Xd = as.matrix(data[,3:46]); Xc = as.matrix(data[,47:53])
rm(data)
n = length(Y); d.Xd = ncol(Xd); d.Xc = ncol(Xc); T.levels = sort(unique(T)); d.T = length(T.levels)

################################################################################
## Nonparametric Estimation - 5-cigarette bins
################################################################################
## Series Basis
knots = c(3,1,3,1,3,1,1)
for (j in 1:d.Xc) {
        if (j==1) R=cbind(Xd, bs(knots=quantile(Xc[,j], probs=seq(1/(knots[j]+1),1-1/(knots[j]+1),length.out=knots[j]),type=1),x=Xc[,j]))
        else R=cbind(R,bs(knots=quantile(Xc[,j], probs=seq(1/(knots[j]+1),1-1/(knots[j]+1),length.out=knots[j]),type=1),x=Xc[,j]))
}
d.R=ncol(R)

## Conditional Expectation (e)
Rt = list(NA)
for (t in 1:d.T) Rt[[t]] = tcrossprod(tcrossprod(solve(qr.R(qr(R[T==T.levels[t],])))),R[T==T.levels[t],])
e.hat = function(m,bt,t) return(R%*%(Rt[[t]]%*%m(Y[T==T.levels[t]],bt)))

## Generalized Propensity Score (p)
gps.fit = multinom(T ~ R, maxit=5000)
#apply(gps.fit$fitted.values,2,summary)
#boxplot(data.frame(gps.fit$fitted.values), notch=T)
#write.csv(gps.fit$fitted.values,file="gps.hat-5cigar.csv",row.names=F)
#gps.hat = read.csv(file="gps.hat.csv")
gps.hat = gps.fit$fitted.values

################################################################################
## Estimation - 5-cigarette bins
################################################################################
d.mt = 7; d.bt = 7; taus = c(.9,.75,.5,1,.25,.1,2)
mt.JID.list = list(function(y,b) ((y<b)-.9),function(y,b) ((y<b)-.75),function(y,b) ((y<b)-.5),function(y,b) (y-b),function(y,b) ((y<b)-.25),function(y,b) ((y<b)-.1),function(y,b) (y^2-b))
mt.JID.matx = function(y,bt) return(cbind((y<bt[1])-.9,(y<bt[2])-.75,(y<bt[3])-.5,y-bt[4],(y<bt[5])-.25,(y<bt[6])-.1),y^2-bt[7])
Table1.rownames=c("0","","1-5","","6-10","","11-15","","16-20","","21+","")
Table1raw.colnames=rep(c("Q.9","Q.75","Q.5","Mean","Q.25","Q.1","2ndMom"),3)
Table1.colnames=rep(c("Q.9","Q.75","Q.5","Mean","Q.25","Q.1","SD"),3)

Table1raw = matrix(NA,nrow=d.T*2,ncol=3*d.bt,
                dimnames=list(Table1.rownames,Table1raw.colnames))
index.coef = seq(1,by=2,length.out=d.T)

## Dummy Regression
for (t in 1:d.T) for (i in 1:length(taus)){
        if (taus[i]<1)  fit = summary(rq(Y~as.numeric(T==T.levels[t])-1,tau=taus[i]),se="ker")
        else fit = summary(lm(Y^taus[i]~as.numeric(T==T.levels[t])-1))
        Table1raw[(2*t-1),i] = round(fit$coef[,1]); Table1raw[(2*t),i] = round(fit$coef[,2])
        cat(paste("\nTreatment",t,"of",d.T,"| Parameter",i,"of",length(mt.JID.list),": Done!"))
}
Table1 = Table1raw; colnames(Table1)=Table1.colnames
Table1[,7] = round(sqrt(Table1raw[,7]-(Table1raw[,4])^2)); Table1[index.coef+1,7] = NA

###########################
## IPW Estimator
###########################
interval.bt = cbind(c(rep(1500,d.bt-1),9000000),c(rep(5500,d.bt-1),13000000))
for (t in 1:d.T) for (i in 1:length(mt.JID.list)){
        IPW.JID = optimize(f=function(b) return(abs(sum(m.IPW(m=mt.JID.list[[i]],b,p=gps.hat[,t],t=t)))),
                           interval=interval.bt[i,])
        Table1raw[(2*t-1),d.bt+i] = round(IPW.JID$minimum);
        cat(paste("\nTreatment",t,"of",d.T,"| Parameter",i,"of",length(mt.JID.list),": Done!"))
}

b.IPW=c(t(Table1raw[index.coef,(d.bt+1):(2*d.bt)]))
# Using numerical derivative
#G = Gamma(m=mt.JID.matx,d.m=6,b=b,d.bt=6,e=e.hat,h=10)
# Using exact formula: density at q or -1
G.IPW = matrix(0,nrow=d.bt*d.T,ncol=d.mt*d.T); h = 1.06*sd(Y)*n^(-1/5)
for (t in 1:d.T) for (i in 1:d.bt) {
        if (taus[i]<1) G.IPW[(t-1)*d.bt+i,(t-1)*d.bt+i]=mean(((T==T.levels[t])/gps.hat[,t])*dnorm((Y-b[(t-1)*d.bt+i])/h)/h)
        else G.IPW[(t-1)*d.bt+i,(t-1)*d.bt+i]=-1
}
V.b.IPW = V.EIF.JID(m=mt.JID.list,b=b,p=gps.hat,e=e.hat)
SPEB.IPW.JID = solve(G.IPW)%*%V.b.IPW%*%solve(G.IPW)/n
Table1raw[index.coef+1,(d.bt+1):(2*d.bt)] = t(matrix(round(sqrt(diag(SPEB.IPW.JID))),nrow=d.bt,ncol=d.T))
write.csv(Table1raw, file="Table1raw-5cigar.csv")

Table1[index.coef,(d.bt+1):(2*d.bt)]=Table1raw[index.coef,(d.bt+1):(2*d.bt)]
Table1[,14] = round(sqrt(Table1raw[,14]-(Table1raw[,11])^2))
Df = diag(1,d.T*d.bt)
diag(Df[seq(4,by=d.bt,length.out=d.T),seq(7,by=7,length.out=d.T)]) = -c(Table1raw[index.coef,11])/c(Table1[index.coef,14])
diag(Df[seq(7,by=d.bt,length.out=d.T),seq(7,by=7,length.out=d.T)]) = 1/(2*c(Table1[index.coef,14]))

SPEB.IPW.JID.f = crossprod(Df,SPEB.IPW.JID)%*%Df
Table1[index.coef+1,(d.bt+1):(2*d.bt)] = t(matrix(round(sqrt(diag(SPEB.IPW.JID.f))),nrow=d.bt,ncol=d.T))
write.csv(Table1, file="Table1-5cigar.csv")

###########################
## EIF Estimator
###########################
for (t in 1:d.T) for (i in 1:length(mt.JID.list)){
        EIF.JID = optimize(f=function(b) return(abs(sum(m.EIF(m=mt.JID.list[[i]],b,p=gps.hat[,t],t=t,e=e.hat)))),
                           interval=interval.bt[i,])
        Table1raw[(2*t-1),2*d.bt+i] = round(EIF.JID$minimum);
        cat(paste("\nTreatment",t,"of",d.T,"| Parameter",i,"of",length(mt.JID.list),": Done!"))
}

b.EIF=c(t(Table1raw[index.coef,(2*d.bt+1):(3*d.bt)]))
G.EIF = matrix(0,nrow=d.bt*d.T,ncol=d.mt*d.T); h = 1.06*sd(Y)*n^(-1/5)
for (t in 1:d.T) for (i in 1:d.bt) {
        if (taus[i]<1) G.EIF[(t-1)*d.bt+i,(t-1)*d.bt+i]=mean(((T==T.levels[t])/gps.hat[,t])*dnorm((Y-b[(t-1)*d.bt+i])/h)/h)
        else G.EIF[(t-1)*d.bt+i,(t-1)*d.bt+i]=-1
}
V.b.EIF = V.EIF.JID(m=mt.JID.list,b=b,p=gps.hat,e=e.hat)
SPEB.EIF.JID = solve(G.EIF)%*%V.b.EIF%*%solve(G.EIF)/n
Table1raw[index.coef+1,(2*d.bt+1):(3*d.bt)] = t(matrix(round(sqrt(diag(SPEB.EIF.JID))),nrow=d.bt,ncol=d.T))
write.csv(Table1raw, file="Table1raw-5cigar.csv")

Table1[index.coef,(2*d.bt+1):(3*d.bt)]=Table1raw[index.coef,(2*d.bt+1):(3*d.bt)]
Table1[,21] = round(sqrt(Table1raw[,21]-(Table1raw[,18])^2))
Df = diag(1,d.T*d.bt)
diag(Df[seq(4,by=d.bt,length.out=d.T),seq(7,by=7,length.out=d.T)]) = -c(Table1raw[index.coef,18])/c(Table1[index.coef,21])
diag(Df[seq(7,by=d.bt,length.out=d.T),seq(7,by=7,length.out=d.T)]) = 1/(2*c(Table1[index.coef,21]))

SPEB.EIF.JID.f = crossprod(Df,SPEB.EIF.JID)%*%Df
Table1[index.coef+1,(2*d.bt+1):(3*d.bt)] = t(matrix(round(sqrt(diag(SPEB.EIF.JID.f))),nrow=d.bt,ncol=d.T))
write.csv(Table1, file="Table1-5cigar.csv")

## Latex Table
Table1tex = Table1
Table1tex[index.coef+1,] = paste("(",Table1[index.coef+1,],")",sep="")
Table1tex[index.coef+1,7] = "n.a."
colnames(Table1tex)=Table1.colnames
rownames(Table1tex)=Table1.rownames
latex(Table1tex, title="", file="table1-5cigar.tex",
      n.group=c(1:7,8:14,15:21),cgroup=c("DRE","IPWE","EIFE"),
      size="scriptsize",
      caption = "Effect of Maternal Smoking Intensity on Birth Weight",
      insert.bottom="Notes: (i) DRE = Dummy Regression Estimator, IPWE = Inverse Probability Weighting Estimator, EIFE = Eficient Influence Function Estimator; (ii) Q.9, Q.75, Q.5, Q.25 and Q.1 are the 90%, 75%, 50%, 25% and 10% quantiles, respectively, and SD is the standard deviation; (iii) standard errors in parentheses."
)

################################################################################
## Figures - 5-cigarette bins
################################################################################
#Table1 = as.matrix(read.csv(file="Table1.csv")[1:(2*d.T),2:22])

## IPW: Figure1
postscript(file="Figure1-IPW-JID-5cigar.ps", horizontal=F)
index.colstoplot = (d.bt+1):(2*d.bt-1)
plotdata.IPW = rep(1:d.T,6)
plotdata.IPW = cbind(plotdata.IPW,c(Table1[index.coef,index.colstoplot]),c(Table1[index.coef+1,index.colstoplot]))
plot(plotdata.IPW[,1:2], axes=F,
       #main="Effect of Maternal Smoking Intensity on Birth Weight \n (Just-identified IPW Estimates)",
       ylab="Birth weight", xlab="Cigarettes-per-day Smoked During Pregnancy")
box(); axis(1, at=1:d.T, lab=Table1.rownames[!(Table1.rownames=="")]); axis(2)
lty.setting=c(2,4,3,1,4,2)
for (i in 1:6) lines(plotdata.IPW[seq(1+(i-1)*d.T,length.out=d.T),1:2], lty=lty.setting[i])
for (i in 1:6) for (j in 1:d.T) lines(c(j,j),c(plotdata.IPW[j+(i-1)*d.T,2]-1.96*plotdata.IPW[j+(i-1)*d.T,3],
                                               plotdata.IPW[j+(i-1)*d.T,2]+1.96*plotdata.IPW[j+(i-1)*d.T,3]))
legend("topright",c("Q.9 & Q.10","Q.75 & Q.25","Median","Mean"),lty=lty.setting[1:4])
dev.off()

################################################################################
## Hypothesis Tests Tables - 5-cigarette bins
################################################################################
# Table 2: Pairwise differences of MMTE
Df = matrix(0,nrow=d.bt*d.T,ncol=choose(d.T,2)); colpos = 1
Table2.names = rep("",choose(d.T,2))
for (i in 1:(d.T-1)) for (j in 1:(d.T-i)){
        Df[c((i-1)*d.bt+4,(i-1+j)*d.bt+4),colpos]=c(-1,1)
        Table2.names[colpos] = paste("T",j+i-1,"-T",i-1,sep="")
        colpos = colpos+1
}

PWDiff.b = crossprod(Df,c(t(Table1[index.coef,(d.bt+1):(2*d.bt)])))
PWDiff.V = crossprod(Df,SPEB.IPW.JID.f)%*%Df

Table2 = matrix("",nrow=choose(d.T,2),ncol=choose(d.T,2), dimnames=list(Table2.names,Table2.names))

diag(Table2) = apply(cbind(PWDiff.b,sqrt(diag(PWDiff.V))),1,putstars)

for (i in 1:(choose(d.T,2)-1)) {
        Df = matrix(0,nrow=choose(d.T,2),ncol=choose(d.T,2)-i); Df[i,] = -1
        if (length((i+1):choose(d.T,2))>1) diag(Df[(i+1):choose(d.T,2),]) = 1
        else Df[(i+1):choose(d.T,2),] = 1
        Table2[i,(i+1):(choose(d.T,2))] = apply(cbind(crossprod(Df,PWDiff.b),sqrt(diag(crossprod(Df,PWDiff.V)%*%Df))),
                                            1,putstars)
}
write.csv(Table2, file="Table2-5cigar.csv")

## Latex Table
Table2tex = Table2
Table2tex[1,2:5]=""; Table2tex[2,3:5]=""; Table2tex[3,4:5]=""; Table2tex[4,5]=""
Table2tex[6,7:9]=""; Table2tex[7,8:9]=""; Table2tex[8,9]="";
Table2tex[10,11:12]=""; Table2tex[11,12]=""
Table2tex[13,14]=""

latex(Table2tex, title="", file="table2-5cigar.tex",
      size="scriptsize",
      caption = "Hypothesis Tests for Pairwise Differences and Difference-in-Differences Effects",
      insert.bottom="Notes: (i) treatments T0, T1, T2, T3, T4 and T5 are $0$, $1$-$5$, $6$-$10$, $11$-$15$, $16$-$20$ and $21$+ cigarettes-per-day smoked, respectively; (ii) pairwise differences are reported on the diagonal, and difference-in-differences are reported outside the diagonal; (iii) in all cases the null hypothesis is zero differential effect; (iv) $^{*}$ significant at $5%$."
)

################################################
### Other Hypothesis Tests
################################################
table3.colnames=c("Number of Restrictions", "Wald Test Statistic", "p-value")
table3.rownames=c("Equal treatment effects (mean, quantiles, spread) for (11-15,16-20,21+)",
                  "Equal treatment effects (mean, quantiles, spread) for (6-10,11-15,16-20,21+)",
                  "Equal treatment effects (mean, quantiles, spread) for (1-5,6-10,11-15,16-20,21+)",
                  "Equal mean and median for each treatment",
                  "Equal mean-median difference (MMD) across treatments",
                  "Equal standard deviation across treatments",
                  "Equal interquartile range (IQR) across treatments",
                  "Equal Q.9-Q.1 range (Q.9-Q.1) across treatments",
                  "Equal MMD, IQR and Q.9-Q.1 across treatments")
table3tex=matrix(NA,nrow=length(table3.rownames),ncol=3, dimnames=list(table3.rownames,table3.colnames))
# Equal Treatment effects on 11-15 to 21+
d.R = 2*7; fila=1
R = matrix(0,nrow=d.R,ncol=d.bt*d.T)
diag(R[,(d.bt*3+1):(d.bt*(d.T-1))])=1
diag(R[,(d.bt*4+1):(d.bt*d.T)])=-1
RB=R%*%b.IPW; Vinv = solve(tcrossprod(R%*%SPEB.IPW.JID,R))
W = crossprod(RB,Vinv%*%RB)
table3tex[fila,1] = d.R
table3tex[fila,2] = round(W, digits=2)
table3tex[fila,3] = round(1-pchisq(W,df=d.R), digits=4)

# Equal Treatment effects on 6-10 to 21+
d.R = 3*7; fila=2
R = matrix(0,nrow=d.R,ncol=d.bt*d.T)
diag(R[,(d.bt*2+1):(d.bt*(d.T-1))])=1
diag(R[,(d.bt*3+1):(d.bt*d.T)])=-1
RB=R%*%b.IPW; Vinv = solve(tcrossprod(R%*%SPEB.IPW.JID,R))
W = crossprod(RB,Vinv%*%RB)
table3tex[fila,1] = d.R
table3tex[fila,2] = round(W, digits=2)
table3tex[fila,3] = round(1-pchisq(W,df=d.R), digits=4)

# Equal Treatment effects on 1-5 to 21+
d.R = 4*7; fila=3
R = matrix(0,nrow=d.R,ncol=d.bt*d.T)
diag(R[,(d.bt*1+1):(d.bt*(d.T-1))])=1
diag(R[,(d.bt*2+1):(d.bt*d.T)])=-1
RB=R%*%b.IPW; Vinv = solve(tcrossprod(R%*%SPEB.IPW.JID,R))
W = crossprod(RB,Vinv%*%RB)
table3tex[fila,1] = d.R
table3tex[fila,2] = round(W, digits=2)
table3tex[fila,3] = round(1-pchisq(W,df=d.R), digits=4)

# Equal mean and median within treatment
d.R = 6; fila=4
R = matrix(0,nrow=d.R,ncol=d.bt*d.T)
for (i in 1:d.R) {
    R[i,3+(i-1)*d.bt]=1
    R[i,4+(i-1)*d.bt]=-1
}
RB=R%*%b.IPW; Vinv = solve(tcrossprod(R%*%SPEB.IPW.JID,R))
W = crossprod(RB,Vinv%*%RB)
table3tex[fila,1] = d.R
table3tex[fila,2] = round(W, digits=2)
table3tex[fila,3] = round(1-pchisq(W,df=d.R), digits=4)

# Equal diff mean-median across treatment
d.R = d.T-1; fila=5
R = matrix(0,nrow=d.R,ncol=d.bt*d.T)
for (i in 1:d.R) {
    R[i,3+(i-1)*d.bt]=1
    R[i,4+(i-1)*d.bt]=-1
    R[i,3+i*d.bt]=-1
    R[i,4+i*d.bt]=1
}
RB=R%*%b.IPW; Vinv = solve(tcrossprod(R%*%SPEB.IPW.JID,R))
W = crossprod(RB,Vinv%*%RB)
table3tex[fila,1] = d.R
table3tex[fila,2] = round(W, digits=2)
table3tex[fila,3] = round(1-pchisq(W,df=d.R), digits=4)

# Equal SD across treatment
d.R = d.T-1; fila=6
R = matrix(0,nrow=d.R,ncol=d.bt*d.T)
for (i in 1:d.R) {
    R[i,7+(i-1)*d.bt]=1
    R[i,7+i*d.bt]=-1
}
RB=R%*%c(t(Table1[index.coef,(d.bt+1):(2*d.bt)])); Vinv = solve(tcrossprod(R%*%SPEB.IPW.JID.f,R))
W = crossprod(RB,Vinv%*%RB)
table3tex[fila,1] = d.R
table3tex[fila,2] = round(W, digits=2)
table3tex[fila,3] = round(1-pchisq(W,df=d.R), digits=4)

# Equal IQR across treatment
d.R = d.T-1; fila=7
R = matrix(0,nrow=d.R,ncol=d.bt*d.T)
for (i in 1:d.R) {
    R[i,2+(i-1)*d.bt]=1
    R[i,5+(i-1)*d.bt]=-1
    R[i,2+i*d.bt]=-1
    R[i,5+i*d.bt]=1
}
RB=R%*%b.IPW; Vinv = solve(tcrossprod(R%*%SPEB.IPW.JID,R))
W = crossprod(RB,Vinv%*%RB)
table3tex[fila,1] = d.R
table3tex[fila,2] = round(W, digits=2)
table3tex[fila,3] = round(1-pchisq(W,df=d.R), digits=4)

# Equal Q.9-Q.1 across treatment
d.R = d.T-1; fila=8
R = matrix(0,nrow=d.R,ncol=d.bt*d.T)
for (i in 1:d.R) {
    R[i,1+(i-1)*d.bt]=1
    R[i,6+(i-1)*d.bt]=-1
    R[i,1+i*d.bt]=-1
    R[i,6+i*d.bt]=1
}
RB=R%*%b.IPW; Vinv = solve(tcrossprod(R%*%SPEB.IPW.JID,R))
W = crossprod(RB,Vinv%*%RB)
table3tex[fila,1] = d.R
table3tex[fila,2] = round(W, digits=2)
table3tex[fila,3] = round(1-pchisq(W,df=d.R), digits=4)

# Equal mean-median diff, IQR, Q.9-Q.1
d.R = 3*(d.T-1); fila=9
R = matrix(0,nrow=d.R,ncol=d.bt*d.T)
for (i in 1:(1*(d.T-1))) {
    R[i,3+(i-1)*d.bt]=1
    R[i,4+(i-1)*d.bt]=-1
    R[i,3+i*d.bt]=-1
    R[i,4+i*d.bt]=1
}
for (i in (1*(d.T-1)+1):(2*(d.T-1))) {
    R[i,2+(i-(1*(d.T-1))-1)*d.bt]=1
    R[i,5+(i-(1*(d.T-1))-1)*d.bt]=-1
    R[i,2+(i-(1*(d.T-1)))*d.bt]=-1
    R[i,5+(i-(1*(d.T-1)))*d.bt]=1
}
for (i in (2*(d.T-1)+1):(3*(d.T-1))) {
    R[i,1+(i-(2*(d.T-1))-1)*d.bt]=1
    R[i,6+(i-(2*(d.T-1))-1)*d.bt]=-1
    R[i,1+(i-(2*(d.T-1)))*d.bt]=-1
    R[i,6+(i-(2*(d.T-1)))*d.bt]=1
}
RB=R%*%b.IPW; Vinv = solve(tcrossprod(R%*%SPEB.IPW.JID,R))
W = crossprod(RB,Vinv%*%RB)
table3tex[fila,1] = d.R
table3tex[fila,2] = round(W, digits=2)
table3tex[fila,3] = round(1-pchisq(W,df=d.R), digits=4)

latex(table3tex, title="Joint Null Hypotheses", file="table3-5cigar.tex",
      size="scriptsize",
      caption = "Joint Hypotheses Tests (IPWE)",
      insert.bottom="Note: all tests have been computed using the IPWE and its corresponding limiting distribution."
)

################################################################################
## Estimation - OVER-IDENTIFIED CASE - Symmetric Distribution
################################################################################
# d.mt = 7; d.bt = 6
# mt.OID = function(y,bt) return(cbind(((y<bt[1])-.9),((y<bt[2])-.75),((y<bt[3])-.5),(y-bt[3]),((y<bt[4])-.25),((y<bt[5])-.1),(y^2-bt[6])))

# Table3raw = matrix(NA,nrow=d.T*2,ncol=2*d.bt,
#                    dimnames=list(Table1.rownames,rep(c("Q.9","Q.75","Q.5/Mean","Q.25","Q.1","2ndMom"),2)))

# index.coef = seq(1,by=2,length.out=d.T);

## IPW: OID Estimator
# b=c(t(Table1raw[index.coef,c(8:10,12:14)]))
# W = Vinv.EIF.OID(m=mt.OID,b=b,p=gps.hat,e=e.hat)
# IPW.OID = optim(par=b,
#                 fn=function(b){M=M.IPW(m=mt.OID,b=b,p=gps.hat); return(M%*%tcrossprod(W,M))})

# Table3raw[index.coef,1:d.bt] = t(matrix(round(IPW.OID$par),nrow=d.bt,ncol=d.T))

# Vinv.b.IPW = Vinv.EIF.OID(m=mt.OID,b=IPW.OID$par,p=gps.hat,e=e.hat)
# index.G = 1:(d.mt*d.T)
# for (i in 1:d.T) index.G = index.G[index.G!=(4+(i-1)*d.mt)]
# G.IPW.OID = G.IPW[,index.G]
# for (i in 1:d.T) G.IPW.OID[4+(i-1)*d.mt,3+(i-1)*d.bt] = -1
# SPEB.IPW.JID = solve(crossprod(G.IPW.OID,V.b.EIF)%*%G.IPW.OID)/n


################################################################################
## Load Data - 2-cigarette bin on [0,10]
################################################################################
data = read.dta("data-2cigar.dta")
#data[data[,2]==10,2]=9
Y = as.vector(data[,1]); T = as.vector(data[,2])
Xd = as.matrix(data[,3:46]); Xc = as.matrix(data[,47:53])
rm(data)
n = length(Y); d.Xd = ncol(Xd); d.Xc = ncol(Xc); T.levels = sort(unique(T)); d.T = length(T.levels)

################################################################################
## Nonparametric Estimation - 2-cigarette bin on [0,10]
################################################################################
## Series Basis
knots = c(3,1,3,1,3,1,1)
for (j in 1:d.Xc) {
        if (j==1) R=cbind(Xd, bs(knots=quantile(Xc[,j], probs=seq(1/(knots[j]+1),1-1/(knots[j]+1),length.out=knots[j]),type=1),x=Xc[,j]))
        else R=cbind(R,bs(knots=quantile(Xc[,j], probs=seq(1/(knots[j]+1),1-1/(knots[j]+1),length.out=knots[j]),type=1),x=Xc[,j]))
}
d.R=ncol(R)

## Conditional Expectation (e)
Rt = list(NA)
for (t in 1:d.T) Rt[[t]] = tcrossprod(tcrossprod(solve(qr.R(qr(R[T==T.levels[t],])))),R[T==T.levels[t],])
e.hat = function(m,bt,t) return(R%*%(Rt[[t]]%*%m(Y[T==T.levels[t]],bt)))

## Generalized Propensity Score (p)
gps.fit = multinom(T ~ R, maxit=5000)
#apply(gps.fit$fitted.values,2,summary)
#boxplot(data.frame(gps.fit$fitted.values), notch=T)
#write.csv(gps.fit$fitted.values,file="gps.hat-1cigar.csv",row.names=F)
#gps.hat = read.csv(file="gps.hat.csv")
gps.hat = gps.fit$fitted.values

################################################################################
## Estimation - 2-cigarette bin on [0,10]
################################################################################
d.mt = 7; d.bt = 7; taus = c(.9,.75,.5,1,.25,.1,2)
mt.JID.list = list(function(y,b) ((y<b)-.9),function(y,b) ((y<b)-.75),function(y,b) ((y<b)-.5),function(y,b) (y-b),function(y,b) ((y<b)-.25),function(y,b) ((y<b)-.1),function(y,b) (y^2-b))
mt.JID.matx = function(y,bt) return(cbind((y<bt[1])-.9,(y<bt[2])-.75,(y<bt[3])-.5,y-bt[4],(y<bt[5])-.25,(y<bt[6])-.1),y^2-bt[7])
Table1.rownames=c("0","","1-2","","3-4","","5-6","","7-8","","9-10","")
Table1raw.colnames=rep(c("Q.9","Q.75","Q.5","Mean","Q.25","Q.1","2ndMom"),3)
Table1.colnames=rep(c("Q.9","Q.75","Q.5","Mean","Q.25","Q.1","SD"),3)

Table1raw = matrix(NA,nrow=d.T*2,ncol=3*d.bt,
                dimnames=list(Table1.rownames,Table1raw.colnames))
index.coef = seq(1,by=2,length.out=d.T)

## Dummy Regression
for (t in 1:d.T) for (i in 1:length(taus)){
        if (taus[i]<1)  fit = summary(rq(Y~as.numeric(T==T.levels[t])-1,tau=taus[i]),se="ker")
        else fit = summary(lm(Y^taus[i]~as.numeric(T==T.levels[t])-1))
        Table1raw[(2*t-1),i] = round(fit$coef[,1]); Table1raw[(2*t),i] = round(fit$coef[,2])
        cat(paste("\nTreatment",t,"of",d.T,"| Parameter",i,"of",length(mt.JID.list),": Done!"))
}
Table1 = Table1raw; colnames(Table1)=Table1.colnames
Table1[,7] = round(sqrt(Table1raw[,7]-(Table1raw[,4])^2)); Table1[index.coef+1,7] = NA

##########################
## IPW Estimator
##########################
interval.bt = cbind(c(rep(1500,d.bt-1),9000000),c(rep(5500,d.bt-1),13000000))
for (t in 1:d.T) for (i in 1:length(mt.JID.list)){
        IPW.JID = optimize(f=function(b) return(abs(sum(m.IPW(m=mt.JID.list[[i]],b,p=gps.hat[,t],t=t)))),
                           interval=interval.bt[i,])
        Table1raw[(2*t-1),d.bt+i] = round(IPW.JID$minimum);
        cat(paste("\nTreatment",t,"of",d.T,"| Parameter",i,"of",length(mt.JID.list),": Done!"))
}

b=c(t(Table1raw[index.coef,(d.bt+1):(2*d.bt)]))
# Using numerical derivative
#G = Gamma(m=mt.JID.matx,d.m=6,b=b,d.bt=6,e=e.hat,h=10)
# Using exact formula: density at q or -1
G.IPW = matrix(0,nrow=d.bt*d.T,ncol=d.mt*d.T); h = 1.06*sd(Y)*n^(-1/5)
for (t in 1:d.T) for (i in 1:d.bt) {
        if (taus[i]<1) G.IPW[(t-1)*d.bt+i,(t-1)*d.bt+i]=mean(((T==T.levels[t])/gps.hat[,t])*dnorm((Y-b[(t-1)*d.bt+i])/h)/h)
        else G.IPW[(t-1)*d.bt+i,(t-1)*d.bt+i]=-1
}
V.b.IPW = V.EIF.JID(m=mt.JID.list,b=b,p=gps.hat,e=e.hat)
SPEB.IPW.JID = solve(G.IPW)%*%V.b.IPW%*%solve(G.IPW)/n
Table1raw[index.coef+1,(d.bt+1):(2*d.bt)] = t(matrix(round(sqrt(diag(SPEB.IPW.JID))),nrow=d.bt,ncol=d.T))
write.csv(Table1raw, file="Table1raw-2cigar.csv")

Table1[index.coef,(d.bt+1):(2*d.bt)]=Table1raw[index.coef,(d.bt+1):(2*d.bt)]
Table1[,14] = round(sqrt(Table1raw[,14]-(Table1raw[,11])^2))
Df = diag(1,d.T*d.bt)
diag(Df[seq(4,by=d.bt,length.out=d.T),seq(7,by=7,length.out=d.T)]) = -c(Table1raw[index.coef,11])/c(Table1[index.coef,14])
diag(Df[seq(7,by=d.bt,length.out=d.T),seq(7,by=7,length.out=d.T)]) = 1/(2*c(Table1[index.coef,14]))

SPEB.IPW.JID.f = crossprod(Df,SPEB.IPW.JID)%*%Df
Table1[index.coef+1,(d.bt+1):(2*d.bt)] = t(matrix(round(sqrt(diag(SPEB.IPW.JID.f))),nrow=d.bt,ncol=d.T))
write.csv(Table1, file="Table1-2cigar.csv")

########################
## EIF Estimator
########################
for (t in 1:d.T) for (i in 1:length(mt.JID.list)){
        EIF.JID = optimize(f=function(b) return(abs(sum(m.EIF(m=mt.JID.list[[i]],b,p=gps.hat[,t],t=t,e=e.hat)))),
                           interval=interval.bt[i,])
        Table1raw[(2*t-1),2*d.bt+i] = round(EIF.JID$minimum);
        cat(paste("\nTreatment",t,"of",d.T,"| Parameter",i,"of",length(mt.JID.list),": Done!"))
}

b=c(t(Table1raw[index.coef,(2*d.bt+1):(3*d.bt)]))
G.EIF = matrix(0,nrow=d.bt*d.T,ncol=d.mt*d.T); h = 1.06*sd(Y)*n^(-1/5)
for (t in 1:d.T) for (i in 1:d.bt) {
        if (taus[i]<1) G.EIF[(t-1)*d.bt+i,(t-1)*d.bt+i]=mean(((T==T.levels[t])/gps.hat[,t])*dnorm((Y-b[(t-1)*d.bt+i])/h)/h)
        else G.EIF[(t-1)*d.bt+i,(t-1)*d.bt+i]=-1
}
V.b.EIF = V.EIF.JID(m=mt.JID.list,b=b,p=gps.hat,e=e.hat)
SPEB.EIF.JID = solve(G.EIF)%*%V.b.EIF%*%solve(G.EIF)/n
Table1raw[index.coef+1,(2*d.bt+1):(3*d.bt)] = t(matrix(round(sqrt(diag(SPEB.EIF.JID))),nrow=d.bt,ncol=d.T))
write.csv(Table1raw, file="Table1raw-2cigar.csv")

Table1[index.coef,(2*d.bt+1):(3*d.bt)]=Table1raw[index.coef,(2*d.bt+1):(3*d.bt)]
Table1[,21] = round(sqrt(Table1raw[,21]-(Table1raw[,18])^2))
Df = diag(1,d.T*d.bt)
diag(Df[seq(4,by=d.bt,length.out=d.T),seq(7,by=7,length.out=d.T)]) = -c(Table1raw[index.coef,18])/c(Table1[index.coef,21])
diag(Df[seq(7,by=d.bt,length.out=d.T),seq(7,by=7,length.out=d.T)]) = 1/(2*c(Table1[index.coef,21]))

SPEB.EIF.JID.f = crossprod(Df,SPEB.EIF.JID)%*%Df
Table1[index.coef+1,(2*d.bt+1):(3*d.bt)] = t(matrix(round(sqrt(diag(SPEB.EIF.JID.f))),nrow=d.bt,ncol=d.T))
write.csv(Table1, file="Table1-2cigar.csv")

## Latex Table
Table1tex = Table1
Table1tex[index.coef+1,] = paste("(",Table1[index.coef+1,],")",sep="")
Table1tex[index.coef+1,7] = "n.a."
colnames(Table1tex)=Table1.colnames
rownames(Table1tex)=Table1.rownames
latex(Table1tex, title="", file="table1-1cigar.tex",
      n.group=c(1:7,8:14,15:21),cgroup=c("DRE","IPWE","EIFE"),
      size="scriptsize",
      caption = "Effect of Maternal Smoking Intensity on Birth Weight",
      insert.bottom="Notes: (i) DRE = Dummy Regression Estimator, IPWE = Inverse Probability Weighting Estimator, EIFE = Eficient Influence Function Estimator; (ii) standard errors in parenthesis."
)


################################################################################
## Figures - 2-cigarette bin on [0,10]
################################################################################
#Table1 = as.matrix(read.csv(file="Table1.csv")[1:(2*d.T),2:22])

## IPW: Figure1
postscript(file="Figure1-IPW-JID-2cigar.ps", horizontal=F)
index.colstoplot = (d.bt+1):(2*d.bt-1)
plotdata.IPW = rep(1:d.T,6)
plotdata.IPW = cbind(plotdata.IPW,c(Table1[index.coef,index.colstoplot]),c(Table1[index.coef+1,index.colstoplot]))
plot(plotdata.IPW[,1:2], axes=F,
       #main="Effect of Maternal Smoking Intensity on Birth Weight \n (Just-identified IPW Estimates)",
       ylab="Birth weight", xlab="Cigarettes-per-day Smoked During Pregnancy")
box(); axis(1, at=1:d.T, lab=Table1.rownames[!(Table1.rownames=="")]); axis(2)
lty.setting=c(2,4,3,1,4,2)
for (i in 1:6) lines(plotdata.IPW[seq(1+(i-1)*d.T,length.out=d.T),1:2], lty=lty.setting[i])
for (i in 1:6) for (j in 1:d.T) lines(c(j,j),c(plotdata.IPW[j+(i-1)*d.T,2]-1.96*plotdata.IPW[j+(i-1)*d.T,3],
                                               plotdata.IPW[j+(i-1)*d.T,2]+1.96*plotdata.IPW[j+(i-1)*d.T,3]))
legend("topright",c("Q.9 & Q.10","Q.75 & Q.25","Median","Mean"),lty=lty.setting[1:4])
dev.off()











################################################################################
## Load Data - 1-cigarette bin on [0,10]
################################################################################
data = read.dta("data-1cigar.dta")
#data[data[,2]==10,2]=9
Y = as.vector(data[,1]); T = as.vector(data[,2])
Xd = as.matrix(data[,3:46]); Xc = as.matrix(data[,47:53])
rm(data)
n = length(Y); d.Xd = ncol(Xd); d.Xc = ncol(Xc); T.levels = sort(unique(T)); d.T = length(T.levels)

################################################################################
## Nonparametric Estimation - 1-cigarette bin on [0,10]
################################################################################
## Series Basis
knots = c(3,1,3,1,3,1,1)
for (j in 1:d.Xc) {
        if (j==1) R=cbind(Xd, bs(knots=quantile(Xc[,j], probs=seq(1/(knots[j]+1),1-1/(knots[j]+1),length.out=knots[j]),type=1),x=Xc[,j]))
        else R=cbind(R,bs(knots=quantile(Xc[,j], probs=seq(1/(knots[j]+1),1-1/(knots[j]+1),length.out=knots[j]),type=1),x=Xc[,j]))
}
d.R=ncol(R)

## Conditional Expectation (e)
Rt = list(NA)
for (t in 1:d.T) Rt[[t]] = tcrossprod(tcrossprod(solve(qr.R(qr(R[T==T.levels[t],])))),R[T==T.levels[t],])
e.hat = function(m,bt,t) return(R%*%(Rt[[t]]%*%m(Y[T==T.levels[t]],bt)))

## Generalized Propensity Score (p)
gps.fit = multinom(T ~ R, maxit=5000)
#apply(gps.fit$fitted.values,2,summary)
#boxplot(data.frame(gps.fit$fitted.values), notch=T)
#write.csv(gps.fit$fitted.values,file="gps.hat-1cigar.csv",row.names=F)
#gps.hat = read.csv(file="gps.hat.csv")
gps.hat = gps.fit$fitted.values

################################################################################
## Estimation - 1-cigarette bin on [0,10]
################################################################################
d.mt = 7; d.bt = 7; taus = c(.9,.75,.5,1,.25,.1,2)
mt.JID.list = list(function(y,b) ((y<b)-.9),function(y,b) ((y<b)-.75),function(y,b) ((y<b)-.5),function(y,b) (y-b),function(y,b) ((y<b)-.25),function(y,b) ((y<b)-.1),function(y,b) (y^2-b))
mt.JID.matx = function(y,bt) return(cbind((y<bt[1])-.9,(y<bt[2])-.75,(y<bt[3])-.5,y-bt[4],(y<bt[5])-.25,(y<bt[6])-.1),y^2-bt[7])
Table1.rownames=c("0","","1","","2","","3","","4","","5","","6","","7","","8","","9","","10","")
Table1raw.colnames=rep(c("Q.9","Q.75","Q.5","Mean","Q.25","Q.1","2ndMom"),3)
Table1.colnames=rep(c("Q.9","Q.75","Q.5","Mean","Q.25","Q.1","SD"),3)

Table1raw = matrix(NA,nrow=d.T*2,ncol=3*d.bt,
                dimnames=list(Table1.rownames,Table1raw.colnames))
index.coef = seq(1,by=2,length.out=d.T)

## Dummy Regression
for (t in 1:d.T) for (i in 1:length(taus)){
        if (taus[i]<1)  fit = summary(rq(Y~as.numeric(T==T.levels[t])-1,tau=taus[i]),se="ker")
        else fit = summary(lm(Y^taus[i]~as.numeric(T==T.levels[t])-1))
        Table1raw[(2*t-1),i] = round(fit$coef[,1]); Table1raw[(2*t),i] = round(fit$coef[,2])
        cat(paste("\nTreatment",t,"of",d.T,"| Parameter",i,"of",length(mt.JID.list),": Done!"))
}
Table1 = Table1raw; colnames(Table1)=Table1.colnames
Table1[,7] = round(sqrt(Table1raw[,7]-(Table1raw[,4])^2)); Table1[index.coef+1,7] = NA

##########################
## IPW Estimator
##########################
interval.bt = cbind(c(rep(1500,d.bt-1),9000000),c(rep(5500,d.bt-1),13000000))
for (t in 1:d.T) for (i in 1:length(mt.JID.list)){
        IPW.JID = optimize(f=function(b) return(abs(sum(m.IPW(m=mt.JID.list[[i]],b,p=gps.hat[,t],t=t)))),
                           interval=interval.bt[i,])
        Table1raw[(2*t-1),d.bt+i] = round(IPW.JID$minimum);
        cat(paste("\nTreatment",t,"of",d.T,"| Parameter",i,"of",length(mt.JID.list),": Done!"))
}

b=c(t(Table1raw[index.coef,(d.bt+1):(2*d.bt)]))
# Using numerical derivative
#G = Gamma(m=mt.JID.matx,d.m=6,b=b,d.bt=6,e=e.hat,h=10)
# Using exact formula: density at q or -1
G.IPW = matrix(0,nrow=d.bt*d.T,ncol=d.mt*d.T); h = 1.06*sd(Y)*n^(-1/5)
for (t in 1:d.T) for (i in 1:d.bt) {
        if (taus[i]<1) G.IPW[(t-1)*d.bt+i,(t-1)*d.bt+i]=mean(((T==T.levels[t])/gps.hat[,t])*dnorm((Y-b[(t-1)*d.bt+i])/h)/h)
        else G.IPW[(t-1)*d.bt+i,(t-1)*d.bt+i]=-1
}
V.b.IPW = V.EIF.JID(m=mt.JID.list,b=b,p=gps.hat,e=e.hat)
SPEB.IPW.JID = solve(G.IPW)%*%V.b.IPW%*%solve(G.IPW)/n
Table1raw[index.coef+1,(d.bt+1):(2*d.bt)] = t(matrix(round(sqrt(diag(SPEB.IPW.JID))),nrow=d.bt,ncol=d.T))
write.csv(Table1raw, file="Table1raw-1cigar.csv")

Table1[index.coef,(d.bt+1):(2*d.bt)]=Table1raw[index.coef,(d.bt+1):(2*d.bt)]
Table1[,14] = round(sqrt(Table1raw[,14]-(Table1raw[,11])^2))
Df = diag(1,d.T*d.bt)
diag(Df[seq(4,by=d.bt,length.out=d.T),seq(7,by=7,length.out=d.T)]) = -c(Table1raw[index.coef,11])/c(Table1[index.coef,14])
diag(Df[seq(7,by=d.bt,length.out=d.T),seq(7,by=7,length.out=d.T)]) = 1/(2*c(Table1[index.coef,14]))

SPEB.IPW.JID.f = crossprod(Df,SPEB.IPW.JID)%*%Df
Table1[index.coef+1,(d.bt+1):(2*d.bt)] = t(matrix(round(sqrt(diag(SPEB.IPW.JID.f))),nrow=d.bt,ncol=d.T))
write.csv(Table1, file="Table1-1cigar.csv")

########################
## EIF Estimator
########################
for (t in 1:d.T) for (i in 1:length(mt.JID.list)){
        EIF.JID = optimize(f=function(b) return(abs(sum(m.EIF(m=mt.JID.list[[i]],b,p=gps.hat[,t],t=t,e=e.hat)))),
                           interval=interval.bt[i,])
        Table1raw[(2*t-1),2*d.bt+i] = round(EIF.JID$minimum);
        cat(paste("\nTreatment",t,"of",d.T,"| Parameter",i,"of",length(mt.JID.list),": Done!"))
}

b=c(t(Table1raw[index.coef,(2*d.bt+1):(3*d.bt)]))
G.EIF = matrix(0,nrow=d.bt*d.T,ncol=d.mt*d.T); h = 1.06*sd(Y)*n^(-1/5)
for (t in 1:d.T) for (i in 1:d.bt) {
        if (taus[i]<1) G.EIF[(t-1)*d.bt+i,(t-1)*d.bt+i]=mean(((T==T.levels[t])/gps.hat[,t])*dnorm((Y-b[(t-1)*d.bt+i])/h)/h)
        else G.EIF[(t-1)*d.bt+i,(t-1)*d.bt+i]=-1
}
V.b.EIF = V.EIF.JID(m=mt.JID.list,b=b,p=gps.hat,e=e.hat)
SPEB.EIF.JID = solve(G.EIF)%*%V.b.EIF%*%solve(G.EIF)/n
Table1raw[index.coef+1,(2*d.bt+1):(3*d.bt)] = t(matrix(round(sqrt(diag(SPEB.EIF.JID))),nrow=d.bt,ncol=d.T))
write.csv(Table1raw, file="Table1raw-1cigar.csv")

Table1[index.coef,(2*d.bt+1):(3*d.bt)]=Table1raw[index.coef,(2*d.bt+1):(3*d.bt)]
Table1[,21] = round(sqrt(Table1raw[,21]-(Table1raw[,18])^2))
Df = diag(1,d.T*d.bt)
diag(Df[seq(4,by=d.bt,length.out=d.T),seq(7,by=7,length.out=d.T)]) = -c(Table1raw[index.coef,18])/c(Table1[index.coef,21])
diag(Df[seq(7,by=d.bt,length.out=d.T),seq(7,by=7,length.out=d.T)]) = 1/(2*c(Table1[index.coef,21]))

SPEB.EIF.JID.f = crossprod(Df,SPEB.EIF.JID)%*%Df
Table1[index.coef+1,(2*d.bt+1):(3*d.bt)] = t(matrix(round(sqrt(diag(SPEB.EIF.JID.f))),nrow=d.bt,ncol=d.T))
write.csv(Table1, file="Table1-1cigar.csv")

## Latex Table
Table1tex = Table1
Table1tex[index.coef+1,] = paste("(",Table1[index.coef+1,],")",sep="")
Table1tex[index.coef+1,7] = "n.a."
colnames(Table1tex)=Table1.colnames
rownames(Table1tex)=Table1.rownames
latex(Table1tex, title="", file="table1-1cigar.tex",
      n.group=c(1:7,8:14,15:21),cgroup=c("DRE","IPWE","EIFE"),
      size="scriptsize",
      caption = "Effect of Maternal Smoking Intensity on Birth Weight",
      insert.bottom="Notes: (i) DRE = Dummy Regression Estimator, IPWE = Inverse Probability Weighting Estimator, EIFE = Eficient Influence Function Estimator; (ii) standard errors in parenthesis."
)


################################################################################
## Figures - 1-cigarette bin on [0,10]
################################################################################
#Table1 = as.matrix(read.csv(file="Table1.csv")[1:(2*d.T),2:22])

## IPW: Figure1
postscript(file="Figure1-IPW-JID-1cigar.ps", horizontal=F)
index.colstoplot = (d.bt+1):(2*d.bt-1)
plotdata.IPW = rep(1:d.T,6)
plotdata.IPW = cbind(plotdata.IPW,c(Table1[index.coef,index.colstoplot]),c(Table1[index.coef+1,index.colstoplot]))
plot(plotdata.IPW[,1:2], axes=F,
       #main="Effect of Maternal Smoking Intensity on Birth Weight \n (Just-identified IPW Estimates)",
       ylab="Birth weight", xlab="Cigarettes-per-day Smoked During Pregnancy")
box(); axis(1, at=1:d.T, lab=Table1.rownames[!(Table1.rownames=="")]); axis(2)
lty.setting=c(2,4,3,1,4,2)
for (i in 1:6) lines(plotdata.IPW[seq(1+(i-1)*d.T,length.out=d.T),1:2], lty=lty.setting[i])
for (i in 1:6) for (j in 1:d.T) lines(c(j,j),c(plotdata.IPW[j+(i-1)*d.T,2]-1.96*plotdata.IPW[j+(i-1)*d.T,3],
                                               plotdata.IPW[j+(i-1)*d.T,2]+1.96*plotdata.IPW[j+(i-1)*d.T,3]))
legend("topright",c("Q.9 & Q.10","Q.75 & Q.25","Median","Mean"),lty=lty.setting[1:4])
dev.off()



