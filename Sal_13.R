Sal_13<-read.csv(file.choose(), head=TRUE)
# Use the next couple of lines if detrended series desired.
# Use the next couple of lines if detrended series desired.
detrSal13<-lm(Sal_13$log_data[1:129]~seq(1,129,1))
N_Sal_13<-exp(detrSal13$residuals)

#Use this line if a natural series desired
#N_Sal_13<-round(Sal_13$mean_copy_fill[1:129])

#Time series plot
par(pty="m")
Sal_13_ts<-plot(N_Sal_13[1:129], ty="b", lwd=2,
                ylab="copy number", xlab="time (Julian day)")
legend("topright", inset=0.015, legend=("Salinas 2013"))


#Phase plot
par(pty="s")
Sal_13_Phpl<-plot(log(N_Sal_13[1:128]),log(N_Sal_13[2:129]),
                  ty="b", lwd=2, xlab="log copy number, time=t",
                  ylab="log copy number, time=t+1")
legend("topleft", inset=0.015, legend=("Salinas 2013"))
lines(log(N_Sal_13[1:129]),log(N_Sal_13[1:129]), col="red3", lty=2)


#ACF function
par(pty="s")
lagmax<-round(0.25*length(N_Sal_13))#Max lag value for acf and MI work
N_Sal_13_acf<-acf(N_Sal_13, type="correlation", plot=TRUE, lag.max=lagmax,
                  ylab="auto-correltation function", xlab="time lag index",
                  main="")
legend("topright", inset=0.015, legend=("Salinas 2013"))

#Average mutual information is supposed to get us to the embedding time delay
Sal_13_AMI<-mutual(N_Sal_13, lag.max = lagmax, main="")
legend("topright", inset=0.015, legend=("Salinas 2013"))

#calculate the average mutual information with lag
# Kcs_13_AMI2<-mutualInformation(N_Kcs_13, lag.max = lagmax, n.partitions = NULL,
#            units = "Nats", do.plot = TRUE, main="")

#According to Huffaker et al p179, the first minima of the AMI, ACF and FNN functions tell
#you first estimtates of the embedding delay (d), Theiler window (tw or t), and embedding dimension (m)
#Suggested code for finding embedding dimension. Adapted from code on p68 of Huffaker et al
#These values are after a bit of experimentation
#plot(Kcs_13_AMI) #to see where first minimum
m.max<-8
d<-5
tw<-3
rt<-10
eps<-sd(N_Sal_13)
fn<-false.nearest(N_Sal_13,m.max,d,tw,rt,eps=eps)
fn
par(pty="s")
plot(fn, main="")
legend("topleft", inset=0.015, legend=("Salinas 2013"))

#Alternative approach for finding m using nonlinearTseries package
Sal_13_ED<-estimateEmbeddingDim(
  N_Sal_13,  number.points = length(N_Sal_13),
  time.lag = 1, max.embedding.dim = 15,
  threshold = 0.95, max.relative.change = 0.1,
  do.plot = TRUE, main = "", xlab = "dimension (d)",
  ylab = "E1(d) & E2(d)", ylim = NULL, xlim = NULL, std.noise=1
)
legend("right",inset=0.015, legend=("Salinas 2013"))

Sal_13_MCLE<-lyap_k(N_Sal_13,m=3,d=5,t=3,k=15,ref=100,s=60,eps=eps)
plot(Sal_13_MCLE) # to identify the linear portion of the MCLE response
locator() #to capture the end points of the linear seciton (See p89 in Huffaker et al)
lyap(Sal_13_MCLE,24,42)
Sal_13lr<-Sal_13_MCLE[24:42]
Sal_13lt<-seq(1,length(Sal_13lr),1)
Sal_13lm_lyap<-lm(Sal_13lr~Sal_13lt)
summary(Sal_13lm_lyap)
Sal_13lm_lyap$coefficients


#Alternative function from nonlinearTseries
#par(pty="m")
#maxlam<-maxLyapunov(
#  N_Sal_13,
#  min.embedding.dim = 2,
#  max.embedding.dim = 7,
#  time.lag = 1,
#  radius = eps,
#  theiler.window = 7,
#  min.neighs = 5,
#  min.ref.points = 30,
#  max.time.steps = 15,
#  number.boxes = NULL,
#  sampling.period = 1,
#  do.plot = TRUE,
# )

#Now calculate the running entropy value
Sal_13ents<-vector(length=120)
range<-range(N_Sal_13)
Sal_13dR1<-discretize(N_Sal_13[1:10], numBins = 15, r=range)
Sal_13ents[1]<-entropy(Sal_13dR1, method = "ML")

for (i in seq(1,119,1)){
  dR_temp<-discretize(N_Sal_13[1:(10+i)], numBins = 15, r=range)
  Sal_13ents[i+1]<-entropy(dR_temp, method = "ML")
}
#plot the entropy with the original time series
#plot the entropy with the original time series
plot.new()
par(mar=c(5,5,2,5))
plot(seq(1:120),Sal_13ents, ylim=c(0,max(Sal_13ents+0.5)), ty="l", 
     lwd=2, ylab="Entropy (nats)", xlab="time point")
par(new = T)
plot(seq(10,129,1),N_Sal_13[10:129], ty="b", col="darkgray",
     lty=2, pch=19, axes=F, xlab=NA, ylab=NA)
axis(side=4)
mtext(side=4, line=3, 'copy number')
legend("topright", inset=0.015,
       legend=c("Entropy (nats)", "copy number"),
       lty=c(1,2), lwd=c(2,1), pch=c(NA, 19), col=c("black", "darkgray"))

#Surrogate tests for non-linear dependence
# set.seed(1)
# S_Sal_13<-Srho.test.ts(N_Sal_13, B=50, lag.max = lagmax, ci.type = c("perm"),
#                       quant = 0.95)
#S_Sal_13
#plot(S_Sal_13, cex.lab=1.2,cex.axis=1.2, lwd=2)

#Alternative method from nonlinearTseries package
st_Sal_13 = surrogateTest(N_Sal_13,significance = 0.02, one.sided = F,
                          FUN = timeAsymmetry, do.plot=F)
plot.new()
par(pty="m")
barplot(c(st_Sal_13$data.statistic, st_Sal_13$surrogates.statistics),
        col=c("red3",rep("lightgray",100)), ylab="test statistic value")
abline(h=sd(st_Sal_13$surrogates.statistics), col="red", lty=2)
abline(h=-sd(st_Sal_13$surrogates.statistics), col="red", lty=2)
legend("bottom", inset=0.015, legend="Salinas, 2013")

#+---------------------+
#Regression model to fit AR model to growth rates
R_Sal_13<-diff(Sal_13$log_data[1:129]) #calculate instantaneous growth rates as differences
#ditch last two values which are NA

length(R_Sal_13)

#The first minimum of the AMI function indicates the embedding delay, so d=4 is a good
#first shot at constructing the linear model of Rt = f(Xt-d)+e
#For Rt, points at Xt, Xt-1, Xt-2 and Xt-3 are needed
lmR_Sal_13<-lm(R_Sal_13[4:128]~Sal_13$log_data[4:128]+Sal_13$log_data[3:127]+Sal_13$log_data[2:126]+Sal_13$log_data[1:125])
lmR_Sal_13
summary(lmR_Sal_13)
plot(lmR_Sal_13$fitted.values, ty="b", #diagnostic plot to take a look at the model fit
     ylim=c(min(R_Sal_13),max(R_Sal_13)))
lines(R_Sal_13[4:128] ,col="darkgray")

#only the first lag is significant - relate this back to the results from the surrogate testing
#refit the model with Lag-1 and -4
lmR_Sal_13<-lm(R_Sal_13[4:128]~Sal_13$log_data[4:128]+Sal_13$log_data[1:125])
lmR_Sal_13
summary(lmR_Sal_13)
plot(lmR_Sal_13$fitted.values, ty="b", #diagnostic plot to take a look at the model fit
     ylim=c(min(R_Sal_13),max(R_Sal_13)))
lines(R_Sal_13[4:128] ,col="darkgray")

#+---------------------------+
#plot the model with the data in good quality
plot(lmR_Sal_13$fitted.values, ty="b", ylim=c(min(R_Sal_13),max(R_Sal_13)),
     ylab="log growth rate", xlab="time (days)")
lines(R_Sal_13[4:128] ,col="darkgray", lwd=2)
abline(h=0, lwd=1, lty=2, col="red3") #add this because if the series is stationary it should center here
legend("bottomleft", inset=0.015, lty=c(1,1), lwd=c(1,2), pch=c(1,NA),
       col=c("black", "darkgrey"), legend=c("fitted values", "observed"))

#+---------------------------+
# Use mean as model and compare with previous
R_Sal_13mn<-mean(R_Sal_13)
R_Sal_13mnv<-rep(R_Sal_13mn,128)
lmR_Sal_13mn<-lm(R_Sal_13~R_Sal_13mnv)

lmR_Sal_13mn
summary(lmR_Sal_13mn)

plot(lmR_Sal_13mn$fitted.values, ty="b", #diagnostic plot to take a look at the model fit
     ylim=c(min(R_Sal_13),max(R_Sal_13)))
lines(R_Sal_13[4:128] ,col="darkgray")

aov_out<-aov(lmR_Sal_13)
aov_out
aov_out2<-aov(lmR_Sal_13mn)
aov_out2

#+-----------------------------+
#Investigation of entropy of binary version of R series
Sal_13Rshan<-vector(length=length(R_Sal_13))
Sal_13Rbits<-as.numeric(R_Sal_13>0)                  
for (i in seq(1,length(Sal_13Rbits),1)){
  Sal_13Rshan[i]<--((cumsum(Sal_13Rbits[1:i])/i)*log((cumsum(Sal_13Rbits[1:i])/i),2))+
    ((1-(cumsum(Sal_13Rbits[1:i])/i))*log((1-(cumsum(Sal_13Rbits[1:i])/i)),2))
}

#+------------------+
#Construct a bootstrap resampling experiment to illustrate departure from
#randomness in observed series

#+---------------------------------+
# This section calculates the entropy of the real series 
Sal_13cum<-cumsum(Sal_13Rbits)
Sal_13p<-vector(length=length(Sal_13Rbits))

#Convert to proportion
for (i in seq(1:128)){
  Sal_13p[i]<-Sal_13cum[i]/i
}
Sal_13q<-1-Sal_13p # q = (1-p)

# plot(Sal_13p, ty="l", ylim=c(0,1), col="red")
#  lines(Sal_13q, col="blue")
#  lines(dum1p)

log2Sal_13p<-log(Sal_13p,2) #calculate log probability
log2Sal_13q<-log(Sal_13q,2) #omit ,2 in log if nats wanted

Shan_Sal_13<--((Sal_13p*log2Sal_13p)+(Sal_13q*log2Sal_13q)) 

#   plot(shan_dum1, ty="l", col="blue")
#    lines(Sal_13Rshan, col="red")


#+------------------+
#This builds the re-sampled data
Sal_13Rbits_mat<-matrix(ncol=1000, nrow=length(Sal_13Rbits))
cede<-1
set.seed(cede)
for (j in 1:1000){     
  Sal_13Rbits_mat[,j]<-sample(Sal_13Rbits)
  cede=cede+1 
}
Sal_13sum_mat<-matrix(nrow=128, ncol=1000)
for (i in 1:1000){
  Sal_13sum_mat[,i]<-cumsum(Sal_13Rbits_mat[,i])
}
# plot(seq(1:128),sum_mat[,1], ty="l")
# for (i in seq(2,999,1)){
#     lines(sum_mat[,i])
#    }
#     lines(cumsum(Sal_13Rbits), lwd=3, col="red")

Sal_13p_mat<-matrix(nrow=128, ncol=1000)
Sal_13q_mat<-matrix(nrow=128, ncol=1000)

for (j in 1:1000){
  for (i in 1:128){
    Sal_13p_mat[i,j]<-Sal_13sum_mat[i,j]/i
  }
}
Sal_13q_mat<-1-Sal_13p_mat

Sal_13log2pmat<-log(Sal_13p_mat,2)
Sal_13log2qmat<-log(Sal_13q_mat,2)
#  qShan[,j]<-q_mat[,j]*log2pmat[,j]
# pShan[,j]<-p_mat[,j]*log2pmat[,j]
Sal_13Shan_ent<-matrix(nrow=128, ncol=1000)
for (j in 1:1000){
  Sal_13Shan_ent[,j]<- -((Sal_13p_mat[,j]*Sal_13log2pmat[,j])+(Sal_13q_mat[,j]*Sal_13log2qmat[,j]))
}
Sal_13var<-vector(length=128)
for (i in seq(2,129,1)){
  Sal_13var[i-1]<-var(Sal_13p[1:i])
}

plot(Sal_13var,Shan_Sal_13, ty="b", col="blue")

#+-----------------------------+
#Diagnostic plot


#+------------------+
#This builds the re-sampled data
Sal_13Rbits_mat<-matrix(ncol=1000, nrow=length(Sal_13Rbits))
cede<-1
set.seed(cede)
for (j in 1:1000){     
  Sal_13Rbits_mat[,j]<-sample(Sal_13Rbits)
  cede=cede+1 
}
Sal_13sum_mat<-matrix(nrow=128, ncol=1000)
for (i in 1:1000){
  Sal_13sum_mat[,i]<-cumsum(Sal_13Rbits_mat[,i])
}
# plot(seq(1:128),sum_mat[,1], ty="l")
# for (i in seq(2,999,1)){
#     lines(sum_mat[,i])
#    }
#     lines(cumsum(Sal_13Rbits), lwd=3, col="red")

Sal_13p_mat<-matrix(nrow=128, ncol=1000)
Sal_13q_mat<-matrix(nrow=128, ncol=1000)

for (j in 1:1000){
  for (i in 1:128){
    Sal_13p_mat[i,j]<-Sal_13sum_mat[i,j]/i
  }
}
Sal_13q_mat<-1-Sal_13p_mat

Sal_13log2pmat<-log(Sal_13p_mat,2)
Sal_13log2qmat<-log(Sal_13q_mat,2)
#  qShan[,j]<-q_mat[,j]*log2pmat[,j]
# pShan[,j]<-p_mat[,j]*log2pmat[,j]
Sal_13Shan_ent<-matrix(nrow=128, ncol=1000)
for (j in 1:1000){
  Sal_13Shan_ent[,j]<- -((Sal_13p_mat[,j]*Sal_13log2pmat[,j])+(Sal_13q_mat[,j]*Sal_13log2qmat[,j]))
}


#+-----------------------------+
par(pty="s")
plot(seq(1:128),log(Sal_13Shan_ent[,1]), ylim=c(-0.8,0),  col="darkgray", cex=0.75,
     xlab="time point", ylab="log(Shannon entropy, bits)")
for (i in seq(2,999,1)){
  points(log(Sal_13Shan_ent[,i]), col="darkgray", cex=0.75)
}
 lines(log(Shan_Sal_13), lwd=2, col="red3")
# lines(log(Shan_Sal_14), lwd=2, lty=2, col="red3")
# lines(log(Shan_Sol_13), lwd=2, lty=2, col="red3")
# lines(log(Shan_Sol_14), lwd=2, lty=2, col="red3")
# lines(log(Shan_Gon_13), lwd=2, lty=2, col="red3")
# lines(log(Shan_Gon_14), lwd=2, lty=2, col="red3")
# lines(log(Shan_Kcn_13), lwd=2, lty=2, col="red3")
#lines(log(Shan_Kcs_13), lwd=2, lty=2, col="red3")
#lines(log(Shan_Kcs_14), lwd=2, lty=2, col="lightseagreen")
legend("bottomright", inset=0.015, legend=("Salinas, 2013"))

#Plot Shannon entropy with the data
# ent_time2<-seq(1,128,1)
plot.new()
par(mar=c(5,5,2,5), pty="m")
plot(ent_time2,log(Shan_Sal_13), ylim=c(-0.4,0), ty="l", 
     lwd=2, ylab="log(Entropy (bits))", xlab="time point")
par(new = T)
plot(seq(1:128),R_Sal_13, ty="b", col="darkgray",
     lty=2, pch=19, axes=F, xlab=NA, ylab=NA)
abline(h=0, col="red3", lty=2)
axis(side=4)
mtext(side=4, line=3, 'copy number instantaneous growth rate')
legend("bottomright", inset=0.015,
       legend=("Salinas, 2013"),
       lty=c(1,2), lwd=c(2,1), pch=c(NA, 19), col=c("black", "darkgray"))

#+------------------------+
#Plot the two entropies against each other just for fun - entropy doodles
par(pty="s")
plot(Sal_13ents,log(Shan_Sal_13[9:128]),ty="l", col="red3", lwd=2,
     ylab="log(Entropy (bits))", xlab="Copy Number entropy (nats)",
     main="Salinas, 2013")
#+-----End of bootstrap analysis of daily rate-----+