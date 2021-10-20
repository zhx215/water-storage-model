# water storage two-season model, runoff first, ET second, with different soil storage volume
rr_list <- c(0.05,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4,5) # S0/P1
s0_d18O_list <- s0_dxs_list <- s2_d18O_list <- s2_dxs_list <- ET1_d18O_list <- ET1_dxs_list <- ET2_d18O_list <- ET2_dxs_list <- 
  vector() # result
f1_list <- f2_list <- spr1_list <- spr2_list <- vector()

# additional
s1_d18O_list <- s1_dxs_list <- s3_d18O_list <- s3_dxs_list <- vector()


for (rr in 1:length(rr_list)){
  smow18 <- 0.0020052
  smow2 <- 155.76*10^(-6)
  k18 <- 0.9723^0.8
  k2 <- 0.9755^0.8
  
  residual_ratio <- rr_list[rr] # S0/P1
  nd1 <- 2/3
  nd2 <- -4
  p1 <- 1
  s0 <- residual_ratio*p1
  et1 <- p1*nd1/(1+nd1)
  et2 <- 0.8
  p2 <- (1+nd2)/nd2*et2
  r1 <- 0.25*p1
  r2 <- 0.25*p2
  
  t1 <- 0.6
  t2 <- 0.6
  h1 <- 0.8
  h2 <- 0.65
  temp1 <- 5+273.15
  temp2 <- 15+273.15
  
  p1_d18O <- seq(-15,0,by=0.2)
  p2_d18O <- seq(-15,0,by=0.2)
  p1_d2H <- seq(-140,20,by=1)
  p2_d2H <- seq(-140,20,by=1)
  
  a18_P_lv <- function(temp){
    exp(1137/temp^2-0.4156/temp-0.00207)
  }
  a18_e_vl <- function(temp){
    1/a18_P_lv(temp)
  }
  a2_P_lv <- function(temp){
    exp(24844/temp^2-76.248/temp+0.05261)
  }
  a2_e_vl <- function(temp){
    1/a2_P_lv(temp)
  }
  
  rayleigh_ET <- function(nd,t,h,delta_initial,alpha,k,smow){
    frac <-  1/(1/nd+1)
    r_initial <- (delta_initial/1000+1)*smow
    ae <- ((1-t)*((alpha*k*r_initial)/(1-h+(1-t)*k*h))+t*r_initial*(1/(1+(1-t)*k*(h/(1-h)))))/r_initial
    delta_ET <- ((delta_initial+1000)*(1-frac)^ae-1000*(1-frac)-delta_initial)/(1-frac-1)*(1-t)+delta_initial*t
    return(c((delta_initial-frac*delta_ET)/(1-frac),delta_ET))
  }
  
  ws_isotope <- function(s0,p1,p2,et1,et2,r1,r2,p1_d,p2_d,s0_d,t1,t2,h1,h2,alpha1,alpha2,k,smow){
    s1_d <- (s0*s0_d+(p1-r1)*p1_d)/(s0+p1-r1)
    nd1_corr <- et1/(s0+p1-r1-et1)
    s2_d <- rayleigh_ET(nd1_corr,t1,h1,s1_d,alpha1,k,smow)[1]
    et1_d <- rayleigh_ET(nd1_corr,t1,h1,s1_d,alpha1,k,smow)[2]
    s3_d <- ((s0+p1-r1-et1)*s2_d+(p2-r2)*p2_d)/(s0+p1-r1-et1+p2-r2)
    nd2_corr <- et2/(s0+p1-r1-et1+p2-r2-et2)
    s4_d <- rayleigh_ET(nd2_corr,t2,h2,s3_d,alpha2,k,smow)[1]
    et2_d <- rayleigh_ET(nd2_corr,t2,h2,s3_d,alpha2,k,smow)[2]
    return(c(s1_d,s2_d,s3_d,s4_d,et1_d,et2_d))
  }
  
  # forcing
  d18O_test <- seq(-30,20,by=0.05)
  d2H_test <- seq(-180,80,by=0.1)
  
  # input
  d18O_wet <- -10
  d18O_id_wet <- which(d18O_test==d18O_wet)
  d2H_wet <- d18O_wet*8+10
  d2H_id_wet <- which(d2H_test==d2H_wet)
  
  d18O_dry <- -6
  d18O_id_dry <- which(d18O_test==d18O_dry)
  d2H_dry <- d18O_dry*8+10
  d2H_id_dry <- which(d2H_test==d2H_dry)
  
  # get S0
  result_d18O <- vector()
  for (m in 1:length(d18O_test)){
    result_d18O <- c(result_d18O,ws_isotope(s0,p1,p2,et1,et2,r1,r2,d18O_wet,d18O_dry,d18O_test[m],t1,t2,h1,h2,a18_e_vl(temp1),a18_e_vl(temp2),k18,smow18)[4])
  }
  residual_d18O <- result_d18O[which.min(abs(result_d18O-d18O_test))]
  offset_d18O <- min(abs(result_d18O-d18O_test))
  
  result_d2H <- vector()
  for (m in 1:length(d2H_test)){
    result_d2H <- c(result_d2H,ws_isotope(s0,p1,p2,et1,et2,r1,r2,d2H_wet,d2H_dry,d2H_test[m],t1,t2,h1,h2,a2_e_vl(temp1),a2_e_vl(temp2),k2,smow2)[4])
  }
  residual_d2H <- result_d2H[which.min(abs(result_d2H-d2H_test))]
  offset_d2H <- min(abs(result_d2H-d2H_test))
  
  # print(c(offset_d18O,offset_d2H))
  
  # using S0 to get ET1, ET2, and S2
  et1_d18O <- ws_isotope(s0,p1,p2,et1,et2,r1,r2,d18O_wet,d18O_dry,residual_d18O,t1,t2,h1,h2,a18_e_vl(temp1),a18_e_vl(temp2),k18,smow18)[5]
  et2_d18O <- ws_isotope(s0,p1,p2,et1,et2,r1,r2,d18O_wet,d18O_dry,residual_d18O,t1,t2,h1,h2,a18_e_vl(temp1),a18_e_vl(temp2),k18,smow18)[6]
  s2_d18O <- ws_isotope(s0,p1,p2,et1,et2,r1,r2,d18O_wet,d18O_dry,residual_d18O,t1,t2,h1,h2,a18_e_vl(temp1),a18_e_vl(temp2),k18,smow18)[2]

  et1_d2H <- ws_isotope(s0,p1,p2,et1,et2,r1,r2,d2H_wet,d2H_dry,residual_d2H,t1,t2,h1,h2,a2_e_vl(temp1),a2_e_vl(temp2),k2,smow2)[5]
  et2_d2H <- ws_isotope(s0,p1,p2,et1,et2,r1,r2,d2H_wet,d2H_dry,residual_d2H,t1,t2,h1,h2,a2_e_vl(temp1),a2_e_vl(temp2),k2,smow2)[6]
  s2_d2H <- ws_isotope(s0,p1,p2,et1,et2,r1,r2,d2H_wet,d2H_dry,residual_d2H,t1,t2,h1,h2,a2_e_vl(temp1),a2_e_vl(temp2),k2,smow2)[2]

  # additional
  s1_d18O <- ws_isotope(s0,p1,p2,et1,et2,r1,r2,d18O_wet,d18O_dry,residual_d18O,t1,t2,h1,h2,a18_e_vl(temp1),a18_e_vl(temp2),k18,smow18)[1]
  s3_d18O <- ws_isotope(s0,p1,p2,et1,et2,r1,r2,d18O_wet,d18O_dry,residual_d18O,t1,t2,h1,h2,a18_e_vl(temp1),a18_e_vl(temp2),k18,smow18)[3]
  
  s1_d2H <- ws_isotope(s0,p1,p2,et1,et2,r1,r2,d2H_wet,d2H_dry,residual_d2H,t1,t2,h1,h2,a2_e_vl(temp1),a2_e_vl(temp2),k2,smow2)[1]
  s3_d2H <- ws_isotope(s0,p1,p2,et1,et2,r1,r2,d2H_wet,d2H_dry,residual_d2H,t1,t2,h1,h2,a2_e_vl(temp1),a2_e_vl(temp2),k2,smow2)[3]
  
  # get the result
  s0_d18O_list <- c(s0_d18O_list,residual_d18O)
  s0_dxs_list <- c(s0_dxs_list,residual_d2H-8*residual_d18O)
  s1_d18O_list <- c(s1_d18O_list,s1_d18O)
  s1_dxs_list <- c(s1_dxs_list,s1_d2H-8*s1_d18O)
  s2_d18O_list <- c(s2_d18O_list,s2_d18O)
  s2_dxs_list <- c(s2_dxs_list,s2_d2H-8*s2_d18O)
  s3_d18O_list <- c(s3_d18O_list,s3_d18O)
  s3_dxs_list <- c(s3_dxs_list,s3_d2H-8*s3_d18O)
  ET1_d18O_list <- c(ET1_d18O_list,et1_d18O)
  ET1_dxs_list <- c(ET1_dxs_list,et1_d2H-8*et1_d18O)
  ET2_d18O_list <- c(ET2_d18O_list,et2_d18O)
  ET2_dxs_list <- c(ET2_dxs_list,et2_d2H-8*et2_d18O)
  
  f1_list <- c(f1_list,(s0+p1-r1-et1)/(s0+p1-r1))
  f2_list <- c(f2_list,(s0+p1-r1-et1+p2-r2-et2)/(s0+p1-r1-et1+p2-r2))
  spr1_list <- c(spr1_list,s0/(s0+p1-r1))
  spr2_list <- c(spr2_list,(s0+p1-r1-et1)/(s0+p1-r1-et1+p2-r2))
}

# plot
par(xpd=T,mar=c(38.5,4.5,1.1,1.7))
plot(NA,NA,xlim=c(0.05,5),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n",xaxs="i",yaxs="i",log="x")
axis(1,c(0.05,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4,5),labels=rep("",13),cex.axis=1)
axis(2,seq(0,1,0.5),labels=c("0","0.5","1"),cex.axis=1)
lines(rr_list,f1_list);points(rr_list,f1_list,cex=1,bg="white",pch=21,lwd=1)
lines(rr_list,f2_list,lty=2);points(rr_list,f2_list,cex=1,bg="white",pch=21,lwd=1)
text(0.055,0.62,expression(F[1]))
text(0.055,0.22,expression(F[2]))
text(0.057,0.85,"(a)")

par(new=T,xpd=T,mar=c(32.1,4.5,7.5,1.7))
plot(NA,NA,xlim=c(0.05,5),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n",xaxs="i",yaxs="i",log="x")
axis(1,c(0.05,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4,5),labels=rep("",13),cex.axis=1)
axis(2,seq(0,1,0.5),labels=c("0","0.5","1"),cex.axis=1)
lines(rr_list,spr1_list);points(rr_list,spr1_list,cex=1,bg="white",pch=21,lwd=1)
lines(rr_list,spr2_list,lty=2);points(rr_list,spr2_list,cex=1,bg="white",pch=21,lwd=1)
text(0.055,0.2,expression(X[1]))
text(0.055,0.6,expression(X[2]))
text(0.057,0.85,"(b)")

par(new=T,xpd=T,mar=c(18.1,4.5,14.1,1.7))
plot(NA,NA,xlim=c(0.05,5),ylim=c(-15,10),xlab="",ylab=expression(paste(delta^18,"O (","\u2030",")")),xaxt="n",yaxt="n",xaxs="i",yaxs="i",log="x")
axis(1,c(0.05,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4,5),labels=rep("",13),cex.axis=1)
axis(2,seq(-15,10,5),cex.axis=1)
lines(rr_list,rep(-10,13),lwd=1.5)
lines(rr_list,rep(-6,13),lwd=1.5,lty=2)
lines(rr_list,s0_d18O_list,col="brown3",lwd=1.5);points(rr_list,s0_d18O_list,cex=1.3,bg="brown3",pch=24,lwd=1)
lines(rr_list,s2_d18O_list,col="brown3",lwd=1.5,lty=2);points(rr_list,s2_d18O_list,cex=1.3,bg="brown3",pch=24,lwd=1)
lines(rr_list,ET1_d18O_list,col="chartreuse3",lwd=1.5);points(rr_list,ET1_d18O_list,cex=1.3,bg="chartreuse3",pch=21,lwd=1)
lines(rr_list,ET2_d18O_list,col="chartreuse3",lwd=1.5,lty=2);points(rr_list,ET2_d18O_list,cex=1.3,bg="chartreuse3",pch=21,lwd=1)
lines(rr_list,s1_d18O_list,col="darkgoldenrod1",lwd=1.5);points(rr_list,s1_d18O_list,cex=1.3,bg="darkgoldenrod1",pch=23,lwd=1)
lines(rr_list,s3_d18O_list,col="darkgoldenrod1",lwd=1.5,lty=2);points(rr_list,s3_d18O_list,cex=1.3,bg="darkgoldenrod1",pch=23,lwd=1)
text(5.6,-10,expression(delta[P[1]]))
text(5.6,-6,expression(delta[P[2]]))
text(0.057,8.5,"(c)")

legend("topright",c(expression(delta[I[1]]),expression(delta[I[2]]),expression(delta[ET[1]]),expression(delta[ET[2]]),expression(delta[S[0]]),expression(delta[S[1]])),
       col=c("darkgoldenrod1","darkgoldenrod1","chartreuse3","chartreuse3","brown3","brown3"),
       lwd=1.5,lty=c(1,2,1,2,1,2),ncol=3,text.width=0.08)
points(0.95,8,cex=1.3,bg="darkgoldenrod1",pch=23,lwd=1);points(0.95,5.95,cex=1.3,bg="darkgoldenrod1",pch=23,lwd=1)
points(1.745,8,cex=1.3,bg="chartreuse3",pch=21,lwd=1);points(1.745,5.95,cex=1.3,bg="chartreuse3",pch=21,lwd=1)
points(3.2,8,cex=1.3,bg="brown3",pch=24,lwd=1);points(3.2,5.95,cex=1.3,bg="brown3",pch=24,lwd=1)

par(new=T,xpd=T,mar=c(4.1,4.5,28.1,1.7))
plot(NA,NA,xlim=c(0.05,5),ylim=c(-35,20),xlab=expression(paste(S[0],"/",P[1])),ylab=expression(paste("d-excess (","\u2030",")")),xaxt="n",yaxt="n",xaxs="i",yaxs="i",log="x")
axis(1,c(0.05,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4,5),labels=c(0.05,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4,5),cex.axis=1)
axis(2,seq(-30,20,10),cex.axis=1)
lines(rr_list,rep(10,13),lwd=1.5)
lines(rr_list,s0_dxs_list,col="brown3",lwd=1.5);points(rr_list,s0_dxs_list,cex=1.3,bg="brown3",pch=24,lwd=1)
lines(rr_list,s2_dxs_list,col="brown3",lwd=1.5,lty=2);points(rr_list,s2_dxs_list,cex=1.3,bg="brown3",pch=24,lwd=1)
lines(rr_list,ET1_dxs_list,col="chartreuse3",lwd=1.5);points(rr_list,ET1_dxs_list,cex=1.3,bg="chartreuse3",pch=21,lwd=1)
lines(rr_list,ET2_dxs_list,col="chartreuse3",lwd=1.5,lty=2);points(rr_list,ET2_dxs_list,cex=1.3,bg="chartreuse3",pch=21,lwd=1)
lines(rr_list,s1_dxs_list,col="darkgoldenrod1",lwd=1.5);points(rr_list,s1_dxs_list,cex=1.3,bg="darkgoldenrod1",pch=23,lwd=1)
lines(rr_list,s3_dxs_list,col="darkgoldenrod1",lwd=1.5,lty=2);points(rr_list,s3_dxs_list,cex=1.3,bg="darkgoldenrod1",pch=23,lwd=1)
text(0.057,16.7,"(d)")
