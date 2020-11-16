# Mean
zero <- c(3,2,1,0,0)
y0   <- zero #in standard space not probability
e1<-1

step     <- 0.002

#Change here
Ecv        <- 1140
deflection <- 1400

while (e1>0.01){
  # CDF in normal space
  U1   <- pnorm(y0[1],0,1)
  U2   <- pnorm(y0[2],0,1)  
  U3   <- pnorm(y0[3],0,1)
  U4   <- pnorm(y0[4],0,1)
  U5   <- pnorm(y0[5],0,1)
  
  #real value
  x1 <- (4.38*(-log(1-U1))^(1/1.7)) #velocity
  x2 <- invcdf(uniformdistribution(8000,260000),U2) #displacement
  x3 <- invcdf(gammadistribution(alpha=1.19,beta=1/0.27),U3) #angle
  x4 <- qnorm(U4,1,0.033)
  x5 <- qnorm(U5,15,5)
  
  #Jacobian
  Jyx11 <- pdf(Weibull(1.7,4.38),x1)/pdf(Normal(0,1),y0[1])
  Jyx22 <- pdf(Uniform(8000,260000),x2)/pdf(Normal(0,1),y0[2])
  Jyx33 <- pdf(Gamma(1.19,1/0.27),x3)/pdf(Normal(0,1),y0[3])
  Jyx44 <- pdf(Normal(1,0.033),x4)/pdf(Normal(0,1),y0[4])
  Jyx55 <- pdf(Normal(15,5),x5)/pdf(Normal(0,1),y0[5])
  
  Jyx <- rbind(c(Jyx11,0,0,0,0),c(0,Jyx22,0,0,0),c(0,0,Jyx33,0,0),c(0,0,0,Jyx44,0),c(0,0,0,0,Jyx55))
  Jxy <- inv(Jyx)
  
  v      <- x1
  M      <- x2
  angles <- x3
  MF   <- x4
  temp <- x5
  
  #limit state function
  
  Lpp        <- ((M)^0.3858)*3.279
  D          <- ((M)^0.2321)*0.9519
  B          <- ((M)^0.33339)*0.8528
  Cb         <- M/(Lpp*B*D*1.02)
  R          <- sqrt((Lpp/4)^2+(B/2)^2)
  K          <- (0.19*Cb+0.11)*Lpp
  TF         <- -3*10^-6*(temp^3)+0.0002*(temp^2)-0.0058*temp+1.0533
  gamma      <- 90-angles-asin(B/(2*R))
  Ce         <- (K^2+R^2*(cosd(gamma))^2)/(K^2+R^2)
  Cm         <- 1+pi*D/(2*Cb*B)
  VF         <- -0.056*log(0.72*deflection/(0.74*v*10))+1.2011
  j          <- (0.75*Lpp)/14+1
  k          <- 2.8495*log(Lpp)-11.996
  n          <- j*exp(-k*angles)+(1.43*log(Lpp)-5.76)
  g          <- Ecv*VF*TF*MF*n-0.5*M*(v/100)^2*Cm*Ce
  
  #gradien
  #gradg1 terhadap velocity
  hv         <- 0.1
  vgrad1     <- v+hv
  VFgrad1    <- -0.056*log(0.72*deflection/(0.74*vgrad1*10))+1.2011
  ggrad1     <- Ecv*VFgrad1*TF*MF*n-0.5*M*(vgrad1/100)^2*Cm*Ce
  gradg1     <- (ggrad1-g)/hv
  
  #gradg2 terhadap mass
  hM              <- 10
  Mgrad2          <- M+hM
  Lppgrad2        <- ((Mgrad2)^0.3858)*3.279
  Dgrad2          <- ((Mgrad2)^0.2321)*0.9519
  Bgrad2          <- ((Mgrad2)^0.33339)*0.8528
  Cbgrad2         <- Mgrad2 /(Lppgrad2*Bgrad2*Dgrad2*1.02)
  Rgrad2          <- sqrt((Lppgrad2/4)^2+(Bgrad2/2)^2)
  Kgrad2          <- (0.19*Cbgrad2+0.11)*Lppgrad2
  gammagrad2      <- 90-angles-asin(Bgrad2/(2*Rgrad2))
  Cegrad2         <- (Kgrad2^2+Rgrad2^2*(cosd(gamma))^2)/(Kgrad2^2+Rgrad2^2)
  Cmgrad2         <- 1+pi*Dgrad2/(2*Cbgrad2*Bgrad2)
  jgrad2          <- (0.75*Lppgrad2)/14+1
  kgrad2          <- 2.8495*log(Lppgrad2)-11.996
  ngrad2          <- jgrad2*exp(-kgrad2*angles)+(1.43*log(Lppgrad2)-5.76)
  ggrad2          <- Ecv*VF*TF*MF*ngrad2-0.5*Mgrad2 *(v/100)^2*Cmgrad2*Cegrad2
  gradg2          <- (ggrad2-g)/hM
  
  #gradg3 terhadap angle
  hang            <- 0.1
  anglesgrad3     <- angles+hang
  gammagrad3      <- 90-anglesgrad3 -asin(B/(2*R))
  Ceanglesgrad3   <- (K^2+R^2*(cosd(gammagrad3))^2)/(K^2+R^2)
  nanglesgrad3    <- j*exp(-k*anglesgrad3)+(1.43*log(Lpp)-5.76)
  ggrad3          <- Ecv*VF*TF*MF*nanglesgrad3 -0.5*M*(v/100)^2*Cm*Ceanglesgrad3
  gradg3          <- (ggrad3-g)/hang
  
  #gradg4 terhadap MF
  hMF        <- 0.1
  MFgrad4    <- MF+hMF
  ggrad4     <- Ecv*VF*TF*MFgrad4*n-0.5*M*(v/100)^2*Cm*Ce
  gradg4     <- (ggrad4-g)/hMF
  
  #gradg5 terhadap temperature
  htemp      <- 0.1
  tempgrad5  <- temp+htemp
  TFgrad5    <- -3*10^-6*(tempgrad5^3)+0.0002*(tempgrad5^2)-0.0058*tempgrad5+1.0533
  ggrad5     <- Ecv*VF*TFgrad5 *MF*n-0.5*M*(v/100)^2*Cm*Ce
  gradg5     <- (ggrad5-g)/htemp
  
  gradgx <- c(gradg1,gradg2,gradg3,gradg4,gradg5)
  
  Grad_G <- gradgx%*%Jxy
  alpha  <- -Grad_G/(sqrt(sum(Grad_G^2)))
  
  d1     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y0)*(alpha[1])-y0[1]
  d2     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y0)*(alpha[2])-y0[2]
  d3     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y0)*(alpha[3])-y0[3]
  d4     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y0)*(alpha[4])-y0[4]
  d5     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y0)*(alpha[5])-y0[5]
  d      <- c(d1,d2,d3,d4,d5)
  
  y0new    <- y0+step*d
  
  U1new    <- pnorm(y0new[1],0,1)
  U2new    <- pnorm(y0new[2],0,1)  
  U3new    <- pnorm(y0new[3],0,1)
  U4new    <- pnorm(y0new[4],0,1)
  U5new    <- pnorm(y0new[5],0,1)
  
  x1new   <- (4.38*(-log(1-U1new))^(1/1.7)) #velocity
  x2new   <- invcdf(uniformdistribution(8000,260000),U2new) #displacement
  x3new   <- invcdf(gammadistribution(alpha=1.19,beta=1/0.27),U3new) #angle
  x4new   <- qnorm(U4new,1,0.033)
  x5new   <- qnorm(U5new,15,5)
  
  vnew      <- x1new
  Mnew      <- x2new
  anglesnew <- x3new
  MFnew     <- x4new
  tempnew   <- x5new
  
  
  Lppnew        <- ((Mnew)^0.3858)*3.279
  Dnew          <- ((Mnew)^0.2321)*0.9519
  Bnew          <- ((Mnew)^0.33339)*0.8528
  Cbnew         <- Mnew/(Lppnew*Bnew*Dnew*1.02)
  Rnew          <- sqrt((Lppnew/4)^2+(Bnew/2)^2)
  Knew          <- (0.19*Cbnew+0.11)*Lppnew
  TFnew         <- -3*10^-6*(tempnew^3)+0.0002*(tempnew^2)-0.0058*tempnew+1.0533
  gammanew      <- 90-anglesnew-asin(Bnew/(2*Rnew))
  Cenew         <- (Knew^2+Rnew^2*(cosd(gammanew))^2)/(Knew^2+Rnew^2)
  Cmnew         <- 1+pi*Dnew/(2*Cbnew*Bnew)
  VFnew         <- -0.056*log(0.72*deflection/(0.74*vnew*10))+1.2011
  jnew          <- (0.75*Lppnew)/14+1
  knew          <- 2.8495*log(Lppnew)-11.996
  nnew          <- jnew*exp(-knew*anglesnew)+(1.43*log(Lppnew)-5.76)
  gnew          <- Ecv*VFnew*TFnew*MFnew*nnew-0.5*Mnew*(vnew/100)^2*Cmnew*Cenew
  
  mold <- 0.5*(sqrt(sum(y0^2)))^2+10*abs(g)
  mnew <- 0.5*(sqrt(sum(y0new^2)))^2+10*abs(gnew)
  
  if(mold<mnew){stop("mold lebih kecil")}
  
  e1<-g
  print(g)
  
  y0 <- y0new}

beta <-(sqrt(sum(y0new^2)))
print(beta)

Designenergy<-0.5*Mnew*(vnew/100)^2*Cmnew*Cenew
print(Designenergy)
