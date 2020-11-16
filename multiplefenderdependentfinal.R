# starting of x trial

x=c(13,220000,0.8,1,15)

step<-0.01

g          <-1
Ecv        <-1140
deflection <-1400

#Weibull parameters
a<-1.7 #k
b<-4.38 #lambda(cm/s)

#Uniform parameters
p<- 8000
q<- 260000

# starting of x trial

while(g>0.00001){
U1   <- cdf(Weibull(a,b),x[1]) #velocity
U2   <- cdf(Uniform(p,q),x[2]) #mass
U21  <- BiCopHfunc1(U1,U2,34,-1.69) #CDF of mass given velocity
U12  <- BiCopHfunc2(U1,U2,34,-1.69) #CDF of vel given mass
U3   <- cdf(Gamma(1.19,1/0.27),x[3]) #angle
U32  <- BiCopHfunc1(U2,U3,33,-1.12) #CDF of angle given mass
U312 <- BiCopHfunc1(U12,U32,5,1.35) #CDF of angle 
U4   <- cdf(Normal(1,0.033),x[4]) #Manufacture factor
U5   <- cdf(Normal(15,5),x[5]) #temperature

#transform to standard normal space
y1 <- qnorm(U1,0,1) 
y2 <- qnorm(U21,0,1)
y3 <- qnorm(U312,0,1)
y4 <- qnorm(U4,0,1)
y5 <- qnorm(U5,0,1)

y<-c(y1,y2,y3,y4,y5)

#Construct Jacobian Matrix
#iterasi mulai dari sini
v     <- x[1]
M     <- x[2]
angles<- x[3]
MF    <- x[4]
temp  <- x[5]

# Construct jacobian matrix
Jyx11 <- pdf(Weibull(a,b),v)/pdf(Normal(0,1),y1)
Jyx22 <- BiCopPDF(U1,U2,family=34,par=-1.69)*pdf(Uniform(p,q),M)/pdf(Normal(0,1),y2)
Jyx33 <- BiCopPDF(U12,U32,family=5,1.35)*BiCopPDF(U2,U3,family=33,-1.12)*pdf(Gamma(1.19,1/0.27),angles)/pdf(Normal(0,1),y3)
Jyx44 <- pdf(Normal(1,0.033),MF)/pdf(Normal(0,1),y4)
Jyx55 <- pdf(Normal(15,5),temp)/pdf(Normal(0,1),y5)

# Find Jyx12
h       <-  1
U1h     <-  cdf(Weibull(a,b),v+h)
U21h    <-  BiCopHfunc1(U1h,U2,family=34,-1.69)
dF21dx1 <-  (U21h-U21)/h
Jyx12   <-  dF21dx1/pdf(Normal(0,1),y2)

# Find Jyx13 
U12h     <- BiCopHfunc2(U1h,U2,family=34,-1.69)
U312h    <- BiCopHfunc1(U12h,U32,family=5,1.35)
df312dx1 <- (U312h-U312)/h
Jyx13    <- df312dx1/pdf(Normal(0,1),y3)

# Find Jyx23 
h        <- 1
U2h      <- cdf(Uniform(p,q),M+h)
U12h     <- BiCopHfunc2(U1,U2h,family=34,-1.69)
U32h     <- BiCopHfunc1(U2h,U3,family=33,-1.12)
U312h    <- BiCopHfunc1(U12h,U32h,family=5,1.35)
df312dx2 <- (U312h-U312)/h
Jyx23    <- df312dx2/pdf(Normal(0,1),y3)

Jyx <- rbind(c(Jyx11,0,0,0,0),c(Jyx12,Jyx22,0,0,0),c(Jyx13,Jyx23,Jyx33,0,0),c(0,0,0,Jyx44,0),c(0,0,0,0,Jyx55))
Jxy <- inv(Jyx)

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
#gradg1 to velocity
hv         <- 0.1
vgrad1     <- v+hv
VFgrad1    <- -0.056*log(0.72*deflection/(0.74*vgrad1*10))+1.2011
ggrad1     <- Ecv*VFgrad1*TF*MF*n-0.5*M*(vgrad1/100)^2*Cm*Ce
gradg1     <- (ggrad1-g)/hv

#gradg2 terhadap mass
hM              <- 1000
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
Grad_G <-  gradgx%*%Jxy
alpha  <- -Grad_G/(sqrt(sum(Grad_G^2)))

gradgx <- c

d1     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y)*(alpha[1])-y[1]
d2     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y)*(alpha[2])-y[2]
d3     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y)*(alpha[3])-y[3]
d4     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y)*(alpha[4])-y[4]
d5     <- (g/(sqrt(sum(Grad_G^2)))+alpha%*%y)*(alpha[5])-y[5]
d      <- c(d1,d2,d3,d4,d5)

ynew   <-  y+step*d

# Find new point xnew
U1new   <- pnorm(ynew[1],0,1)
U21new  <- pnorm(ynew[2],0,1)
U2new   <- BiCopHinv1(U1new,U21new,family=34,-1.69)
U312new <- pnorm(ynew[3],0,1)
U12new  <- BiCopHfunc2(U1new,U2new,family=34,-1.69)
U32new  <- BiCopHinv1(U12new,U312new,family=5,1.35)
U3new   <- BiCopHinv1(U2new,U32new,family=33,-1.12)
U4new   <- pnorm(ynew[4],0,1)
U5new   <- pnorm(ynew[5],0,1)

#Find X

x1new <- (b*(-log(1-U1new))^(1/a))
x2new <- invcdf(uniformdistribution(p,q),U2new)
x3new <- invcdf(gammadistribution(1.19,1/0.27),U3new)
x4new <- qnorm(U4new,1,0.033)
x5new <- qnorm(U5new,15,5)

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

designenergy <- 0.5*Mnew*(vnew/100)^2*Cmnew*Cenew

ci    <- (sqrt(sum(y^2))/sqrt(sum(Grad_G ^2)))*2+10

mold <- 0.5*(sqrt(sum(y^2)))^2+ci*abs(g)
mnew <- 0.5*(sqrt(sum(ynew^2)))^2+ci*abs(gnew)

if(mold<mnew){stop("stop")}
print(g)
#print(mold)
#print(mnew)
y  <- ynew

x<-c(x1new,x2new,x3new,x4new,x5new)
}

beta=(sqrt(sum(ynew^2)))


#alphacorrected
alpha1<-qnorm(U1new,0,1)/beta #velocity
alpha2<-qnorm(U2new,0,1)/beta
alpha3<-qnorm(U3new,0,1)/beta
alpha4<-qnorm(U4new,0,1)/beta
alpha5<-qnorm(U5new,0,1)/beta
alphacorrected<-c(alpha1,alpha2,alpha3,alpha4,alpha5)

sigmahat<- Jxy%*%(Jxy)
diagonal<-diag(sigmahat)
dakar<-sqrt(diagonal)
diagonal<- rbind(c(dakar[1],0,0,0,0),c(0,dakar[2],0,0,0),c(0,0,dakar[3],0,0),c(0,0,0,dakar[4],0),c(0,0,0,0,dakar[5]))

atas<- alpha%*%Jyx%*%diagonal

gamma1<-atas/(sqrt(sum(atas^2)))

#printresult
print(beta)
print(alpha)
print(alphacorrected)
print(designenergy)
print(gamma1)

