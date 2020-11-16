% untuk Robert
clear all;
close all;

rng('default');

DWT_cat=xlsread('forsimulation',1,'A2:A27');
disp_cat=xlsread('forsimulation',1,'B2:B27');
LOA_cat=xlsread('forsimulation',1,'C2:C27');
LPP_cat=xlsread('forsimulation',1,'D2:D27');
B_cat=xlsread('forsimulation',1,'E2:E27');
LADEN_cat=xlsread('forsimulation',1,'F2:F27');
CB_cat=xlsread('forsimulation',1,'G2:G27');
TEU_cat=xlsread('forsimulation',1,'H2:H27');

m = 1;
velocity= icdf('Weibull',rand(m,1),4.38,1.7); % in cm/s
angle= icdf('Gamma',rand(m,1), 1.19,0.27);
displacement = 8000+(260000-8000).*rand(m,1); %in ton

% velocity= icdf('Weibull',velcdf(1:m),4.38,1.7); % in cm/s
% angle= icdf('Gamma',angcdf(1:m), 1.19,0.27);
% displacement = 8000+(260000-8000).*masscdf(1:m); %in ton

rateddeflection = xlsread('forsimulation',4,'A1:A23');
ratedenergy = xlsread('forsimulation',4,'B1:B23');
Ratedeflection = 0.72*rateddeflection(17); %mm for cone fender the rated deflection is 72%
Ratedenergy = ratedenergy(17); %kNm

%preallocating
LPP=zeros(m,1);
LOA=zeros(m,1);
B=zeros(m,1);
CB=zeros(m,1);
LADEN=zeros(m,1);
Cm = zeros(m,1);
K=zeros(m,1);
R=zeros(m,1);
gamma=zeros(m,1);

for j = 1:m
    for i = 1:1:25
    if displacement(j,1)<=disp_cat(i) && displacement(j,1)>disp_cat(i+1)
        LOA(j,1) = LOA_cat(i);
        LPP(j,1) = LPP_cat(i);
        B(j,1) = B_cat(i);
        CB(j,1) = CB_cat(i);
        LADEN(j,1) = LADEN_cat(i);
    end
    end
end

%fender coordinate
Fender = 0:14:500;
n = length(Fender);
h = (Ratedeflection/0.72)/1000; %tinggi fender

for i = 1:n
    k(i) = Fender(i); % coordinate x fender
    p(i) = h; % coordinate y fender
    l(i) = h*(1-0.72); % maximum deflection of fender
end
fender(:,1) = k;
fender(:,2) = p;

pusatx = 140; %centre of bow radius
pusaty = h;

tic
for o = 1:m
LOA1 = LOA(o);
Lpp1 = LPP(o);
B1   = B(o);
alpha = angle(o);
z = LOA1*0.25;

Rb= 0.5*(B1/2+LOA1^2/(8*B1));
theta(o)= {0:0.5:(asind(z/Rb)+0.2)};
theta1 = cell2mat(theta(o));
theta1 = sort(theta1,'descend');

clear x y sisa xtrans ytrans ytransmin kapal
for i = 1:length(theta1) % bagian depan kapal
     x(i,1)=pusatx-Rb*sind(theta1(i));
     y(i,1)=Rb - Rb*cosd(theta1(i))+ pusaty;
end

sisa= pusatx:1:(pusatx+LOA1-z); % bagian lurusnya kapal

for i = 1:length(sisa)
    x(i+length(theta1),1) = sisa(i);
    y(i+length(theta1),1) = pusaty;
end

xtrans = pusatx +cosd(alpha)*(x-pusatx)- sind(alpha)*(y-pusaty);
ytrans = pusaty +sind(alpha)*(x-pusatx)+ cosd(alpha)*(y-pusaty);
ytransmin =min(ytrans);


nof = zeros(n,1);
% ngecek fender yang activated pada awalnya
for i = 1:n
    for j = 1:length(xtrans)
    if fender(i,1)-0.2 < xtrans(j,1)+0.2 && fender(i,1)+0.2 > xtrans(j,1) && fender(i,2)>=ytrans(j,1)
       nof(i) = 1;
    end
    end
end

sumnof=sum(nof);
if sumnof == 1
   ytrans = ytrans-0.72*h;
else
   ytrans = ytrans-(ytransmin-0.28*h);
end
%

if any(ytrans(:)<0)
    ytrans = ytrans-min(ytrans);
end
    
% coba itung fender activated
kapal = zeros(length(xtrans),2);
kapal(:,1) = xtrans;
kapal(:,2) = ytrans;

clear index index1 selisihx
for i = 1:n
    for j = 1:length(kapal)
    selisihx(j,i)= abs(fender(i,1)-kapal(j,1));
    minimal = min(selisihx);
      index = min(find(abs(selisihx(:,i)) == minimal(1,i), 1 ));
    end
    index1(i,1)=index;
end

nof = zeros(n,1);
for i = 1:n
    for j = 1:length(kapal)
    if fender(i,1)-0.5 <kapal(j,1)+0.5 && fender(i,1)+0.5>kapal(j,1) && fender(i,2)>kapal(j,2) || kapal(j,2)<0
       nof(i) = 1;
    end
    end
end

deflection=zeros(n,1);
for i = 1:n
    if nof(i)==1
       deflection(i,1)= h-kapal(index1(i),2);
    end
end

rd = deflection/h*100;
energyratio = min(-0.0002826*(rd).^3+0.03703*(rd).^2+0.2029*rd-0.9977,100);

numberoffenderactivated = sum(nof);
% summary(o,1) = numberoffenderactivated;

for s = 1:n
        if energyratio(s,1)>0
            energyratio(s,1)=energyratio(s,1);
        else
            energyratio(s,1)=0;
        end
end
total=sum(energyratio)/100;
total1(o,1)=total;

figure
hold on
plot(xtrans,ytrans)
ylim([-10 200])
scatter(k,p)
plot(k,l)
hold off
end
toc

%% dengan sudut
LOA = 418;
Lpp = 395;
B = 56.4;
Fender = 0:14:500;
n = length(Fender);
h = 1.3; %tinggi fender
alpha = 0.0287;

z = LOA*0.25;
Rb= 0.5*(B/2+LOA^2/(8*B));
theta(:,1) = 0:0.5:(asind(z/Rb)+0.2);
theta(:,1) = sort(theta(:,1),'descend');

pusatx = 150; %centre of bow radius
pusaty = h;

for p = 1:length(alpha) %% berbagai macam alpha

for i = 1:length(theta); % bagian depan kapal
    x(i,1)=pusatx-Rb*sind(theta(i));
    y(i,1)=Rb-Rb*cosd(theta(i))+pusaty;
end

sisa= pusatx:1:(pusatx+LOA-z); % bagian lurusnya kapal

for i = 1:length(sisa)
    x(i+length(theta),1) = sisa(i);
    y(i+length(theta),1) = 0+pusaty;
end

L = sqrt(x.^2+y.^2);

xtrans = pusatx +cosd(alpha(p))*(x-pusatx)- sind(alpha(p))*(y-pusaty);
ytrans = pusaty +sind(alpha(p))*(x-pusatx)+ cosd(alpha(p))*(y-pusaty);

ytransmin = min(ytrans);
ytrans = ytrans-(ytransmin-(1-0.72)*h);

%fender coordinate
for i = 1:n
    k(i) = Fender(i);
    m(i) = h;
    l(i) = h*(1-0.72);
end

% coba itung fender activated
kapal(:,1) = xtrans;
kapal(:,2) = ytrans;
fender(:,1) = k;
fender(:,2) = m;

for i = 1:n
    for j = 1:length(kapal)
    selisihx(j,i)= abs(fender(i,1)-kapal(j,1));
    minimal = min(selisihx);
    index = find(abs(selisihx(:,i)) == minimal(1,i));
    end
    index1(i,1)=index;
end

nof = zeros(n,1);

for i = 1:n
    for j = 1:length(kapal)
    if fender(i,1)-0.5 <kapal(j,1)+0.5 && fender(i,1)+0.5>kapal(j,1) && fender(i,2)>kapal(j,2) || kapal(j,2)<0
       nof(i) = 1;
    end
    end
end

deflection=zeros(n,1)
for i = 1:n
    if nof(i)==1
       deflection(i,1)= h-kapal(index1(i),2);
    end
end

deflectionall(:,p) = deflection;
rd = deflectionall/h*100;
energyratio = min(-0.0002826*(rd).^3+0.03703*(rd).^2+0.2029*rd-0.9977,100);

numberoffenderactivated = sum(nof);
summary(p) = numberoffenderactivated;

figure
hold on
plot(xtrans,ytrans)
ylim([-10 200])
scatter(k,m)
plot(k,l)
hold off
end

for s = 1:n
    for t = 1:length(alpha)
        if energyratio(s,t)>0
            energyratio(s,t)=energyratio(s,t);
        else
            energyratio(s,t)=0;
        end
    end
end

figure
scatter(alpha,summary);
xlabel('alpha')
ylabel('number of fender')

sum=sum(energyratio)/100;
figure
scatter(alpha,sum);
xlabel('berthing angle')
ylabel('energy ratio of fender')
%% dengan sudut 1 aja
clear all;
close all

LOA = 400;
Lpp = 385;
B = 59;
Fender = 0:14:500;
n = length(Fender);
h = 1.4; %tinggi fender
alpha = 1.5;

z = LOA*0.25;
Rb= 0.5*(B/2+LOA^2/(8*B));
theta(:,1) = 0:0.1:(asind(z/Rb)+0.2);
theta(:,1) = sort(theta(:,1),'descend')

pusatx = 150; %centre of bow radius
pusaty = 0.72*h;

for i = 1:length(theta);
    x(i,1)=pusatx-Rb*sind(theta(i));
    y(i,1)=Rb-Rb*cosd(theta(i))+0.72*h;
end

sisa= pusatx:1:(pusatx+LOA-z);

for i = 1:length(sisa)
    x(i+length(theta),1) = sisa(i);
    y(i+length(theta),1) = 0+0.72*h;
end

L = sqrt(x.^2+y.^2);

xtrans = pusatx +cosd(alpha)*(x-pusatx)- sind(alpha)*(y-pusaty);
ytrans = pusaty +sind(alpha)*(x-pusatx)+ cosd(alpha)*(y-pusaty);

%fender coordinate
for i = 1:n
    k(i) = Fender(i);
    m(i) = h;
end

xmirror = x-(2*(y(1)-y)*cosd(90-alpha));
ymirror = y+(2*(y(1)-y)*sind(90-alpha));

% coba itung fender activated
kapal(:,1) = xtrans;
kapal(:,2) = ytrans;
fender(:,1) = k;
fender(:,2) = m;

nof = zeros(n,1);
for i = 1:n;
    for j = 1:length(kapal);
    if fender(i,1)-0.5 <kapal(j,1)+0.5 && fender(i,1)+0.5>kapal(j,1) && fender(i,2)>kapal(j,2) || kapal(j,2)<0;;
        nof(i) = 1;
    end
    end
end

numberoffenderactivated = sum(nof)

figure
hold on
plot(xtrans,ytrans)
ylim([-10 200])
plot(k,m)
scatter(k,m)
hold off

%% the last dapet daritest gambar
close all;
clear total1;
clear sudut;
clear nof1;
m=1;

for o = 1:m
LOA1 = LOA(o);
Lpp1 = LPP(o);
B1   = B(o);
Fender = 0:14:500;
n = length(Fender);
h = (1400)/1000; %tinggi fender
alpha = angle(o);

z = LOA1*0.25;
Rb= 0.5*(B1/2+LOA1^2/(8*B1));
theta(o)= {0:0.5:(asind(z/Rb)+0.2)};
theta1 = cell2mat(theta(o));
theta1 = sort(theta1,'descend');

pusatx = 140; %centre of bow radius
pusaty = h;

for i = 1:length(theta1) % bagian depan kapal
     x(i,1)=pusatx-Rb*sind(theta1(i));
     y(i,1)=Rb - Rb*cosd(theta1(i))+ pusaty;
end

sisa= pusatx:1:(pusatx+LOA-z); % bagian lurusnya kapal

for i = 1:length(sisa)
    x(i+length(theta1),1) = sisa(i);
    y(i+length(theta1),1) = 0+pusaty;
end

L = sqrt(x.^2+y.^2);

xtrans1 = pusatx +cosd(alpha)*(x-pusatx)- sind(alpha)*(y-pusaty);
ytrans1 = pusaty +sind(alpha)*(x-pusatx)+ cosd(alpha)*(y-pusaty);

xtrans = pusatx +cosd(alpha)*(x-pusatx)- sind(alpha)*(y-pusaty);
ytrans = pusaty +sind(alpha)*(x-pusatx)+ cosd(alpha)*(y-pusaty);
ytransmin =min(ytrans);

%fender coordinate
for i = 1:n
    k(i) = Fender(i);
    m(i) = h;
    l(i) = h*(1-0.72);
end

fender(:,1) = k;
fender(:,2) = m;

nof = zeros(n,1);
for i = 1:n
    for j = 1:length(xtrans)
    if fender(i,1)-0.2 < xtrans(j,1)+0.2 && fender(i,1)+0.2 > xtrans(j,1) && fender(i,2)>=ytrans(j,1)
       nof(i) = 1;
    end
    end
end

sumnof=sum(nof);
if sumnof == 1;
   ytrans = ytrans-0.72*h;
else
   ytrans = ytrans-(ytransmin-(1-0.72)*h);
end
%

if any(ytrans(:)<0)
    ytrans = ytrans-min(ytrans);
end
    
% coba itung fender activated
kapal = zeros(length(xtrans),2);
kapal(:,1) = xtrans;
kapal(:,2) = ytrans;


for i = 1:n
    for j = 1:length(kapal)
    selisihx(j,i)= abs(fender(i,1)-kapal(j,1));
    minimal = min(selisihx);
    index = min(find(abs(selisihx(:,i)) == minimal(1,i), 1 ));
    end
    index1(i,1)=index;
end

nof = zeros(n,1);
for i = 1:n
    for j = 1:length(kapal)
    if fender(i,1)-0.5 <kapal(j,1)+0.5 && fender(i,1)+0.5>kapal(j,1) && fender(i,2)>kapal(j,2) || kapal(j,2)<0
       nof(i) = 1;
    end
    end
end

deflection=zeros(n,1);
for i = 1:n
    if nof(i)==1
       deflection(i,1)= h-kapal(index1(i),2);
    end
end

rd = deflection/h*100;
energyratio = min(-0.0002826*(rd).^3+0.03703*(rd).^2+0.2029*rd-0.9977,100);

numberoffenderactivated = sum(nof);
summary(o,1) = numberoffenderactivated;

for s = 1:n
        if energyratio(s,1)>0
            energyratio(s,1)=energyratio(s,1);
        else
            energyratio(s,1)=0;
        end
end
total=sum(energyratio)/100;
total1(o,1)=total;
Rb1(o,1)=Rb;

figure
hold on
plot(xtrans,ytrans)
xlim([0 350])
ylim([-10 100])
scatter(k,m)
plot(k,l)
hold off
end


