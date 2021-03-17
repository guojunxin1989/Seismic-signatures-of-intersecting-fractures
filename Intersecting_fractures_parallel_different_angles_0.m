clear;
close;
clc;

Km=26;%bulk modulus of the background medium,GPa
Gm=31;%shear modulus of the background medium
por=0.1;%porosity of the background medium
Ks=37;%bulk modulus of the solid 
Kf=2.25;%bulk modulus of the fluid
dens=2.65;%density of the solid
denf=1;%density of the fluid
vis=0.001*10^(-9);%viscosity of fluid, Gpa*s
perm=10^(-15);%permeability of background, m^2
d=1;%diameter of the fracture,m
r=d/2;%radius of the fracture,m
h=0.01;%thickness
c=h/2;%fracture half-thickness
cona=4*10^(-10);%dimensionless
con=cona*c/vis;

ep=0.000001;
angN=[ep pi/6 pi/4 pi/3 pi/2-ep];
AN=length(angN);
Iang=ep;
f=-3:0.1:5;
f=10.^f;
Df=length(f);
v=zeros(Df,AN);
inQ=zeros(Df,AN);

for AI=1:AN
    
ang1=angN(AI);
ang2=pi-ang1-Iang;

cden1=0.025;%fracture density of first set fracture
cden2=0.025;%fracture density of second set fracture

L=Km+4/3*Gm;
a=1-Km/Ks;
M=((a-por)/Ks+por/Kf)^(-1);
C=a*M;
H=L+a^2*M;
den=dens*(1-por)+denf*por;

parfor FN=1:Df
    
FN

ft=f(FN);
w=2*pi*ft;
 
abk=(2*pi)/d*5;
DM=5000;
k=abk/DM:abk/DM:abk;

wb=por*vis/(denf*perm*10^(-6));
permb=perm*(sqrt(1-1i*w/(2*wb))-1i*w/wb)^(-1);
denb=(1i*vis)/(w*permb)*10^6;
b=(den*M+denb*H-2*denf*C)/(2*(M*H-C^2));
%s1=sqrt(b-sqrt(b^2-(den*denb-denf^2)/(M*H-C^2)));
s1=sqrt((den*denb-denf^2)/(M*H-C^2)/(b+sqrt(b^2-(den*denb-denf^2)/(M*H-C^2))));
s2=sqrt(b+sqrt(b^2-(den*denb-denf^2)/(M*H-C^2)));
s3=sqrt((den*denb-denf^2)/(Gm*denb));

k1=w*s1*10^(-3);
k2=w*s2*10^(-3);
k3=w*s3*10^(-3);

x1=-(H*s1^2-den)/(C*s1^2-denf);
x2=-(H*s2^2-den)/(C*s2^2-denf);
x3=-denf/denb;

kN=length(k);

et1=zeros(1,kN);
et2=zeros(1,kN);
et3=zeros(1,kN);

for lk=1:kN
    
kk=k(lk);
if kk<=k1
    et1(lk)=-1i*sqrt(k1^2-kk^2);
else
    et1(lk)=sqrt(kk^2-k1^2);
end
    et2(lk)=-1i*sqrt(k2^2-kk^2);
if kk<=k3
    et3(lk)=-1i*sqrt(k3^2-kk^2);
else
    et3(lk)=sqrt(kk^2-k3^2);
end

end


% normal fracture discontinuity 

D1=(H-C+C*x1-M*x1)/Gm;
D2=(H-C+C*x2-M*x2)/Gm;
E1=2*(x2-x3)*Gm*k1^2*(D1-1);
E2=2*(x1-x3)*Gm*k2^2*(D2-1);
E=E1-E2;

a11=-(1+x1)*et1+(1+x3)*(2*et1.*(k.^2))./(2*k.^2-k3^2)+(1i*con/w+c/Kf)*k1^2*(C+M*x1);
a12=-1i*con/w*k1^2*(C+M*x1)*ones(1,kN);
b11=(1+x2)*et2-(1+x3)*(2*et2.*(k.^2))./(2*k.^2-k3^2)-(1i*con/w+c/Kf)*k2^2*(C+M*x2);
b12=1i*con/w*k2^2*(C+M*x2)*ones(1,kN);

c11=zeros(1,kN);
c12=zeros(1,kN);

for lk=1:kN

AM=[a11(lk),a12(lk);a12(lk),a11(lk)];
BM=[b11(lk),b12(lk);b12(lk),b11(lk)];
CM=BM\AM;
c11(lk)=CM(1,1);
c12(lk)=CM(1,2);
end

d11=(2*Gm*k.^2-Gm*k1^2*D1-(4*Gm*((k.^2).*et3).*et1)./(2*k.^2-k3^2))+(2*Gm*k.^2-Gm*k2^2*D2-(4*Gm*((k.^2).*et3).*et2)./(2*k.^2-k3^2)).*c11;
d12=(2*Gm*k.^2-Gm*k2^2*D2-(4*Gm*((k.^2).*et3).*et2)./(2*k.^2-k3^2)).*c12;

Hf1=((2*k.^2*(x3-x2)+k3^2*(1+x2))./(et1.*k*E)).*(d11+d12)-1;
Hf2=((2*k.^2*(x3-x2)+k3^2*(1+x2))./(et1.*k*E)).*(d11-d12)-1;

N=50;
ds=r/N;
Ma1=zeros(N,N);
Ma2=zeros(N,N);
V1=zeros(N,1);
V2=zeros(N,1);

Mx=0;
sumM1=0;
sumM2=0;

for mf=0:Mx

if(mf==0)
    Em=1;
else
    Em=2;
end

for I=1:N
    for J=1:N
       si=(I-1)*ds;
       sj=(J-1)*ds; 
       
       HK1=((k.*Hf1).*besselj(mf+0.5,si*k)).*besselj(mf+0.5,sj*k);
       My1=trapz(k,HK1)*(si*sj)^0.5;
       HK2=((k.*Hf2).*besselj(mf+0.5,si*k)).*besselj(mf+0.5,sj*k);
       My2=trapz(k,HK2)*(si*sj)^0.5;
       if(I==J)
           Ma1(I,J)=My1*ds+1;
           Ma2(I,J)=My2*ds+1;
       else
           Ma1(I,J)=My1*ds;
           Ma2(I,J)=My2*ds;
       end
       
    end
    
    
    if(si==0)
      if(mf>0)
          Tmv1=0;
          Tmv2=0;
      else
          Tmv1=sqrt(2/pi);
          Tmv2=sqrt(2/pi);
      end
    else
       Tmv1=besselj(mf+0.5,k1*si*sin(ang1))/sqrt(k1*si*sin(ang1)); 
       Tmv2=besselj(mf+0.5,k1*si*sin(ang2))/sqrt(k1*si*sin(ang2)); 
    end
    
    p01=Gm*k1*(D1-2*sin(ang1)^2)*Em*(1i)^(mf+1)*sqrt(pi/2)*Tmv1;
    p02=Gm*k1*(D1-2*sin(ang2)^2)*Em*(1i)^(mf+1)*sqrt(pi/2)*Tmv2;
    p1=p01+p02;
    p2=p01-p02;
    V1(I)=-p1*si;
    V2(I)=-p2*si;
end

P1=Ma1\V1;
P2=Ma2\V2;

sum1=0;
sum2=0;
sum3=0;
sum4=0;
for I=1:N
    si=(I-1)*ds;
    sum1=sum1+si^0.5*P1(I)*besselj(mf+0.5,si*k1*sin(ang1))*ds;
    sum2=sum2+si^0.5*P2(I)*besselj(mf+0.5,si*k1*sin(ang1))*ds;
    sum3=sum3+si^0.5*P1(I)*besselj(mf+0.5,si*k1*sin(ang2))*ds;
    sum4=sum4+si^0.5*P2(I)*besselj(mf+0.5,si*k1*sin(ang2))*ds;
end

T1=sqrt(2*k1*sin(ang1)/pi)*sum1;
T2=sqrt(2*k1*sin(ang1)/pi)*sum2;
A1s=T1*(2*(k1*sin(ang1))^2*(x3-x2)+k3^2*(1+x2))/(-1i*k1*cos(ang1)*k1*sin(ang1)*E);
A2s=T2*(2*(k1*sin(ang1))^2*(x3-x2)+k3^2*(1+x2))/(-1i*k1*cos(ang1)*k1*sin(ang1)*E);
A1=(A1s+A2s)/2;

T3=sqrt(2*k1*sin(ang2)/pi)*sum3;
T4=sqrt(2*k1*sin(ang2)/pi)*sum4;
A3s=T3*(2*(k1*sin(ang2))^2*(x3-x2)+k3^2*(1+x2))/(-1i*k1*cos(ang2)*k1*sin(ang2)*E);
A4s=T4*(2*(k1*sin(ang2))^2*(x3-x2)+k3^2*(1+x2))/(-1i*k1*cos(ang2)*k1*sin(ang2)*E);
A2=(A3s-A4s)/2;

sumM1=sumM1+(-1i)^mf*cos(ang1)*A1;
sumM2=sumM2+(-1i)^mf*cos(ang2)*A2;

end

f0=k1^2*sumM1;
f1=k1^2*sumM2;

% shear fracture discontinuity

E=(C^2-(H-Gm)*M)*(x1-x2)*k1^2/(Gm*(C+M*x2));

F1=2*et1-(2*k.^2-k3^2)./et3+(k1^2*(2*k.^2-k3^2)*(H+C*x1))./((2*Gm*et3).*(k.^2));
F2=(k1^2*(C+M*x1))/(k2^2*(C+M*x2))*(2*et2-(2*k.^2-k3^2)./et3+(k2^2*(2*k.^2-k3^2)*(H+C*x2))./((2*Gm*et3).*(k.^2)));
Hf=(F1-F2).*k/E-1;

N=50;
ds=r/N;
Ma=zeros(N,N);
V1=zeros(N,1);
V2=zeros(N,1);

Mx=0;
sumM1=0;
sumM2=0;

for mf=0:Mx

if(mf==0)
    Em=1;
else
    Em=2;
end

for I=1:N
    for J=1:N
       si=(I-1)*ds;
       sj=(J-1)*ds; 
       
       HK=((k.*Hf).*besselj(mf+0.5,si*k)).*besselj(mf+0.5,sj*k);
       My=trapz(k,HK)*(si*sj)^0.5;
       if(I==J)
           Ma(I,J)=My*ds+1;
       else
           Ma(I,J)=My*ds;
       end
       
    end
    
    if(si==0)
      if(mf>0)
          Tmv1=0;
      else
          Tmv1=sqrt(2/pi);
      end
    else
       Tmv1=besselj(mf+0.5,k1*si*sin(ang1))/sqrt(k1*si*sin(ang1)); 
    end
    
    p1=2*cos(ang1)*Em*(1i)^mf*sqrt(pi/2)*Tmv1;
    V1(I)=p1*si;
    
    if(si==0)
      if(mf>0)
          Tmv2=0;
      else
          Tmv2=sqrt(2/pi);
      end
    else
       Tmv2=besselj(mf+0.5,k1*si*sin(ang2))/sqrt(k1*si*sin(ang2)); 
    end
    
    p2=2*cos(ang2)*Em*(1i)^mf*sqrt(pi/2)*Tmv2;
    V2(I)=p2*si;
end

P1=Ma\V1;
P2=Ma\V2;
sum1=0;
sum2=0;
for I=1:N
    si=(I-1)*ds;
    sum1=sum1+si^0.5*P1(I)*besselj(mf+0.5,si*k1*sin(ang1))*ds;
    sum2=sum2+si^0.5*P2(I)*besselj(mf+0.5,si*k1*sin(ang2))*ds;
end

sum1=sum1*sqrt(2*k1*sin(ang1)/pi);
sum2=sum2*sqrt(2*k1*sin(ang2)/pi);
sumM1=sumM1+sum1*(-1i)^mf;
sumM2=sumM2+sum2*(-1i)^mf;

end

f2=(k1^3*sin(2*ang1))/(2*E)*sumM1;
f3=(k1^3*sin(2*ang2))/(2*E)*sumM2;

fT1=f0+f2;
fT2=f1+f3;

n01=cden1/r^3;
n02=cden2/r^3;
keff=k1*(1+(2*pi*n01)/k1^2*fT1+(2*pi*n02)/k1^2*fT2);
v(FN,AI)=(2*pi*ft)/real(keff);
inQ(FN,AI)=2*imag(keff)/real(keff);

end

end

figure(1)
semilogx(f,v(:,1),'m','linewidth',1.3);
hold on;
semilogx(f,v(:,2),'r','linewidth',1.3);
hold on;
semilogx(f,v(:,3),'g','linewidth',1.3);
hold on;
semilogx(f,v(:,4),'b','linewidth',1.3);
hold on;
semilogx(f,v(:,5),'k','linewidth',1.3);

xlabel('Frequency(Hz)','FontSize',12);
ylabel('{\itv_{p}}(m/s)','FontSize',12);
legend('\theta = 0','\theta = 30','\theta = 45','\theta = 60','\theta = 90');
xlim([0.001 100000]);

figure(2)
loglog(f,inQ(:,1),'m','linewidth',1.3);
hold on;
loglog(f,inQ(:,2),'r','linewidth',1.3);
hold on;
loglog(f,inQ(:,3),'g','linewidth',1.3);
hold on;
loglog(f,inQ(:,4),'b','linewidth',1.3);
hold on;
loglog(f,inQ(:,5),'k','linewidth',1.3);

xlabel('Frequency(Hz)','FontSize',12);
ylabel('1/{\itQ_{p}}','FontSize',12);
legend('\theta = 0','\theta = 30','\theta = 45','\theta = 60','\theta = 90');
xlim([0.001 100000]);

save parallel_case f v inQ;


