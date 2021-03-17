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
h=0.01;%thickness,m
c=h/2;%fracture half-thickness
cden1=0.025;%fracture density of first set fracture
cden2=0.025;%fracture density of second set fracture

ar=c/r;
p1=(4*pi*ar*cden1)/3;
p2=p1;

Em=(9*Km*Gm)/(3*Km+Gm);%Young's modulus of the background material
vm=Em/(2*Gm)-1;%Possion's ratio of the background material
den=(1-por-p1-p2)*dens+(por+p1+p2)*denf;%density of the sample

Kms=Km+(1-Km/Ks)^2/((1-por)/Ks+por/Kf-Km/Ks^2);
Ems=(9*Kms*Gm)/(3*Kms+Gm);
vms=Ems/(2*Gm)-1;
Las=Kms-2/3*Gm;%Lame constant for the saturated background material
Lms=Kms+4/3*Gm;%Compressional modulus for the saturated background material 

ang1A=0:pi/100:(2*pi);
NA=length(ang1A);
Iang=pi/2;

v1=zeros(NA,3);
N=[1,0,0];

for AN=1:NA

ang1=ang1A(AN);
ang2=pi-ang1-Iang;

Kmf=0.000001;
Gmf=0.000001;
a1=1;
a2=1;
a3=ar;
opt=3;

%low frequency

H1=General_Eshelby_model(Km,Gm,Kmf,Gmf,a1,a2,a3,p1,opt);
Znd1=H1(3,3);
Ztd1=H1(4,4);

N1=[-cos(ang1),sin(ang1),0];
Se1=Excess_compliance(Znd1,Ztd1,N1);
N2=[-cos(ang2),-sin(ang2),0];
Se2=Excess_compliance(Znd1,Ztd1,N2);

B1=1/Em;
B2=-vm/Em;
Sd=[B1, B2,B2,0,0,0;
    B2,B1,B2,0,0,0;
    B2,B2,B1,0,0,0;
    0,0,0,2*(B1-B2),0,0;
    0,0,0,0,2*(B1-B2),0;
    0,0,0,0,0,2*(B1-B2)];
S=Sd+Se1+Se2;
C0=inv(S);
C1=General_anisotropy_Gassmann(C0,por+p1+p2,Ks,Kf);
[D,V]=Phase_velocity_general_anisotropy(N,C1,den);
v1(AN,1)=V(3,3)*1000;

%intermediate frequency 

H2=General_Eshelby_model(Kms,Gm,Kmf,Gmf,a1,a2,a3,p1,opt);
Znd2=H2(3,3);
Ztd2=H2(4,4);

N1=[-cos(ang1),sin(ang1),0];
Se1=Excess_compliance(Znd2,Ztd2,N1);
N2=[-cos(ang2),-sin(ang2),0];
Se2=Excess_compliance(Znd2,Ztd2,N2);

B1=1/Ems;
B2=-vms/Ems;
Ss=[B1, B2,B2,0,0,0;
    B2,B1,B2,0,0,0;
    B2,B2,B1,0,0,0;
    0,0,0,2*(B1-B2),0,0;
    0,0,0,0,2*(B1-B2),0;
    0,0,0,0,0,2*(B1-B2)];
S=Ss+Se1+Se2;
C0=inv(S);
C2=General_anisotropy_Gassmann(C0,p1+p2,Kms,Kf);
[D,V]=Phase_velocity_general_anisotropy(N,C2,den);
v1(AN,2)=V(3,3)*1000;

%high frequency

Kmf=2.25;
Gmf=0.00001;

H3=General_Eshelby_model(Kms,Gm,Kmf,Gmf,a1,a2,a3,p1,opt);
Zns=H3(3,3);
Zts=H3(4,4);

N1=[-cos(ang1),sin(ang1),0];
Se1=Excess_compliance(Zns,Zts,N1);
N2=[-cos(ang2),-sin(ang2),0];
Se2=Excess_compliance(Zns,Zts,N2);

S=Ss+Se1+Se2;
C3=inv(S);
[D,V]=Phase_velocity_general_anisotropy(N,C3,den);
v1(AN,3)=V(3,3)*1000;

end

ang=ang1A*180/pi;
plot(ang,v1(:,1));
hold on;
plot(ang,v1(:,2));
hold on;
plot(ang,v1(:,3));

save fra_intersect_pi_2 ang v1










