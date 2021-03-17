function Se=Excess_compliance(Zn,Zt,N)
% calculate the excess compliances of the fractures at any orientation
% angles. 
% Zn: normal compliances of the fractures;
% Zt: shear compliance of the fractures;
% N: fracture normal vector
Se=zeros(6,6);
n1=N(1);
n2=N(2);
n3=N(3);
Se(1,1)=Zt*n1^2+(Zn-Zt)*n1^4;
Se(1,2)=(Zn-Zt)*n1^2*n2^2;
Se(1,3)=(Zn-Zt)*n1^2*n3^2;
Se(1,4)=2*(Zn-Zt)*n1^2*n2*n3;
Se(1,5)=Zt*n1*n3+2*(Zn-Zt)*n1^3*n3;
Se(1,6)=Zt*n1*n2+2*(Zn-Zt)*n1^3*n2;
Se(2,2)=Zt*n2^2+(Zn-Zt)*n2^4;
Se(2,3)=(Zn-Zt)*n2^2*n3^2;
Se(2,4)=Zt*n2*n3+2*(Zn-Zt)*n2^3*n3;
Se(2,5)=2*(Zn-Zt)*n2^2*n1*n3;
Se(2,6)=Zt*n1*n2+2*(Zn-Zt)*n2^3*n1;
Se(3,3)=Zt*n3^2+(Zn-Zt)*n3^4;
Se(3,4)=Zt*n2*n3+2*(Zn-Zt)*n3^3*n2;
Se(3,5)=Zt*n1*n3+2*(Zn-Zt)*n3^3*n1;
Se(3,6)=2*(Zn-Zt)*n3^2*n1*n2;
Se(4,4)=Zt*(n2^2+n3^2)+4*(Zn-Zt)*n2^2*n3^2;
Se(4,5)=Zt*n1*n2+4*(Zn-Zt)*n1*n2*n3^2;
Se(4,6)=Zt*n1*n3+4*(Zn-Zt)*n1*n2^2*n3;
Se(5,5)=Zt*(n1^2+n3^2)+4*(Zn-Zt)*n1^2*n3^2;
Se(5,6)=Zt*n2*n3+4*(Zn-Zt)*n1^2*n2*n3;
Se(6,6)=Zt*(n1^2+n2^2)+4*(Zn-Zt)*n1^2*n2^2;
Se(2,1)=Se(1,2);
Se(3,1)=Se(1,3);
Se(3,2)=Se(2,3);
Se(4,1)=Se(1,4);
Se(4,2)=Se(2,4);
Se(4,3)=Se(3,4);
Se(5,1)=Se(1,5);
Se(5,2)=Se(2,5);
Se(5,3)=Se(3,5);
Se(5,4)=Se(4,5);
Se(6,1)=Se(1,6);
Se(6,2)=Se(2,6);
Se(6,3)=Se(3,6);
Se(6,4)=Se(4,6);
Se(6,5)=Se(5,6);