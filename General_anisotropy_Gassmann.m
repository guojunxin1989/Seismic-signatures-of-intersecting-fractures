function C=General_anisotropy_Gassmann(C0,por,Ks,Kf)

C=zeros(6,6);
a=[1 1 1 0 0 0];
sum=0;
for i=1:3
    sum=sum+a(i)-(C0(1,i)+C0(2,i)+C0(3,i))/(3*Ks);
end
for i=1:6
    for j=1:6
      b1=a(i)-(C0(1,i)+C0(2,i)+C0(3,i))/(3*Ks);
      b2=a(j)-(C0(1,j)+C0(2,j)+C0(3,j))/(3*Ks);
      a1=1/(por*(1/Kf-1/Ks));
      D=1+a1/(3*Ks)*sum;
      C(i,j)=C0(i,j)+a1/D*b1*b2;
    end
end
