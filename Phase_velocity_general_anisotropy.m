function [D,V]=Phase_velocity_general_anisotropy(N,C,den)
% This function is used to calculate the phase velocities for anisotropic
% media in any propagation directions.
% N: the wave propagation direction vector 
% C: stiffness matrix of the material
% den: density of material
% D: wave polorization direction vector 
% V: wave velocity 
M=zeros(3,3);
M(1,1)=C(1,1)*N(1)^2+2*C(1,6)*N(1)*N(2)+2*C(1,5)*N(1)*N(3)+2*C(5,6)*N(2)*N(3)+C(6,6)*N(2)^2+C(5,5)*N(3)^2;
M(1,2)=C(1,6)*N(1)^2+(C(1,2)+C(6,6))*N(1)*N(2)+(C(1,4)+C(5,6))*N(1)*N(3)+(C(4,6)+C(2,5))*N(2)*N(3)+C(2,6)*N(2)^2+C(4,5)*N(3)^2;
M(1,3)=C(1,5)*N(1)^2+(C(1,4)+C(5,6))*N(1)*N(2)+(C(1,3)+C(5,5))*N(1)*N(3)+(C(3,6)+C(4,5))*N(2)*N(3)+C(4,6)*N(2)^2+C(3,5)*N(3)^2;
M(2,2)=C(6,6)*N(1)^2+2*C(2,6)*N(1)*N(2)+2*C(4,6)*N(1)*N(3)+2*C(2,4)*N(2)*N(3)+C(2,2)*N(2)^2+C(4,4)*N(3)^2;
M(2,3)=C(5,6)*N(1)^2+(C(4,6)+C(2,5))*N(1)*N(2)+(C(3,6)+C(4,5))*N(1)*N(3)+(C(2,3)+C(4,4))*N(2)*N(3)+C(2,4)*N(2)^2+C(3,4)*N(3)^2;
M(3,3)=C(5,5)*N(1)^2+2*C(4,5)*N(1)*N(2)+2*C(3,5)*N(1)*N(3)+2*C(3,4)*N(2)*N(3)+C(4,4)*N(2)^2+C(3,3)*N(3)^2;
M(2,1)=M(1,2);
M(3,1)=M(1,3);
M(3,2)=M(2,3);

[D,V]=eig(M);
V=(V/den).^(1/2);

