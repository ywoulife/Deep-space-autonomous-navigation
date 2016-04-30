function [r,v]=coetostate(coe)
mu=398600.4418;
a=coe(1);
e=coe(2);
b=coe(3);
i=coe(4);
u=coe(5);
TA=coe(6);
p=a*(1-e^2);
r=p/(1+e*cos(TA));
h=sqrt(mu*p);
%-----------------------------calculate rp and vp 在近焦点坐标系
rp=r*[cos(TA);sin(TA);0];
vp=(mu/h)*[-sin(TA);e+cos(TA);0];
%----------------------------------
% R3_W=[cos(RA) sin(RA) 0;-sin(RA) cos(RA) 0;0 0 1];
% %-----------------------------------------
% R1_i=[1 0 0;0 cos(incl) sin(incl);0 -sin(incl) cos(incl)];
%   %--------------------------------------------
% R3_w=[cos(w) sin(w) 0;-sin(w) cos(w) 0;0 0 1];
B=[cos(b)*cos(u)-sin(b)*sin(u)*cos(i),-cos(b)*sin(u)-sin(b)*cos(u)*cos(i),...
  sin(b)*sin(i);sin(b)*cos(u)+cos(b)*sin(u)*cos(i),-sin(b)*sin(u)+cos(b)*cos(u)*cos(i),...
   -cos(b)*sin(i);sin(u)*sin(i),cos(u)*sin(i),cos(i)];%相对至绝对
%----------------------------------------

%-------------------------------------
r=B*rp;
v=B*vp;
%----------------------------------
r=r';
v=v';

end 