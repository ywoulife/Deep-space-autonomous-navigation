function R=absolute2relative(r_va,r_vb)
%由位置和速度矢量求轨道根数
r0=r_va(1:3);v0=r_va(4:6);
r1=r_vb(1:3);v1=r_vb(4:6);
h=cross(r0,v0);     %位置和速度的叉乘
hm=norm(h); %三个坐标的距离
rm=norm(r0);
i=acos(h(3)/hm);
b=atan2(h(1),-h(2)); 
u=atan2(r0(3),(r0(2)*sin(b)+r0(1)*cos(b))*sin(i));
B=[cos(b)*cos(u)-sin(b)*sin(u)*cos(i),-cos(b)*sin(u)-sin(b)*cos(u)*cos(i),...
  sin(b)*sin(i);sin(b)*cos(u)+cos(b)*sin(u)*cos(i),-sin(b)*sin(u)+cos(b)*cos(u)*cos(i),...
   -cos(b)*sin(i);sin(u)*sin(i),cos(u)*sin(i),cos(i)];
B1=B';
dr=r1-r0;dv=v1-v0;   %相对位置和速度，即x1-x0
R(1:3)=B1*dr;
R(4:6)=B1*dv-(cross([0,0,hm/rm^2],R(1:3)))';
