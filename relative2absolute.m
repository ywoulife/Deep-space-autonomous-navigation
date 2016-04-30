function [r2,v2]=relative2absolute(r_v0,r_v1)
%已知参考卫星的轨道根数，那就已经知道三个旋转角度，即旋转矩阵
%输入为参考卫星的绝对信息和两颗卫星之间的相对信息，首先要将相对量转为绝对坐标系下的量
%然后再通过矢量加减求得伴随卫星的绝对位置速度
r0=r_v0(1:3);v0=r_v0(4:6);
r1=r_v1(1:3);v1=r_v1(4:6);
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
r2=B1\r1+r0;%将相对量转到惯性系中
v2=B1\(v1+(cross([0,0,hm/rm^2],r1'))')+v0;
  %伴随卫星绝对位置位置和速度，即x1+x0
