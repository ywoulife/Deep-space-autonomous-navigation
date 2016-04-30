%考虑地球非球形引力J2项，大气阻力摄动力
%进行轨道的状态估计
%趁这个机会把之前本科的滤波的书籍都复习一遍，进一步了解更新的滤波算法
%首先使用ＥＫＦ进行
clear all
clc
Virgin;         %初始化本次相关的参数
% InitSimFun;     %初始化全局参数
a=6878.136;
e=0.0003;
incl=0.001;          %inclination
%----------------根据给定的设计参数设计构形----------------------------------
[r3]=formationdesign(a,e,RA,incl,omega,f,R1);%R 为绕飞半径
xt(:,1)=r3;
%--------------------龙格库塔轨道递推----------------------------------
OldSat=r3;

%---------------------------------------------------------------
for i=2:N
  pnoise=w_k(:,i);
  onoise=v_k(:,i);
  NewSat = RKFixed4(T,OldSat,@diff_eq)+pnoise;%,para1,para2,para3,...
  yt(:,i)=H*NewSat+onoise;
  xt(:,i)=NewSat;
  OldSat=NewSat;
end
%----------------------------extended kalman filter----------------------
xe=zeros(12,N);
xe(:,1)=xt(:,1)+sqrt(P0)*randn(12,1);
P=P0;
p11sp=zeros(1,N);p22sp=zeros(1,N);p33sp=zeros(1,N);p44sp=zeros(1,N);p55sp=zeros(1,N);
p66sp=zeros(1,N);p77sp=zeros(1,N);p88sp=zeros(1,N);p99sp=zeros(1,N);p1010sp=zeros(1,N);
p1111sp=zeros(1,N);p1212sp=zeros(1,N);
p11sn=zeros(1,N);p22sn=zeros(1,N);p33sn=zeros(1,N);p44sn=zeros(1,N);p55sn=zeros(1,N);
p66sn=zeros(1,N);p77sn=zeros(1,N);p88sn=zeros(1,N);p99sn=zeros(1,N);p1010sn=zeros(1,N);
p1111sn=zeros(1,N);p1212sn=zeros(1,N);
p11sp(1)=sqrt(P0(1,1))*1e3;
p11sn(1)=(-1)*sqrt(P0(1,1))*1e3;
p22sp(1)=sqrt(P0(2,2))*1e3;
p22sn(1)=(-1)*sqrt(P0(2,2))*1e3;
p33sp(1)=sqrt(P0(3,3))*1e3;
p33sn(1)=(-1)*sqrt(P0(3,3))*1e3;
p44sp(1)=sqrt(P0(4,4))*1e3;
p44sn(1)=(-1)*sqrt(P0(4,4))*1e3;
p55sp(1)=sqrt(P0(5,5))*1e3;
p55sn(1)=(-1)*sqrt(P0(5,5))*1e3;
p66sp(1)=sqrt(P0(6,6))*1e3;
p66sn(1)=(-1)*sqrt(P0(6,6))*1e3;
p77sp(1)=sqrt(P0(7,7))*1e3;
p77sn(1)=(-1)*sqrt(P0(7,7))*1e3;
p88sp(1)=sqrt(P0(8,8))*1e3;
p88sn(1)=(-1)*sqrt(P0(8,8))*1e3;
p99sp(1)=sqrt(P0(9,9))*1e3;
p99sn(1)=(-1)*sqrt(P0(9,9))*1e3;
p1010sp(1)=sqrt(P0(10,10))*1e3;
p1010sn(1)=(-1)*sqrt(P0(10,10))*1e3;
p1111sp(1)=sqrt(P0(11,11))*1e3;
p1111sn(1)=(-1)*sqrt(P0(11,11))*1e3;
p1212sp(1)=sqrt(P0(12,12))*1e3;
p1212sn(1)=(-1)*sqrt(P0(12,12))*1e3;
h=T;
for k=2:N
   %state  transmition matrix
    x_k=xe(:,k-1);%x_k is only the intermediate value
    r_a=sqrt(x_k(1)^2+x_k(2)^2+x_k(3)^2);%|0，I |
    r_b=sqrt(x_k(7)^2+x_k(8)^2+x_k(9)^2);%|FA，0|
    F_A=[-u*(r_a^2-3*x_k(1)^2)/r_a^5,u*(3*x_k(1)*x_k(2))/r_a^5,u*(3*x_k(1)*x_k(3)/r_a^5);...
    u*(3*x_k(1)*x_k(2))/r_a^5,-u*(r_a^2-3*x_k(2)^2)/r_a^5,u*3*x_k(2)*x_k(3)/r_a^5;...
    u*(3*x_k(1)*x_k(3))/r_a^5,u*3*x_k(2)*x_k(3)/r_a^5,-u*(r_a^2-3*x_k(3)^2)/r_a^5];
    Fa=[zeros(3),eye(3);F_A,zeros(3)];%|0，I  |
                                      %|F_A，0|
    F_B=[-u*(r_b^2-3*x_k(7)^2)/r_b^5,u*(3*x_k(7)*x_k(8))/r_b^5,u*(3*x_k(7)*x_k(9)/r_b^5);...
    u*(3*x_k(7)*x_k(8))/r_b^5,-u*(r_b^2-3*x_k(8)^2)/r_b^5,u*3*x_k(8)*x_k(9)/r_b^5;...
    u*(3*x_k(7)*x_k(9))/r_b^5,u*3*x_k(8)*x_k(9)/r_b^5,-u*(r_b^2-3*x_k(9)^2)/r_b^5];
    Fb=[zeros(3),eye(3);F_B,zeros(3)];
    A=[Fa,zeros(6);zeros(6),Fb];%12x12状态转移矩阵|Fa,0| 雅可比矩阵
    %状态转移矩阵
    phi=eye(12)+A*T+1/2*A*A*T^2;
    %--------------------------------
    OldSat=x_k;
    x_m = RKFixed4(T,OldSat,@diff_eq);
%观测方程，计算出实测的观测值
    P_m=phi*P*phi'+Q;
    %observation 
    y_m=H*x_m;
    %measurement update
    K=P_m*H'/(H*P_m*H'+R);
    xe(:,k)=x_m+K*(yt(:,k)-y_m);
    P=(eye(12)-K*H)*P_m;
    x_err(:,k)=(xt(:,k)-xe(:,k))*1e3;
    
%     
%   perr=(xt(:,k)-xe(:,k))*1e3;
%   pekf(k)=sqrt(perr(1)^2+perr(2)^2+perr(3)^2);
    p11sp(k)=sqrt(P(1,1))*1e3;
    p11sn(k)=(-1)*sqrt(P(1,1))*1e3;
    p22sp(k)=sqrt(P(2,2))*1e3;
    p22sn(k)=(-1)*sqrt(P(2,2))*1e3;
    p33sp(k)=sqrt(P(3,3))*1e3;
    p33sn(k)=(-1)*sqrt(P(3,3))*1e3;
    p44sp(k)=sqrt(P(4,4))*1e3;
    p44sn(k)=(-1)*sqrt(P(4,4))*1e3;
    p55sp(k)=sqrt(P(5,5))*1e3;
    p55sn(k)=(-1)*sqrt(P(5,5))*1e3;
    p66sp(k)=sqrt(P(6,6))*1e3;
    p66sn(k)=(-1)*sqrt(P(6,6))*1e3;
    p77sp(k)=sqrt(P(7,7))*1e3;
    p77sn(k)=(-1)*sqrt(P(7,7))*1e3;
    p88sp(k)=sqrt(P(8,8))*1e3;
    p88sn(k)=-1*sqrt(P(8,8))*1e3;
    p99sp(k)=sqrt(P(9,9))*1e3;
    p99sn(k)=(-1)*sqrt(P(9,9))*1e3;
    p1010sp(k)=sqrt(P(10,10))*1e3;
    p1010sn(k)=(-1)*sqrt(P(10,10))*1e3;
    p1111sp(k)=sqrt(P(11,11))*1e3;
    p1111sn(k)=(-1)*sqrt(P(11,11))*1e3;
    p1212sp(k)=sqrt(P(12,12))*1e3;
    p1212sn(k)=(-1)*sqrt(P(12,12))*1e3;

end
set(0,'DefaultFigureWindowStyle','docked');
t=T*(1:N);
figure(1);%绘制图像，误差，标准差
subplot(3,1,1),plot(t,x_err(1,:),'b-',t,p11sp,'r-.',t,p11sn,'r-.');grid on;
xlabel('Time(sec)');ylabel('x Position Error(m)');title('Position Estimation Error of Satellite A');axis([0,10000,-3000,3000]);
subplot(3,1,2),plot(t,x_err(2,:),'b-',t,p22sp,'r-.',t,p22sn,'r-.');
xlabel('Time(sec)');ylabel('y Position Error(m)');axis([0,10000,-2000,2000]);grid on;
subplot(3,1,3),plot(t,x_err(3,:),'b-',t,p33sp,'r-.',t,p33sn,'r-.');
xlabel('Time(sec)');ylabel('z Position Error(m)');axis([0,10000,-2000,2000]);grid on;
figure(2);
subplot(3,1,1),plot(t,x_err(7,:),'b-',t,p77sp,'r-.',t,p77sn,'r-.');grid on;
xlabel('Time(sec)');ylabel('x Position Error(m)');title('Position Estimation Error of Satellite B');axis([0,10000,-3000,3000]);
subplot(3,1,2),plot(t,x_err(8,:),'b-',t,p88sp,'r-.',t,-1*p88sp,'r-.');
xlabel('Time(sec)');ylabel('y Position Error(m)');axis([0,10000,-2000,2000]);grid on;
subplot(3,1,3),plot(t,x_err(9,:),'b-',t,p99sp,'r-.',t,p99sn,'r-.');
xlabel('Time(sec)');ylabel('z Position Error(m)');axis([0,10000,-2000,2000]);grid on;
figure(3);
subplot(3,1,1),plot(t,x_err(4,:),'b-',t,p44sp,'r-.',t,p44sn,'r-.');grid on;
xlabel('Time(sec)');ylabel('x Velocity Error(m/s)');title('Velocity Estimation Error of Satellite A');axis([0,10000,-0.5,0.5]);
subplot(3,1,2),plot(t,x_err(5,:),'b-',t,p55sp,'r-.',t,p55sn,'r-.');
xlabel('Time(sec)');ylabel('y Velocity Error(m/s)');axis([0,10000,-0.5,0.5]);grid on;
subplot(3,1,3),plot(t,x_err(6,:),'b-',t,p66sp,'r-.',t,p66sn,'r-.');
xlabel('Time(sec)');ylabel('z Velocity Error(m/s)');axis([0,10000,-0.5,0.5]);grid on;
figure(4);
subplot(3,1,1),plot(t,x_err(10,:),'b-',t,p1010sp,'r-.',t,p1010sn,'r-.');grid on;
xlabel('Time(sec)');ylabel('x Velocity Error(m/s)');title('Position Estimation Error of Satellite B');axis([0,10000,-0.5,0.5]);
subplot(3,1,2),plot(t,x_err(11,:),'b-',t,p1111sp,'r-.',t,p1111sn,'r-.');
xlabel('Time(sec)');ylabel('y Velocity Error(m/s)');axis([0,10000,-0.5,0.5]);grid on;
subplot(3,1,3),plot(t,x_err(12,:),'b-',t,p1212sp,'r-.',t,p1212sn,'r-.');
xlabel('Time(sec)');ylabel('z Velocity Error(m/s)');axis([0,10000,-0.5,0.5]);
grid on;
figure(5);
plot3(xe(1,:),xe(2,:),xe(3,:),'k');hold on
plot3(xe(7,:),xe(8,:),xe(9,:),'k');
grid on;
