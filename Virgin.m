%定义全局变量 初始化
global u  RA  omega f R1 T J Re
global N H  Q R x P0 w_k v_k x_err xt yt
% global A2m Cd JED Mat

u=398600.4418;
RA=0;             %Right ascending
omega=0;          %近地点幅角
f=0;              %True Anomaly
R1=10;%单位km
T=5;
N=4000;%仿真次数
H=[eye(3),zeros(3),-eye(3),zeros(3)];%观测矩阵
%------------------------------------------------
%step2:define noise assumption
Q=diag([10^(-12) 10^(-12) 10^(-12) 10^(-16) 10^(-16) 10^(-16)...
    10^(-12) 10^(-12) 10^(-12) 10^(-16) 10^(-16) 10^(-16)]);%过程噪声方差
R=diag([10^(-8) 10^(-8) 10^(-8)]);%观测噪声方差
%-------------------------------------------------------------------------
%step3:initiate state and covariance

P0=diag([100 100 100 10^(-6) 10^(-6) 10^(-6)...
    100 100 100 100^(-6) 10^(-6) 10^(-6)]);%初始方差
x=zeros(12,N);

%------------------------------------------------------------------------

%simulation only calculate the true state trajectory for comparison
%also calculate the measurement vector
w_k=[10^(-6)*wgn(3,N,0,'real');10^(-8)*wgn(3,N,0,'real');10^(-6)*wgn(3,N,0,'real');10^(-8)*wgn(3,N,0,'real')];
v_k=10^(-4)*wgn(3,N,0,'real');%观测噪声
%-------------------------------------------
x_err=zeros(12,N);
xt=zeros(12,N);%位置速度状态
yt=zeros(3,N); %观测状态
%--------------计算状态更新，动力学模型包含大气摄动-------------
Re=6378.14;
J=1.0826e-3;
