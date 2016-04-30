function NewSat = RKFixed4(dt,OldSat,DfunHandle,varargin)
%功能：利用Runge_Kutta四级四阶定步长方法求解Dx = f(x)形式的微分方程；
%          适用于求解非刚性微分方程的定步长法，精度为O5(dt).
%Inputs: dt        ---积分步长
%        OldSat    ---上一步的状态,n*1的列向量
%        DfunHandle---计算导数的右函数句柄
%        varargin  ---变长度输入宗量，用于传递计算右导数需要的参数，1*n的元胞数组
%Output: NewSat    ---新积分步的状态，n*1的列向量
%Note 1：函数调用方法为：
%        NewSat = RungeKuttaFix4(dt,OldSat,@Dfun,para1,para2,para3,...)
%        其中,Dfun为计算右导数的函数名，para1,para2,para3,...为计算右导数时需要的参数
%Note2： Dfun的函数原型为： DY = Dfun(OldSat,varargin);
%Note3： 参数传递后的调用方法为：         
%         para1,para2,para3,...等传到Dfun函数中后，分别成为varargin中的第i个元胞，因此用如下形式本地化
%         para1 = varargin{1};
%         para2 = varargin{2};
%         para3 = varargin{3};
%         .......
%         之后就可以在Dfun中使用para1,para2,para3,...等参数
%Note4：  计算表明，对于定步长单步算法，四阶Runge-Kutta算法的精度已基本达到最高精度，六阶的精度提高不多
%
%First Edited by Zhang Hongbo on 2008-03-18
%Modified by Zhang Hongbo on 2011-03-28: 修改了变长度输入宗量的传递方法，可直接调用本函数


%1. 数据准备
C(1) = dt/2;    C(2) = dt/2;    C(3) = dt;  C(4) = dt;  C(5) = dt/2;
[Row,Col] = size(OldSat);
for cnt = 1:Row
    Tmp2Sat(cnt,1) = OldSat(cnt);
    NewSat (cnt,1) = OldSat(cnt);
end

%2. 四阶积分
for Step = 1:4
    DiffSat = feval(DfunHandle,Tmp2Sat,varargin{:});
    for cnt = 1:Row
        Tmp2Sat(cnt,1) = OldSat(cnt)+C(Step  )*DiffSat(cnt);
        NewSat (cnt,1) = NewSat(cnt)+C(Step+1)*DiffSat(cnt)/3;
    end
end
return;
        
        
        
        
        
        
        
        
        
        
        
