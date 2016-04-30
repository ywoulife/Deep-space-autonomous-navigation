function f=diff_state(x)
% load w_k.mat;
% r = randi([1 1000],1,1);
% w=w_k(:,r);
R=6378.14;
J=1.0826e-3;
r_a=sqrt(x(1)^2+x(2)^2+x(3)^2);
r_b=sqrt(x(7)^2+x(8)^2+x(9)^2);
u=398600.4418;
f(1,1)=x(4);
f(2,1)=x(5);
f(3,1)=x(6);
f(4,1)=-u*x(1)/r_a^3*(1+1.5*J*(R/r_a)^2*(1-5*(x(3)/r_a)^2));
f(5,1)=-u*x(2)/r_a^3*(1+1.5*J*(R/r_a)^2*(1-5*(x(3)/r_a)^2));
f(6,1)=-u*x(3)/r_a^3*(1+1.5*J*(R/r_a)^2*(1-5*(x(3)/r_a)^2));
f(7,1)=x(10);
f(8,1)=x(11);
f(9,1)=x(12);
f(10,1)=-u*x(7)/r_b^3*(1+1.5*J*(R/r_b)^2*(1-5*(x(3)/r_b)^2));
f(11,1)=-u*x(8)/r_b^3*(1+1.5*J*(R/r_b)^2*(1-5*(x(3)/r_b)^2));
f(12,1)=-u*x(9)/r_b^3*(1+1.5*J*(R/r_b)^2*(1-5*(x(3)/r_b)^2));
end