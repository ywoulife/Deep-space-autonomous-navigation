function  dx=diff_eq(OldSat)
u=398600.4418;
R=6378.14;
J=1.0826e-3;
r_a=sqrt(OldSat(1)^2+OldSat(2)^2+OldSat(3)^2);
r_b=sqrt(OldSat(7)^2+OldSat(8)^2+OldSat(9)^2);
%u=398600.4418;
dx(1,1)=OldSat(4);
dx(2,1)=OldSat(5);
dx(3,1)=OldSat(6);
dx(4,1)=-u*OldSat(1)/r_a^3*(1+1.5*J*(R/r_a)^2*(1-5*(OldSat(3)/r_a)^2));
dx(5,1)=-u*OldSat(2)/r_a^3*(1+1.5*J*(R/r_a)^2*(1-5*(OldSat(3)/r_a)^2));
dx(6,1)=-u*OldSat(3)/r_a^3*(1+1.5*J*(R/r_a)^2*(3-5*(OldSat(3)/r_a)^2));
dx(7,1)=OldSat(10);
dx(8,1)=OldSat(11);
dx(9,1)=OldSat(12);
dx(10,1)=-u*OldSat(7)/r_b^3*(1+1.5*J*(R/r_b)^2*(1-5*(OldSat(3)/r_b)^2));
dx(11,1)=-u*OldSat(8)/r_b^3*(1+1.5*J*(R/r_b)^2*(1-5*(OldSat(3)/r_b)^2));
dx(12,1)=-u*OldSat(9)/r_b^3*(1+1.5*J*(R/r_b)^2*(3-5*(OldSat(3)/r_b)^2));

end