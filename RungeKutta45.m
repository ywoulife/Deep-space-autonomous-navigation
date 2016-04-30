function [xt1]=RungeKutta45(state,step,Drag12state,pnoise)
    h=step;
    k1=diff_eq(state,Drag12state,pnoise);
    k2=diff_eq(state+h*k1/2,Drag12state,pnoise);
    k3=diff_eq(state+h*k2/2,Drag12state,pnoise);
    k4=diff_eq(state+h*k3,Drag12state,pnoise);
    xt1=state+h*(k1+2*k2+2*k3+k4)/6;
end