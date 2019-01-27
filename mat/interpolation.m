% Symbol Computation For Interpolation Coeffiences
syms h;
x=-[-h/2 h/2 3/2*h 5/2*h];
A=[x.^0;x.^1;x.^2;x.^3];
b=[0;-1;0;0];
b=sym(b);
A\b