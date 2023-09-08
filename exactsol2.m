clearvars, close all

% Finding the exact solution with syms tool.
syms ro(x) a(x) f(x) mu alfa
Dro=diff(ro,1);
f=1;
a=mu*exp(x*alfa/mu);
roex=dsolve(-diff((a*Dro),1)==f, ro(0)==0, ro(1)==0);
roex=simplify(roex);
roex=vectorize(roex);