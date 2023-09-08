clearvars, close all

% Section 1: computation of the exact solution.

syms u(x) psi(x) H(x) mu alfa
psi(x)=alfa*x;
Du=diff(u,1);
Dpsi=diff(psi,1);
H(x)=mu*Du-Dpsi*u;
f=1;
uex=dsolve(-diff(H,1)==f, u(0)==0, u(1)==0);
uex=simplify(uex);
uex=vectorize(uex);
