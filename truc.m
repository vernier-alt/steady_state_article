clear all
syms Omegax real;
syms Omegay real;
syms Omegaz real; 
syms T1 real positive;
syms T2 real positive;
syms X;
syms a b c;
syms x y z;

syms u(t) v(t) w(t)

A = [-1/T2-X,         -Omegaz,   Omegay;
    
    Omegaz,       -1/T2-X,        -Omegax;
    
    -Omegay,     Omegax,         -1/T1-X];

d = -(T1*Omegax^2*T2^2*X + T1*Omegax^2*T2 + T1*Omegay^2*T2^2*X + T1*Omegay^2*T2 + T1*Omegaz^2*T2^2*X + Omegaz^2*T2^2 + T1*T2^2*X^3 + T2^2*X^2 + 2*T1*T2*X^2 + 2*T2*X + T1*X + 1)/(T1*T2^2);
 
B = [-1/T2,         -Omegaz,   Omegay;
    
    Omegaz,       -1/T2,        -Omegax;
    
    -Omegay,     Omegax,         -1/T1];
E = [0,         -Omegaz,   Omegay;
    
    Omegaz,       0,        -Omegax;
    
    -Omegay,     Omegax,         0];

C = [0, 0, -x;
    -1, 0, -y;
    0, -1, -z];

ode1 = diff(u) == -x*w;
ode2 = diff(v) == -u -y*w;
ode3 = diff(w) == -v -z*w;
odes = [ode1; ode2;ode3];

cond1 = u(0) == 0;
cond2 = diff(v(0)) == 0;
conds = [cond1];


S = dsolve(odes,conds);
% S.w
%  
% ans =
%  
% C1*exp(t*root(#X^3 + #X^2*z - #X*y + x, #X, 1)) + C2*exp(t*root(#X^3 + #X^2*z - #X*y + x, #X, 2)) + C3*exp(t*root(#X^3 + #X^2*z - #X*y + x, #X, 3))
%  
 
r1 =                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ((((T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) + (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) - ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^2 - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))^3)^(1/2) - (T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) - (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) + ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^(1/3) + ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))/((((T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) + (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) - ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^2 - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))^3)^(1/2) - (T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) - (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) + ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^(1/3) - (T2^2 + 2*T1*T2)/(3*T1*T2^2);
r2 =  - (3^(1/2)*(((((T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) + (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) - ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^2 - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))^3)^(1/2) - (T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) - (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) + ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^(1/3) - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))/((((T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) + (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) - ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^2 - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))^3)^(1/2) - (T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) - (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) + ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^(1/3))*1i)/2 - ((((T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) + (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) - ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^2 - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))^3)^(1/2) - (T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) - (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) + ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^(1/3)/2 - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))/(2*((((T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) + (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) - ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^2 - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))^3)^(1/2) - (T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) - (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) + ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^(1/3)) - (T2^2 + 2*T1*T2)/(3*T1*T2^2);
r3 =   (3^(1/2)*(((((T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) + (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) - ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^2 - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))^3)^(1/2) - (T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) - (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) + ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^(1/3) - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))/((((T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) + (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) - ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^2 - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))^3)^(1/2) - (T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) - (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) + ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^(1/3))*1i)/2 - ((((T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) + (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) - ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^2 - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))^3)^(1/2) - (T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) - (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) + ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^(1/3)/2 - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))/(2*((((T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) + (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) - ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^2 - ((T2^2 + 2*T1*T2)^2/(9*T1^2*T2^4) - (T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1)/(3*T1*T2^2))^3)^(1/2) - (T2^2 + 2*T1*T2)^3/(27*T1^3*T2^6) - (T1*Omegax^2*T2 + T1*Omegay^2*T2 + Omegaz^2*T2^2 + 1)/(2*T1*T2^2) + ((T2^2 + 2*T1*T2)*(T1*Omegax^2*T2^2 + T1*Omegay^2*T2^2 + T1*Omegaz^2*T2^2 + 2*T2 + T1))/(6*T1^2*T2^4))^(1/3)) - (T2^2 + 2*T1*T2)/(3*T1*T2^2);
 
% s1 = C1*exp(t*(-r2))*(z*(-r2) - y + (-r2)^2) - exp(t*(-r1))*(C1*(-r2)^2 - C2*y - C1*y + C2*(-r3)^2 + C1*z*(-r2) + C2*z*(-r3)) + C2*exp(t*(-r3))*(z*(-r3) - y + (-r3)^2)
 
% s2 = (exp(t*(-r1))*(z + (-r1))*(C1*(-r2)^2 - C2*y - C1*y + C2*(-r3)^2 + C1*z*(-r2) + C2*z*(-r3)))/(z*(-r1) - y + (-r1)^2) - C2*exp(t*(-r3))*(z + (-r3)) - C1*exp(t*(-r2))*(z + (-r2))
 
% s3 =C1*exp(t*(-r2)) + C2*exp(t*(-r3)) - (exp(t*(-r1))*(C1*(-r2)^2 - C2*y - C1*y + C2*(-r3)^2 + C1*z*(-r2) + C2*z*(-r3)))/(z*(-r1) - y + (-r1)^2)
 
  % voir racine polynôme avec déterminant , rajouter une condition
  % initiale
 