function W = lambertWhybrid(a,b)
% Lambert W function W{x} = W{a*exp(b)}, where x>=0
% Hybrid calculation formula
% Stratis Batzelis - 5 Mar 2019

% Initialization
logx = log(a)+b;
log9 = log(9);

% Asymptotic expansion for x>=9
L1 = logx(logx>=log9);
L2 = log(L1);
W(logx>=log9) = L1 -L2 +L2./L1.*( 1 +( -2+L2 +( 6-9*L2+2*L2.*L2 +( -12+36*L2-22*L2.*L2+3*L2.*L2.*L2 +(60-300*L2+350*L2.*L2-125*L2.*L2.*L2+12*L2.*L2.*L2.*L2)./(5*L1) )./(2*L1))./(3*L1) )./(2*L1));
% Series expansion for 0<=x<9
u = exp( logx(logx<log9) -1);
p = 1-u;
r = (1+u).*(1+u);
W(logx<log9) = u + u./(1+u).*p.*(1 +p./r/2.*(1 + p./r/3.*(-2*u+1 +p./r/4.*(6*u.*u-8*u+1 +p./r/5.*(-24*u.*u.*u+58*u.*u-22*u+1) ) ) ) );

end








