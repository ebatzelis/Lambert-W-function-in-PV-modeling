function W = lambertWasymp7(a,b)
% Lambert W function W{x} = W{a*exp(b)}, where x>=3
% Asymptotic formula - 7 terms
% Stratis Batzelis - 5 Mar 2019

% Initialization
logx = log(a)+b;
log3 = log(3);
L1 = logx(logx>=log3);
L2 = log(L1);

% Compute Lambert W{x} for x>=3
W(logx>=log3) = L1 -L2 +L2./L1.*( 1 +( -2+L2 +( 6-9*L2+2*L2.*L2 +( -12+36*L2-22*L2.*L2+3*L2.*L2.*L2 +(60-300*L2+350*L2.*L2-125*L2.*L2.*L2+12*L2.*L2.*L2.*L2)./(5*L1) )./(2*L1))./(3*L1) )./(2*L1));
% Set W{x}=NaN for x<3
W(logx<log3) = NaN;

end