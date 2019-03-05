function W = lambertWasymp4(a,b)
% Lambert W function W{x} = W{a*exp(b)}, where x>=3
% Asymptotic formula - 4 terms
% Stratis Batzelis - 5 Mar 2019

% Initialization
logx = log(a)+b;
log3 = log(3);
L1 = logx(logx>=log3);
L2 = log(L1);

% Compute Lambert W{x} for x>=3
W(logx>=log3) = L1 -L2 +L2./L1.*( 1 + (-2+L2)./(2*L1));
% Set W{x}=NaN for x<3
W(logx<log3) = NaN;

end

