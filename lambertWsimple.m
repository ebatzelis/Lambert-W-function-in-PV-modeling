function W = lambertWsimple(a,b)
% Lambert W function W{x} = W{a*exp(b)}, where x>=2
% Simple approximation formula
% Stratis Batzelis - 5 Mar 2019

% Initialization
logx = log(a)+b;
log2 = log(2);
L = logx(logx>=log2);

% Compute Lambert W{x} for x>=2
W(logx>=log2) = L -L.*log(L)./(L+1);
% Set W{x}=NaN for x<2
W(logx<log2) = NaN;

end





