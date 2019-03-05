function W = lambertWanalyt(x)
% Lambert W function W{x}, where x>=3e-5
% Analytical approximation formula
% Stratis Batzelis - 5 Mar 2019

% Initialization
E = 0.4586887;
lb = 3e-5;
z = x(x>=lb);

% Compute Lambert W{x} for x>=3e-5
W(x>=lb) = (1+E)*log(6/5*z./log(12/5*z./log(1+12/5*z)) ) -E*log(2*z./log(1+2*z));
% Set W{x}=NaN for x<e3-5
W(x<lb) = NaN;

end
