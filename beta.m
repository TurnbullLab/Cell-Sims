function b = beta(t,par)
%beta(t) is the probability per unit time that a gcp divides into 2 gcs
% t is time, measured in (postnatal) days.
% Td, lambda, tau are passed to function via the "par" structure:
% Td, lambda are constant, tau is variable: different for each clone.
%   Td     = dividing time of gcps, constant
%   lambda = "tuning constant" for sigmoidal function 
%   tau    = time of differentiation for a given clone, variable

Td=par.Td;
tau=par.tau;
lam=par.lambda;

arg=(t/tau)^lam;
b=(log(2)/Td)*arg/(1+arg);

end
