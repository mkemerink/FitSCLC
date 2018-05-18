function gamma = Gill(gammaB,gammaT0,T,lowest)
%{
Empirical Gill expression for gamma(T) = B(1/kT-1/kT0)

Parameter lowest enables enforcing gamma>=0 or not.

Linköping University
Martijn Kemerink, May 18, 2018
%}

kT = 1.380650324e-23*T/1.602176487e-19; %[eV] thermal energy
gamma = max(lowest, gammaB./kT.*(1-T/gammaT0) ); %[eV(V/m)^-0.5] field enhancement factor

end