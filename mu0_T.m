function mu0 = mu0_T(muStar,C,E0,T)
%{
Function that generalizes the T-dependence of mu0 in the GDM and Arrhenius
models.
GDM:        mu0(T) = muStar*exp(-(C*sigma/kT)^2)
Arrhenius:  mu0(T) = muStar*exp(-Eact/kT),

generic:    mu0(T) = muStar*exp(-(C(1)*E0/kT)^C(2))

Linköping University
Martijn Kemerink, May 18, 2018
%}

kT = 1.380650324e-23*T/1.602176487e-19; %[eV] thermal energy
mu0 = muStar*exp(-(C(1)*E0./kT).^C(2)); %[m2/Vs] zero field mobility

end