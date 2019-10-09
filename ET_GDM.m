function MuEnhance = ET_GDM(In,gamma,c,Fnorm)
%{
Implementation of the effective-temperature version of a density- and
field-dependent mobility for VRH hopping on either a regular or a random
lattice with an uncorrelated Gaussian DOS.

This version works by quadratic interpolation using coefficients in a
look-up table since finding EF by brute force at each data point takes
forever. The look-up table is generated in MakeLookup.m and gives
log10(mobility) since mobility is a roughly exponential function of
density, field and temperature.

The value of In.mu0 is not really important since it is normalized out.

Linköping University
Martijn Kemerink, October 03, 2019
%}

%effective T, log(concentration)
% Fnorm = F*In.alpha*Par.q/Par.kB
Teff = sqrt(In.T^2+(gamma*Fnorm).^2); %[K] allways >=T
log_c = max(In.LookUp.minC,log10(c)); %protection against underflow
log_c = min(In.LookUp.maxC,log_c); %protection against overflow

%indices of nearest c,T points on mesh
iC = int16(In.LookUp.preC*(min(In.LookUp.maxC,log_c)-In.LookUp.minC)+1);
iT = int16(In.LookUp.preT*(min(In.LookUp.maxT,Teff) -In.LookUp.minT)+1);
ind = sub2ind([In.LookUp.nC,In.LookUp.nT],iC,iT);

%mobility
dc = log_c-In.LookUp.c(iC); %distance from nearest mesh point
dT = Teff-In.LookUp.T(iT)';
mu = In.LookUp.MuTable(ind) + ...
    In.LookUp.dc(ind) .*dc + ... 1st derivative
    In.LookUp.ddc(ind).*dc.^2 + ... 2nd derivative
    In.LookUp.dT(ind) .*dT + ...
    In.LookUp.ddT(ind).*dT.^2; %mind that mu table is log10(mu)
MuEnhance = 10.^mu/In.mu0; %[1] ODDD expects an enhancement factor

end