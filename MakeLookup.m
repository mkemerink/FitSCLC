function In = MakeLookup(In,dataT,Vmax)
%{
Makes a lookup table for mobility enhancement factor vs concentration and
(effective) temperature in case of ET-GDM model in which the field is
mapped on an increased effective temperature of the charge carrier
population.

The routine consists of near-identical copies of corresponding code in
ODDD and ET_GDM_old.

Linköping University
Martijn Kemerink, May 18, 2018
%}

function conc = CalcC(sHat,sHat2,H,K,E_Fermi)
%{
Calculates concentration for given EF in Gaussian DOS, uses Paasch JAP 107,
104501 (2010), i.e. the same method as used in ODDD for Gaussian traps.
%}
zeta = E_Fermi; %[1] eq.3, E_Fermi is already in units of kT
degen = zeta>-sHat2; %true implies degenerate limit
Gnondeg = exp(0.5*sHat2+zeta)./(exp(K*(zeta+sHat2))+1);
Gdeg = 0.5*erfc(-H*zeta/(1.41421*sHat));
conc = ~degen.*Gnondeg+degen.*Gdeg;
end

% physical constants (see ODDD)
Tref = min(dataT); %[K] reference T
Par.q = 1.602176487e-19; %[Coulomb]
Par.kB = 1.380650324e-23; %[J/K]
Par.E_th = Par.kB*Tref/Par.q; %[eV] thermal energy

% mobility stuff (see ODDD)
In.N0 = (In.aNN)^-3; %[m-3] density of states
Par.sigma_hat = In.sigma/Par.E_th; %[1] reduced disorder
Par.mu_star = In.nu0*In.aNN^2/In.sigma; %[m2/Vs]
Par.B = 0.47; %[1] from Cottaar, fcc + MA rates (cf c1 in other models)
Par.Ecrit =-0.491*Par.sigma_hat; %[1] from Cottaar, fcc + MA rates
Par.gamma = 0.67; %[1] prefactor in Shklovskii-expression for Teff
Par.lambda = 0.875; %[1] critical exponent (0.97 in Cottaar for fcc/MA)


Fnorm = 1.5*(Vmax/In.L)*In.alpha*Par.q/Par.kB;
Teff_max = sqrt(max(dataT)^2+(Par.gamma*Fnorm).^2); %[K] allways >=T

%Init T, c axes: extended to enable 2nd derivative
In.LookUp.minT = min(dataT); %[K] lowest (effective) T to consider
In.LookUp.maxT = Teff_max; %[K] highest Teff, higher seems unphysical (1000)
In.LookUp.nT = max(1,ceil((In.LookUp.maxT-In.LookUp.minT)/10)+1); %[1] nr of T points, 10 K interval is empirical
In.LookUp.preT = (In.LookUp.nT-1)/(In.LookUp.maxT-In.LookUp.minT); %prefactor to find index
dT = 1/In.LookUp.preT;
In.LookUp.T = linspace(In.LookUp.minT-dT,In.LookUp.maxT+dT,In.LookUp.nT+2); %[K] T-axis of lookup table

In.LookUp.minC =-8; %[1] log10 of lowest concentration
In.LookUp.maxC = -1; %[1] log10 of highest concentration
In.LookUp.nC = 10*round(In.LookUp.maxC-In.LookUp.minC)+1; %10 points per decade is empirical
In.LookUp.preC = (In.LookUp.nC-1)/(In.LookUp.maxC-In.LookUp.minC);
dC = 1/In.LookUp.preC;
In.LookUp.c = linspace(In.LookUp.minC-dC,In.LookUp.maxC+dC,In.LookUp.nC+2)'; %[K] concentration-axis of lookup table


%Init
In.LookUp.MuTable = zeros(In.LookUp.nC+2,In.LookUp.nT+2); %[1]
relT = Tref./In.LookUp.T;
sh = Par.sigma_hat*relT; %[1] shorthand to save flops
sh2 = sh.^2; %[1] shorthand to save flops
Hs = 1.02964-0.50216./sh-0.12889./sh2; %eq.14
Ks =-0.08566+1.65879./sh-0.74604./sh2; %eq.14

%find Fermi energy for each c, T combination and calc mu
for iT=1:In.LookUp.nT+2 %loop over (effective) T
    mu0_T = Par.B*Par.mu_star*sh(iT).^(1-Par.lambda).*...
        exp(-0.5*sh2(iT) - Par.Ecrit); %[m2/Vs] Eq. 5a
    EF0 =-sh2(iT); %trial EF=-sigma^2/kT in units of kT
    for iC=1:In.LookUp.nC+2 %loop over concentration
        f = @(x,sh_k,sh2_k,Hs_k,Ks_k)...
            log10(CalcC(sh_k,sh2_k,Hs_k,Ks_k,x))-In.LookUp.c(iC);
        sh_k = sh(iT);
        sh2_k = sh2(iT);
        Hs_k = Hs(iT);
        Ks_k = Ks(iT);
        fun = @(x)f(x,sh_k,sh2_k,Hs_k,Ks_k);
        EF = fzero(fun,EF0); %[1]
        EF0 = EF; %use this EF as trial for next concentration
        In.LookUp.MuTable(iC,iT) = log10( mu0_T*10^-In.LookUp.c(iC)*exp(EF + 0.5*sh2(iT)) );
    end
end

%derivatives for interpolation
ce = In.LookUp.MuTable;
up = circshift(In.LookUp.MuTable,-1,1);
do = circshift(In.LookUp.MuTable,1,1);
le = circshift(In.LookUp.MuTable,-1,2);
ri = circshift(In.LookUp.MuTable,1,2);

In.LookUp.dc = 0.5*(up-do)*In.LookUp.preC;
In.LookUp.dT = 0.5*(le-ri)*In.LookUp.preT;
In.LookUp.ddc = 0.5*(up-2*ce+do)*In.LookUp.preC^2;
In.LookUp.ddT = 0.5*(le-2*ce+ri)*In.LookUp.preT^2;
%keep middle part
In.LookUp.MuTable = In.LookUp.MuTable(2:end-1,2:end-1);
In.LookUp.dc = In.LookUp.dc(2:end-1,2:end-1);
In.LookUp.dT = In.LookUp.dT(2:end-1,2:end-1);
In.LookUp.ddc = In.LookUp.ddc(2:end-1,2:end-1);
In.LookUp.ddT = In.LookUp.ddT(2:end-1,2:end-1);
In.LookUp.T = In.LookUp.T(2:end-1);
In.LookUp.c = In.LookUp.c(2:end-1);

end