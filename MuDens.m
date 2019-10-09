function mu = MuDens(In)
%{
Note to self: Essentially this is a stripped-down version of the function
MuDens_fun in the directory MuDens. See there for a heading with full
details. Mind, however, that some parameters have changed name for
consistency with SCLCfit.

Martijn Kemerink, LiU, October 03, 2019
%}

%% Moved user input
nRough = 51; %[1] nr of points in array to find rough hopping optimum (51)

%% Setting of constants & derived input
q = 1.60217646263e-19; %[Coulomb]
kB = 1.380650324e-23; %[J/K]

%% Preallocation of output arrays
mu = zeros(In.LookUp.nC+2,In.LookUp.nT+2); %[m2/Vs] effective mobility

%% actual program

for ic = 1:In.LookUp.nC+2 %concentration loop
    c = 10^In.LookUp.c(ic); %In.LookUp.c = log10(c)
    
    for iT = 1:In.LookUp.nT+2 %temperature loop
        
        E_th = kB*In.LookUp.T(iT)/q; %[eV] thermal energy
        mu0 = In.B*In.aNN^2*In.nu0/E_th*(E_th/In.sigma)^In.lambda1*(In.alpha/In.aNN)^In.lambda2; %[m2/Vs] mobility prefactor for VRH
        
        EF = In.LookUp.EF(ic,iT)*E_th; %[eV] Fermi energy
        
        %find most probable hop and corresponding mobility etc.
        Efinal = linspace(EF,2*In.sigma,nRough); %[eV] final state energies
        n_dE0 = erf(EF/(In.sigma*sqrt(2)));
        n_dE = In.N0*0.5*(erf(Efinal/(In.sigma*sqrt(2)))-n_dE0); %[m-3] nr of states between EF and Efinal
        R = (In.Bc./((4/3)*pi*n_dE)).^(1/3); %[m] get R* from (4/3)pi x R*^3 x n_dE = Bc
        R(R==Inf) = 1; %[m]
        %corresponding hopping probability p according to Millar/Abrahams:
        p = exp(-2*R/In.alpha).*exp(-(Efinal-EF)/E_th); %[1]
        %find final site energy as maximum of hopping probability (i.e. at zero derivative)
        dp = diff(p); %[1]
        dp = ([dp(1) dp] + [dp dp(end)])/2; %derivative, element addition to keep same size
        rough_max = find(p==max(p)); %index of approximate maximum
        if rough_max==nRough %safety measure
            E_star = Efinal(end);
        else
            deriv = @(Ei) pchip(Efinal,dp,Ei); %use interpolation for higher accuracy
            E_star = fzero(deriv,[Efinal(rough_max-1) Efinal(rough_max+1)]); %[eV] E*
        end
        %mu with aNN as length scale in the prefactor
        mu(ic,iT) = mu0*pchip(Efinal,p,E_star)/c; %[m2/Vs] mobility
%         %mu with R as length scale - old method for mobility (assumes <R^2> = <R>^2
%         R_star = pchip(Efinal,R,E_star); %[m] R* hopping distance
%         mu(ic,iT) = mu0*(R_star/In.aNN)^2*pchip(Efinal,p,E_star)/c; %[m2/Vs] mobility
%         % - new method for mobility, based on <R^2> (SLOW and little difference)
%         num = integral(@(Ei)pchip(Efinal,R.^2.*p,Ei),Efinal(1),Efinal(end));
%         den = integral(@(Ei)pchip(Efinal,p,Ei),Efinal(1),Efinal(end));
%         R2 = num/den;
%         mu(ic,iT) = mu0*(R2/In.aNN^2)*pchip(Efinal,p,E_star)/c; %[m2/Vs] mobility
    end %end of temperature loop
    
end %end of concentration loop

mu = exp(2*In.aNN/In.alpha)*mu; %[m2/Vs] normalize to tunneling probability to nearest neighbor site

end