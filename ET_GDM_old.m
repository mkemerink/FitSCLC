function MuEnhance = ET_GDM(Par,mu0,T,c,Fnorm)
%{
Implementation of the effective-temperature version of a density- and
field-dependent mobility for nearest neighbor hopping in an uncorrelated
Gaussian DOS.

Commented-out statements are shown for clarity, parameters are set in the
calling routine (ODDD or ...).

Linköping University
Martijn Kemerink, April 26, 2018
%}

function conc = CalcC(sHat,sHat2,H,K,E_Fermi)
%{
Calculates concentration for given EF in Gaussian DOS, uses Paasch JAP 107,
104501 (2010), i.e. the same method as used in ODDD for Gaussian traps.
%}

zeta = E_Fermi; %[1] eq.3 INTEGRATE IN EQ FOR SPEED
degen = zeta>-sHat2; %true implies degenerate limit
Gnondeg = exp(0.5*sHat2+zeta)./(exp(K*(zeta+sHat2))+1);
Gdeg = 0.5*erfc(-H*zeta/(1.41421*sHat));
conc = ~degen.*Gnondeg+degen.*Gdeg;
            
end

% Par.sigma_hat = In.sigma/Par.E_th; %[1] reduced disorder
% Par.mu_star = In.nu0*In.aNN^2/In.sigma; %[m2/Vs]
% Par.B = 0.47; %[1] from Cottaar, fcc + MA rates (cf c1 in other models)
% Par.Ecrit =-0.491*Par.sigma_hat; %[1] from Cottaar, fcc + MA rates
% Par.gamma = 0.67; %[1] prefactor in Shklovskii-expression for Teff
% Par.lambda = 0.875; %[1] critical exponent (0.97 in Cottaar for fcc/MA)
% In.mu0 = Par.B*Par.mu_star*Par.sigma_hat^(1-Par.lambda)*...
%     exp(-0.5*Par.sigma_hat^2-Par.Ecrit); %[m2/Vs] Eq. 5a
% Par.mu_pre = Par.B*Par.nu0/(Par.E_th*In.aNN)*Par.sigma_hat^Par.lambda; %

%Init
c = max(1e-14,c); %to prevent 1/0 errors (1e-8)
log_c = log10(c);
nPts = size(c,1);
EF = zeros(nPts,1); %[1] EF in units of kT

%effective T
% Fnorm = F*In.alpha*Par.q/Par.kB
Teff = sqrt(T^2+(Par.gamma*Fnorm).^2); %[K]
relT = T./Teff;

%Fermi energy
sh = Par.sigma_hat*relT; %[1] sigma_hat at Teff, shorthand to save flops
sh2 = sh.^2; %[1] sigma^2/kTeff in units of kTeff
Hs = 1.02964-0.50216./sh-0.12889./sh2; %eq.14
Ks =-0.08566+1.65879./sh-0.74604./sh2; %eq.14

for k=1:nPts
    f = @(x,sh_k,sh2_k,Hs_k,Ks_k) log10(CalcC(sh_k,sh2_k,Hs_k,Ks_k,x))-log_c(k);
    sh_k = sh(k);
    sh2_k = sh2(k);
    Hs_k = Hs(k);
    Ks_k = Ks(k);
    fun = @(x)f(x,sh_k,sh2_k,Hs_k,Ks_k);
    EF0 =-sh2(k); %trial EF=-sigma^2/kT in units of kTeff
    EF(k) = fzero(fun,EF0); %[1] EF in units of kTeff
end

% mu = Par.mu_pre*relT.^(Par.lambda+1).*exp((EF-Par.Ecrit).*relT)./(In.N0*c); ????
mu0_T = Par.B*Par.mu_star*(Par.sigma_hat.*relT).^(1-Par.lambda).*...
    exp(-0.5*(Par.sigma_hat.*relT).^2 - Par.Ecrit.*relT); %[m2/Vs] Eq. 5a
mu = mu0_T./c.*exp(EF + 0.5*(Par.sigma_hat.*relT).^2);

MuEnhance = mu/mu0; %[1] ODDD expects an enhancement factor

end