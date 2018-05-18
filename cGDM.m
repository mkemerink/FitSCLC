function MuEnhance = cGDM(Par,mu0,c,Fred)
%{
Implementation of the Bouhassoune parametrization of a density- and
field-dependent mobility for nearest neighbor hopping in a correlated
Gaussian DOS.

Commented-out statements are shown for clarity, parameters are set in the
calling routine ODDD.

Equation numbers refer to Bouhassoune et al., Org. El. 10 (2009), 437-445

Linköping University
Martijn Kemerink, May 18, 2018
%}

% Par.c1 = 0.49; %[1] in paper c1=1.0e-9=0.49*exp(-20)
% Par.c2 = 1; %[1] in paper c2=exp(-20)
% Par.Chi = 0.29; %[1] prefactor
% Par.mu_star = In.nu0*In.aNN^2/In.sigma; %[m2/Vs] (called mu0 in paper)
% Par.Q = 2.4/(1-Par.sigma_hat); %Eq. 12
% In.mu0 = Par.c1*Par.mu_star*exp(-Par.Chi*Par.sigma_hat^2); %[m2/Vs] Eq. 9
% Par.delta = 2.3*(log(0.5*Par.sigma_hat^2+1.4*Par.sigma_hat)-0.327)/Par.sigma_hat^2; %Eq. A3
% Par.r = 0.7*Par.sigma_hat^-0.7; %Eq. A7

%high-field regime
mu_hi = Par.c2./Fred*Par.mu_star.*(1-c); %[m2/Vs] Eq. 10

%low-field regime (note that Fred_star=0.16, above Eq. A6)
h = ones(size(Fred)); %below Eq. A7
use_h2 = Fred<=0.16/2; %Eq. A8
h2 = (4/3)*Fred/0.16; %Eq. A8
h(use_h2) = h2(use_h2);
use_h3 = (Fred>0.16/2) & (Fred<=0.16); %Eq. A9
h3 = 1-(4/3)*((Fred/0.16)-1).^2; %Eq. A9
h(use_h3) = h3(use_h3);
f = exp(h.*(1.05-1.2*c.^Par.r)*(Par.sigma_hat^(3/2)-2).*(sqrt(1+2*Fred)-1)); %Eq. A6

c_min = min(0.025,c); %Eq. A2
g = exp((0.25*Par.sigma_hat.^2+0.7*Par.sigma_hat).*(2*c_min).^Par.delta); %Eq. A1

mu_lo = mu0*g.*f; %[m2/Vs] Eq. 8

%total:
mu_tot = (mu_lo.^Par.Q + mu_hi.^Par.Q).^(1/Par.Q); %[m2/Vs] Eq. 11
MuEnhance = mu_tot/mu0; %[1] ODDD expects an enhancement factor

end