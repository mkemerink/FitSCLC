function MuEnhance = eGDM(sigma_hat,delta,c,Fred)
%{
Implementation of the Pasveer parametrization of a density- and
field-dependent mobility for nearest neighbor hopping in an uncorrelated
Gaussian DOS.

Commented-out statements are shown for clarity, parameters are set in the
calling routine ODDD.

Linköping University
Martijn Kemerink, May 18, 2018
%}

% mu0 = aNN^2*nu0/sigma; %called mu_star in the main code
% c1 = 1;
% c2 = 0.44;
% sigma_hat = sigma/Par.E_th;
% delta = 2*(log(sigma_hat^2-sigma_hat)-log(log(4)))/sigma_hat^2;

% mu0_T = mu0*c1*exp(-c2*sigma_hat^2); %called mu0 in the main code
c_enhance = exp(0.5*(sigma_hat^2-sigma_hat)*(2*c).^delta);
% Fred = F*In.aNN/In.sigma;
F_enhance = exp(0.44*(sigma_hat^1.5-2.2)*(sqrt(1+0.8*Fred.^2)-1));

% mu = c_enhance.*F_enhance*mu0;
MuEnhance = c_enhance.*F_enhance;

end