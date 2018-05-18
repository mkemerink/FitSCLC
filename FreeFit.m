function [jFit,mu0Fit,gammaFit] = FreeFit(L,epsR,bias,T,MinGamma,j0,par,C)
%{
Performs a free fit of mu0 and gamma to best match j0(bias) at the given
temperature T.

This must be a separate function to allow running it from both FitSCLC and
FitSCLCerror

Linköping University
Martijn Kemerink, May 18, 2018
%}

%% help function
    function error = FreeFitError(p)
        %has the same function as FitSCLCerror has for constrained fitting.
        mu0 = 10^p(1);
        gamma = max(MinGamma,p(2)); %prevent unallowed gamma<0 from being used
        j = SCLC(mu0,gamma,L,epsR,bias);
        error = sum( (log10(abs(j0)) - log10(abs(j))).^2 );
    end

%% main function

%Set p0 with initial (guess) values; use log to get values O(1)
p0 = [log10( mu0_T(10^par(2,2),C,par(3,2),T) ) ... mobility
    Gill(10^par(4,2),par(5,2),T,MinGamma)]; %field enhancement factor gamma

%Fit
p = fminsearch(@(p0) FreeFitError(p0),p0); %unconstrained=faster

mu0Fit = 10^p(1);
gammaFit = p(2);
jFit = SCLC(mu0Fit,gammaFit,L,epsR,bias);

end