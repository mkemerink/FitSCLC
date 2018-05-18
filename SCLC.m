function j = SCLC(mu0,gamma,L,epsR,bias) %SCLC(mu0,gamma,In)
%{
Calculates j(bias) for the given parameters from the Murgatroyd expression
for space charge limited conduction. Uses SI units

Linköping University
Martijn Kemerink, May 18, 2018
%}

j = (9/8)*8.854187816e-12*epsR*mu0*L^-3*(bias.^2).*...   normal SCLC
    exp(0.891*gamma*sqrt(abs(bias/L))); %field enhancement

end