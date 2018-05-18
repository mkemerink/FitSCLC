function [Error] = FitSCLCerror(p,DATA,In,par,C)
%{
FitSCLCerror is called by fmincon that is invoked in the main program.
It...
1) Converts fit parameters in p back to 'intelligable' entities
2a) Runs the function SCLC to calculate the model jV according to p
(constrained fit) or starts a nested optimization of gamma and mu (free
fit) or
2b) Runs the function ODDD to calculate the model jV using drift-diffusion
3) Calculates the error as the summed squared differences between the jV in
DATA and model. The 10log of j's is used to properly cover the full jV.

Linköping University
Martijn Kemerink, May 18, 2018
%}

%Init/transfer of parameters
if In.FitModel==1 %Murgatroyd/Gill
    Vbi = p(1);
    if In.Constrained %constrained fit
        muStar = 10^p(2); %p(2-5) do not exist for FitType 2
        E0 = p(3); %sigma or E_act
        C = In.C; %prefactor and exponent in mu0(T)
        gammaB = 10^p(4);
        gammaT0 = p(5);
    end
else %drift-diffusion + e/c/ET-GDM
    In.phi0 = [p(1),p(2)];
    In.sigma = p(3);
    In.nu0 = 10^p(4);
    In.aNN = 1e-9*p(5);
    if In.muModel==3 %ET-GDM
        In.alpha = 1e-9*p(6);
    end
end

%Init lookup table in case of ET-GDM
if In.FitModel==2 && In.muModel==3
    Vmax = zeros(DATA.NrTemp,1);
    for k = 1:DATA.NrTemp %stupid loop to find Vmax...
        In.bias0 = DATA.jVT(DATA.ind(:,k),1,k);
        Vmax(k) = max(In.bias0);
    end
    In.L = min(DATA.L); %worst case, L is overruled below
    In = MakeLookup(In,DATA.T,max(Vmax));
end

%jV calculation
Error = 0;
for k = 1:DATA.NrTemp %loop over all loaded measurements
    j0 = DATA.jVT(DATA.ind(:,k),2,k); %[A/m2] current to fit
    In.bias0 = DATA.jVT(DATA.ind(:,k),1,k); %[V] corresponding applied bias
    In.T = DATA.T(k); %[K]
    In.L = DATA.L(k); %[m]
    
    if In.FitModel==1 %Murgatroyd/Gill
        if In.Constrained %constrained fit
            gamma = Gill(gammaB,gammaT0,In.T,In.MinGamma);
            mu0 = mu0_T(muStar,C,E0,In.T);
            j = SCLC(mu0,gamma,In.L,In.epsR,In.bias0-Vbi); %i.e. SCLC(mu0,gamma,L,epsR,effective bias)
        else %free fit: run a fit of mu and gamma for each T
            [j,~,~] = FreeFit(In.L,In.epsR,In.bias0-Vbi,In.T,In.MinGamma,j0,par,C); %[j,mu,gamma], last 2 not needed here
        end
    else %drift-diffusion + e/c/ET-GDM
        if In.SpeedMode %work on reduced bias axis
            Vmin = min(In.bias0); Vmax = max(In.bias0);
            nV = ceil(In.OutputPar(3)*...
                log10(Vmax/Vmin)+1); %nr of output points
            In.bias = logspace(log10(Vmin),log10(Vmax),nV)'; %[V]
        else %full bias (all experimental data points in range)
            In.bias = In.bias0;
        end
        
        %run ODDD
        [~,Out] = ODDD(In); %[In,Out], In is not needed

        %extract j from Out structure
        if In.SpeedMode %reduced bias axis
            if sum(~isnan(Out.j))>=2 %protect against no solution
                j = interp1(In.bias,Out.j,In.bias0,'pchip');
            else
                j = NaN(size(In.bias0));
            end
        else %full bias (all experimental data points in range)
            j = Out.j;
        end
        
    end
    
    %update error
    Error = Error + sum( (log10(abs(j0)) - log10(abs(j))).^2 );
end
Error = Error/DATA.NrTemp; %allows comparing errors with different numbers of T

end

