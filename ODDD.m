function [In,Out] = ODDD(In)

%{
Function ODDD (One-Dimensional Drift-Diffusion) calculates what the name
promises.
This version is for single carrier devices and has been stripped
down a bit for use with the SCLC fitting GUI. The calculation controls are
taken from the table In.DDpar.

To include trapping effects in the calculation, select trapModel 1-3 (don't
forget to switch back to 0 when done!) and set the corresponding
parameters. These are not optimized in the fitting. For doping, just set
doping>0. You can change parameters while the GUI is running, just press
'Save'.

Equations are described in the document 'Equations for ODDD.docx'.

Bias is applied on left contact (x=0), the right contact (x=L) is grounded.
Potentials and energies are in eV (e=-q with q the positive elementary
charge)

In case of non-convergence, using a bit lower In.mix can help. The image
potential modification of the contacts lowers the stability a bit, and can
be cut out for ohmic contacts

Linköping University
Martijn Kemerink, May 18, 2018
%}

%% Initialization

% some further device parameters
In.Egap = 100; %[eV] effective bandgap
% In.muModel = 2; %2=Pasveer parametrization of eGDM;
In.trapModel = 0; %0=none; 1=single trap; 2=Gaussian; 3=exponential
    %general & trapModel=1 - single shallow trap in Boltzmann limit
    In.Nt = 5e22; %[m-3] total density of trap states (1e23|-|1e25))
    In.Et = 0.10; %[eV] trap depth, >0 (0.15|0.55)
    %trapModel=2 - Gaussian trap DOS
    In.sigmat = 0.1; %[eV] trap width
    %trapModel=3 - exponential trap DOS (Et not used)
    In.Tt = 900; %[K] characteristic temperature Tt>T
In.doping = 0; %[m-3] doping concentration

% physical constants
Par.q = 1.602176487e-19; %[Coulomb]
Par.eps0 = 8.854187816e-12; %[]
Par.kB = 1.380650324e-23; %[J/K]
Par.E_th = Par.kB*In.T/Par.q; %[eV] thermal energy
Par.C = Par.eps0*In.epsR/In.L; %[C/Vm2] areal capacitance
In.Vbi = In.phi0(1)-In.phi0(2); %[eV] built-in voltage

% calculation control
In.dx = In.DDpar(1); %[m] mesh size (0.1...1e-9)
In.MinIt = In.DDpar(2); %[1] min nr of iterations in self consistent loop (5)
In.MaxIt = In.DDpar(3); %[1] max nr. of iterations in self-consistent loop (200)
In.dxCont = In.DDpar(4); %[m] thickness of contact region (5e-9)
In.Qexcess = In.DDpar(5); %[1] factor for exceeding SCLC charge concentration (10)
In.mix = In.DDpar(6); %[1] weight of newest solution in mixing with old (0.1...0.2)
In.MaxChange = In.DDpar(7); %[1] convergence criterion: max relative change in potential (1e-5)
In.UseImPot = In.DDpar(8); %[1] toggles use of image potential; false is more stable (true)

% mobility stuff
In.N0 = (In.aNN)^-3; %[m-3] density of states
Par.sigma_hat = In.sigma/Par.E_th; %[1] reduced disorder
Par.mu_star = In.nu0*In.aNN^2/In.sigma; %[m2/Vs]
if In.muModel==1 %eGDM (Pasveer) - uncorrelated Gaussian DOS
    Par.c1 = 0.87; %[1] in paper c1=1.8e-9=0.87*exp(-20)
    Par.c2 = 0.44; %[1] prefactor (also called Chi)
    Par.delta = 2*(log(Par.sigma_hat^2-Par.sigma_hat)-log(log(4)))/Par.sigma_hat^2;
    In.mu0 = Par.c1*Par.mu_star*exp(-Par.c2*Par.sigma_hat^2); %[m2/Vs]
elseif In.muModel==2 %cGDM (Bouhassoune) - correlated Gaussian DOS
    Par.c1 = 0.49; %[1] in paper c1=1.0e-9=0.49*exp(-20)
    Par.c2 = 1; %[1] in paper c2=exp(-20)
    Par.Chi = 0.29; %[1] prefactor
    Par.Q = 2.4/(1-Par.sigma_hat); %Eq. 12
    In.mu0 = Par.c1*Par.mu_star*exp(-Par.Chi*Par.sigma_hat^2); %[m2/Vs] Eq. 9
    Par.delta = 2.3*(log(0.5*Par.sigma_hat^2+1.4*Par.sigma_hat)-0.327)/Par.sigma_hat^2; %Eq. A3
    Par.r = 0.7*Par.sigma_hat^-0.7; %Eq. A7
elseif In.muModel==3 %ET-GDM on lattice (Cottaar PRL and Baranovskii for Teff)
    Par.B = 0.47; %[1] from Cottaar, fcc + MA rates (cf c1 in other models)
    Par.Ecrit =-0.491*Par.sigma_hat; %[1] from Cottaar, fcc + MA rates
    Par.gamma = 0.67; %[1] prefactor in Shklovskii-expression for Teff
    Par.lambda = 0.875; %[1] critical exponent (0.97 in Cottaar for fcc/MA)
    In.mu0 = Par.B*Par.mu_star*Par.sigma_hat^(1-Par.lambda)*...
        exp(-0.5*Par.sigma_hat^2-Par.Ecrit); %[m2/Vs] Eq. 5a in Cottaar
    Par.mu_pre = Par.B*In.nu0/(Par.E_th*In.aNN)*Par.sigma_hat^Par.lambda; %
elseif In.muModel==4 %ET-GDM on random grid
    Par.gamma = 0.67; %[1] prefactor in Shklovskii-expression for Teff
    Par.Ecrit =-0.491*Par.sigma_hat; %[1] from Cottaar, fcc + MA rates
    %rough estimate for mu0, the value is not important as it gets overruled in ET_GDM
    In.mu0 = (In.aNN^2*In.nu0/Par.E_th)*...
        exp(-(Par.Ecrit-In.sigma^2/Par.E_th)); %[m2/Vs] mobility prefactor for VRH
end
In.D0 = In.mu0*Par.E_th; %[m2/s] diffusion constant

% meshing + BC
In.nV = int32(size(In.bias,1)); %nr of points in j(V) curve
In.nx = round(In.L/In.dx+1); %nr of mesh points
In.MaxIt = int32(In.MaxIt); %int32() is for faster indexing
In.nxCont = round(In.dxCont/In.dx); %[1] thickness of contact region in unit cells
Out.mesh = (0:In.dx:double(In.nx-1)*In.dx); %[m] x-mesh

% convergence
Out.Qfactor = zeros(In.MaxIt,1); %[1] reduction factor of space charge density
Out.RelChange = zeros(In.MaxIt,1); %[1] relative change in V per iteration
Out.converged = false(In.nV,1); %[1] set true if convergence is reached

% miscellanious
Out.j = zeros(In.nV,1); %[A/m2] current density corresponding to input applied bias
Out.jSCLC = (9/8)*Par.eps0*In.epsR*In.mu0*abs(In.bias-In.Vbi).^2*In.L^-3;

%% unit conversion (energy in units of kT, density in concentration)
Work.phi0 = In.phi0/Par.E_th; %[1] injection barriers
Work.Vbi = Work.phi0(1)-Work.phi0(2); %[1] built-in voltage
% Work.n0 = min(1,exp(-Work.phi0)); %[1] BC for charge concentrations
Work.bias = In.bias/Par.E_th; %[1] applied bias in kT
Work.Et = In.Et/Par.E_th; %[1] trap energy
Work.Nt = In.Nt/In.N0; %[1]
Work.doping = In.doping/In.N0; %[1] doping acts as positive background charge density

%various prefactors
Work.preJ = Par.q/In.dx*In.N0; %current calculation
Work.Qpar = In.dx*Par.q*In.N0; %to calculate total areal charge density
if In.muModel==1 || In.muModel==2
    Work.preF = Par.E_th/In.dx*In.aNN/In.sigma; %to SI * normalization of field
elseif In.muModel==3 || In.muModel==4
    Work.preF = Par.E_th/In.dx*In.alpha*Par.q/Par.kB; %to SI * normalization of field
end
if In.UseImPot
    Work.preIm = (1/Par.E_th/In.dx)*(Par.q/(4*pi*Par.eps0*In.epsR)); %Vim @ contacts
else %will give incorrect field dependence for non-Ohmic contacts!
    Work.preIm = 0;
end

%% trapping stuff

if In.trapModel==0 %no traps
    Par.preTrap = 0; %switches off trapping
    Par.expTrap = 0;
elseif In.trapModel==1 %single trap
    Par.preTrap = exp(Work.Et);
    Par.expTrap = 1;
elseif In.trapModel==2 %Gaussian trap DOS
    Par.preTrap = 1; %because of unit of Nt
    Par.sigma_hatt = In.sigmat/Par.E_th; %[1] reduced disorder of trap DOS
    Par.sh2 = Par.sigma_hatt^2; %[1] shorthand to save flops
    Par.H = 1.02964-0.50216/Par.sigma_hatt-0.12889/Par.sigma_hatt^2; %eq.14
    Par.K =-0.08566+1.65879/Par.sigma_hatt-0.74604/Par.sigma_hatt^2; %eq.14
else %exponential trap distribution
    Par.preTrap = 1;
    Par.expTrap = In.T/In.Tt; %i.e. 1/l with l=Tt/T
end

%% continuity equation + density matrices

%indices of elements for matrix corresponding to continuity equation
%C*n=U, n=concentration vector, U=effective generation vector
Ci  = (2:In.nx-1)'; %row index = column index of diagonal
CjL = (1:In.nx-2)'; %column index of lower diagonal
CjU = (3:In.nx)'; %column index of upper diagonal
Work.Ci = [1;Ci;Ci;Ci;In.nx];
Work.Cj = [1;CjL;Ci;CjU;In.nx];
Work.Cv = zeros(size(Work.Ci)); %to contain elements of continuity eq.
Work.Cv(1) = 1; Work.Cv(end) = 1; %[1] enables setting BC for concentration n

%indices of elements for matrix to calculate interpollated densities
%D*n=interpolated density
Ci  = (1:In.nx-1)'; %row index = column index of diagonal
CjU = (2:In.nx)'; %column index of upper diagonal
Work.Di = [Ci;Ci];
Work.Dj = [Ci;CjU];
Work.Dn = zeros(size(Work.Di)); %to contain elements of electron interpol.
Work.Dnv = zeros(size(Work.Di));
clear Ci CjL CjU

Work.D = zeros(In.nx-1,1); %[m2/s] diffusion constant vs position
Out.nf_SG = zeros(In.nx-1,1); %[m-3] interpolated electron density

%U = G-R = generation-recombination
Work.U = zeros(In.nx,1); %[s-1]

%% Poisson matrices
%{
P*V=Rho; 1st and last elements of Rho contain BC for potential V
Inversion is done using QR decomposition (from MatLab help):
If A*x = b then x = R\(R'\(A'*b))
This is preferred over x = inv(A)*b and equally fast
%}
P = spalloc(In.nx,In.nx,3*In.nx+2);
for i=2:In.nx-1
    P(i,i-1) = 1; P(i,i+1) = 1; %off-diagonal elements
    P(i,i) =-2; %diagonal
end
P = (Par.E_th/In.dx^2)*(Par.eps0*In.epsR/Par.q/In.N0)*P;
P(1,1) = 1; P(In.nx,In.nx) = 1; %BC at left, right contact
R = qr(P);
Pp = P';
Rp = R';   

Work.Rho = zeros(In.nx,1); %[1] init charge concentration (dope and free compensate)

%% Calculation

Work.nf = Work.doping*ones(In.nx,1);
Work.nt = zeros(In.nx,1);

for iBias=1:In.nV %loop over requested biases
    Out.Qfactor(:) = 0;
    Out.RelChange(:) = 0;
    %estimate total areal charge density in device
    Work.Qtyp = Par.C*max(1,abs(In.bias(iBias))); %[C/m2]
    %comment-out next line to generate trial V on basis of previous charge concentration
%     Work.Rho = zeros(In.nx,1);
    Work.Rho(1) =-(Work.bias(iBias)-Work.Vbi); Work.Rho(In.nx) = 0;%[1] boundary conditions
    Work.Vold = R\(Rp\(Pp*Work.Rho));
    Work.F = diff(Work.Vold); %[1/dx] field in kT per dx
    
    for It=1:In.MaxIt %self-consistency loop
        
        %BC for concentration from image force-modified injection barriers
        Work.phi = Work.phi0-sqrt(max(0,Work.preIm.*[-Work.F(1),Work.F(end)]));
        Work.n0 = min(1,exp(-Work.phi)); %[1] BC for charge concentrations
        Work.U(1) = Work.n0(1);  Work.U(In.nx) = Work.n0(2); %[1]

        %assemble matrices for continuity equation
        %use B(x) = x/(1-exp(x)) = x/2-(x/2)/tanh(x/2) = B1-B2 where
        %x=(mu/D)*(V(i+1)-V(i)) and use -x/tanh(-x)=x/tanh(x)
        Work.B1 =-0.5*diff(Work.Vold); %x/2 -sign is for eV
        hlp = 1./tanh(Work.B1); %1/tanh(x/2)
        Work.B2 = Work.B1.*hlp; %(x/2)/tanh(x/2)
        if In.muModel==0
            Work.D(:) = In.D0; %no mobility enhancement
        elseif In.muModel==1 %eGDM (Pasveer)
            Work.D(:) = In.D0*eGDM(Par.sigma_hat,Par.delta,...
                0.5*(Work.nf(1:end-1)+Work.nf(2:end)),Work.preF*Work.F);
        elseif In.muModel==2 %cGDM (Bouhassoune)
            Work.D(:) = In.D0*cGDM(Par,In.mu0,...
                0.5*(Work.nf(1:end-1)+Work.nf(2:end)),abs(Work.preF*Work.F));
        elseif In.muModel==3 || In.muModel==4 %ET-GDM (Cottaar/Baranovskii)
            Work.D(:) = In.D0*ET_GDM(In,Par.gamma,...
                0.5*(Work.nf(1:end-1)+Work.nf(2:end)),abs(Work.preF*Work.F));
        end
        Work.Bp = Work.D.*( Work.B1-Work.B2); %D*B(x)
        Work.Bm = Work.D.*(-Work.B1-Work.B2); %D*B(-x)
        Work.Cv(2:end-1) = In.dx^-2*...
            [-Work.Bm(1:end-1);Work.Bp(1:end-1)+Work.Bm(2:end);-Work.Bp(2:end)];
        Work.C = sparse(Work.Ci,Work.Cj,Work.Cv);
        %solve continuity, left divide solves Ax=B for x (x=A\B)
        Work.nf = abs(Work.C\Work.U); %[1] abs() suppresses small imaginary parts
        
        %calculate trapped charge densities
        if In.trapModel==2 %Gaussian trap, follows Paasch JAP 107, 104501 (2010)
            Work.EF = log(Work.nf); %[1] Fermi energy w.r.t. band
            zeta = Work.EF+Work.Et; %[1] eq.3
            degen = zeta>-Par.sh2; %true implies degenerate limit
            Gnondeg = exp(0.5*Par.sh2+zeta)./(exp(Par.K*(zeta+Par.sh2))+1);
            Gdeg = 0.5*erfc(-Par.H*zeta/(sqrt(2)*Par.sigma_hatt));
            Work.nt = Par.preTrap*(~degen.*Gnondeg+degen.*Gdeg);
        else %all other cases
            Work.nt = Par.preTrap*Work.nf.^Par.expTrap; %[1] in units of In.Nt
        end
        Work.nt = min(Work.nt,1); %avoids over-filling (esp. for single trap)
        
        %calculate total charge density
        Work.Rho(2:end-1) = -(Work.nf(2:end-1) + ... free charges +
            Work.Nt*Work.nt(2:end-1) - ... trapped charges -
            Work.doping); %doping, in units of N0
           
        %avoid overshoot due to too unphysical charge concentration
        Work.Qtot = Work.Qpar*sum(abs(Work.Rho(In.nxCont+1:end-In.nxCont))); %[C/m2]
        Out.Qfactor(It) = max(In.Qexcess,Work.Qtot/Work.Qtyp)/In.Qexcess; %[1] >1 means unrealistic Q
        Work.Rho(2:end-1) = (1/Out.Qfactor(It))*Work.Rho(2:end-1); %so reduce occupations accordingly

        %solve Poisson & make new V & F
        Work.Vnew = R\(Rp\(Pp*Work.Rho));
        Out.RelChange(It) = sum(abs(Work.Vnew-Work.Vold))/max(1,sum(abs(Work.Vold)));
        Work.Vold = (1-In.mix)*Work.Vold + In.mix*Work.Vnew;
        Work.F = diff(Work.Vold); %[1/dx] field in kT per dx
        
        %convergence
        hlp_RC = mean(Out.RelChange(max(1,It-In.MinIt+1):It)); %noise filter, averages last points
        hlp_Qf = mean(Out.Qfactor(max(1,It-In.MinIt+1):It)); %artificial charge reduction is not allowed
        if hlp_RC<In.MaxChange && hlp_Qf==1 && It>=In.MinIt
            Out.converged(iBias) = true;
            break
        end
        
    end %end of self-consistency loop
    if ~Out.converged(iBias)
        disp(['Calculation ', num2str(iBias), ' not converged. RelChange=',...
            num2str(hlp_RC), '; Qfactor=', num2str(hlp_Qf), '.'])
    end
    
    %assemble matrices for density interpolation
    %use 1/(1-exp(x)) = 1/2-(1/2)/tanh(x/2) and tanh(-x)=-tanh(x)
    Work.Dnv(:) = 0.5*[1+hlp;1-hlp];
    Work.Dn = sparse(Work.Di,Work.Dj,Work.Dnv);
    Out.nf_SG = Work.Dn*Work.nf; %[1] solve SG-interpolated electron density
    Out.jn = Work.preJ*Work.D.*Out.nf_SG.*Work.F; %[A/m2] electron j vs. x
    Out.j(iBias) = mean(Out.jn); %[A/m2]
    
end %end of loop over bias

%% post processing

%bring Work-arrays to SI units in Out structure
Out.V = Par.E_th*Work.Vold; %[V] electrostatic potential vs x
Out.nf = In.N0*Work.nf; %[m-3] free density vs x
Out.nt = In.Nt*Work.nt; %[m-3] trapped density vs x
Out.EF = Out.V+Par.E_th*log(Work.nf); %[eV] Fermi energy vs x

end