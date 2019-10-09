function EF = findEf(E2,DOS2,E_th,c)
%determines the Fermi energy in the effective DOS defined by the lookup
%table E2,DOS2 (x,y)

    function dc = residue(Efermi)
        %calculates the difference between the charge carrier concentration
        %corresponding to Efermi and the actual concentration.
        
        function occ = Occupation(E)
            occ = pchip(E2,DOS2,E)./(1+exp((E-Efermi)/E_th));
        end
        
        dc = integral(@Occupation,Emin,Emax) - c; %[1] is zero for correct Fermi energy
    end

Emin = min(E2); %[eV] lower integration edge
Emax = max(E2); %[eV] upper integration edge

EF = fzero(@residue,0);

% residue(EF) %check the convergence

end