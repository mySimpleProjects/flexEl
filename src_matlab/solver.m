function [U, reactions] = solver(data)
    % SOLVER: assembla K, forze, applica BC e risolve K U = F
    % OUTPUT:
    %   U          : vettore spostamenti nodali [nDOF×1]
    %   reactions  : vettore reazioni globali [nDOF×1]

    % Assemblaggio
    [Kglob, Fglob, Ae] = assembleKSBQP( ...
        data.coords, data.conn, data.thickness, ...
        data.E, data.nu, data.xi, data.w, data.q );
    
    
    fprintf("\nArea: A = %.1f mm\n", sum(Ae));
    % Combino forze puntuali e distribuite
    F = data.BC.force + Fglob;

    % Impongo Dirichlet
    nDOF   = size(Kglob,1);
    U      = zeros(nDOF,1);
    fixed  = data.BC.dirichlet.dofs(:);
    U(fixed) = data.BC.dirichlet.values(:);

    % Risolvo sistema ridotto
    allDOF  = (1:nDOF).';
    freeDOF = setdiff(allDOF, fixed);
    Kff     = Kglob(freeDOF, freeDOF);
    Kfi     = Kglob(freeDOF, fixed);
    Ff      = F(freeDOF) - Kfi*U(fixed);
    U(freeDOF) = Kff \ Ff;

    % Calcolo reazioni su tutti i DOF: R = K*U - F
    reactions = Kglob*U - F;
end
