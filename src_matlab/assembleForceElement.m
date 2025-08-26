function Fe = assembleForceElement(coords, xi, w, q)
    % ASSEMBLEFORCEELEMENT: vettore forze nodali per carico superficiale costante q
    nGP = length(xi);
    Fe  = zeros(12,1);

    % Costruisci Cinv per trasformazione shape→dofs
    C = zeros(12);
    for i = 1:4
        C((i-1)*3 + (1:3), :) = P_matrix(coords(i,1), coords(i,2));
    end
    Cinv = inv(C);

    for i = 1:nGP
        for j = 1:nGP
            xi_i = xi(i); eta_j = xi(j);
            wt   = w(i)*w(j);

            % Jacobiano del mapping isoparametrico
            dN = 0.25 * [-(1-eta_j), -(1-xi_i);
                          (1-eta_j), -(1+xi_i);
                          (1+eta_j),  (1+xi_i);
                         -(1+eta_j),  (1-xi_i)];
            J    = dN' * coords;
            detJ = abs(det(J));

            % Shape isoparametriche
            N_iso = 0.25 * [(1-xi_i)*(1-eta_j);
                            (1+xi_i)*(1-eta_j);
                            (1+xi_i)*(1+eta_j);
                            (1-xi_i)*(1+eta_j)];
            % Corrispondenza (x,y)
            x = N_iso' * coords(:,1);
            y = N_iso' * coords(:,2);

            % Shape displacement→dof
            Pval = P_matrix(x,y);      % 3×12
            Ntot = Pval * Cinv;        % 3×12
            Nw   = Ntot(1,:);          % 1×12 → spostamento W

            % Accumula Fe
            Fe = Fe + (Nw' * q) * detJ * wt;
        end
    end
end
