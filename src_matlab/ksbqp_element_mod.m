function [Ke, Ae] = ksbqp_element_mod(coords, t, E, nu, xi, w)
    % Calcola Ke e Ae usando mapping isoparametrico con shape-func del plate

    % rigidità materiale
    D = E * t^3 / (12 * (1 - nu^2)) * ...
        [1, nu, 0;
         nu, 1, 0;
         0,  0, (1 - nu) / 2];

    % matrice di cambio base C
    C = zeros(12);
    for i = 1:4
        Pi = P_matrix(coords(i,1), coords(i,2));
        C((i-1)*3 + (1:3), :) = Pi;
    end
    Cinv = inv(C);

    % Gauss
    [Xi, Eta] = meshgrid(xi, xi);
    gps    = [Xi(:), Eta(:)];
    weights = kron(w, w);

    K0 = zeros(12);
    Ae = 0;

    for gp = 1:length(weights)
        xi_i  = gps(gp,1);
        eta_i = gps(gp,2);
        wgt   = weights(gp);

        % ottengo shape + derivate parametriche
        [Nw_i, Ntx_i, Nty_i, dNw_dxi_i, dNw_deta_i] = ...
            shape_functions_at_gauss_points(xi_i, eta_i);

        % Jacobiana isoparametrica
        J = [ dNw_dxi_i  * coords(:,1),   dNw_deta_i * coords(:,1);
              dNw_dxi_i  * coords(:,2),   dNw_deta_i * coords(:,2) ];
        detJ = det(J);
        Ae = Ae + detJ * wgt;

        % Q‐matrix
        Q = zeros(3,12);
        for a = 1:4
            Q(1, (a-1)*3 + 2) = Nw_i(a);
            Q(2, (a-1)*3 + 3) = Nw_i(a);
            Q(3, (a-1)*3 + 2) = Nty_i(a);
            Q(3, (a-1)*3 + 3) = Ntx_i(a);
        end

        K0 = K0 + (Q' * D * Q) * detJ * wgt;
    end

    % trasformo nella base fisica
    Ke = Cinv' * K0 * Cinv;
end
