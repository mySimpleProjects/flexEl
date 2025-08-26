%% ksbqp_element.m
function [Ke, Ae] = ksbqp_element(coords, t, E, nu, xi, w)
    % Compute element stiffness (Ke) and area (Ae) for 4-node plate
    D = E*t^3/(12*(1-nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    C = zeros(12);
    for i = 1:4
        Pi = P_matrix(coords(i,1), coords(i,2));
        C((i-1)*3 + (1:3), :) = Pi;
    end
    Cinv = inv(C);
    Ae = 0; K0 = zeros(12);
    for i = 1:length(xi)
      for j = 1:length(xi)
        xi_i = xi(i); eta_j = xi(j);
        wgt = w(i)*w(j);
        N_iso = 0.25*[(1-xi_i)*(1-eta_j);
                     (1+xi_i)*(1-eta_j);
                     (1+xi_i)*(1+eta_j);
                     (1-xi_i)*(1+eta_j)];
        x = N_iso'*coords(:,1);
        y = N_iso'*coords(:,2);
        dN = 0.25*[-(1-eta_j), -(1-xi_i);
                    (1-eta_j), -(1+xi_i);
                    (1+eta_j),  (1+xi_i);
                   -(1+eta_j),  (1-xi_i)];
        detJ = det(dN'*coords);
        Ae = Ae + detJ*wgt;
        Q = Q_matrix(x,y);
        K0 = K0 + (Q'*D*Q)*detJ*wgt;
      end
    end
    Ke = Cinv' * K0 * Cinv;
end