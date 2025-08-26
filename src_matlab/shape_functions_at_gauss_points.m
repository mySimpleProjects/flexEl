function [Nw, Ntx, Nty, dNw_dxi, dNw_deta] = shape_functions_at_gauss_points(xi, eta)
    % Calcola le funzioni di forma dell'elemento Kirchhoff e le derivate
    % parametriche di Nw nel singolo punto di Gauss [xi, eta].

    % Nodi di riferimento (locali)
    coords = [-1, -1;
               1, -1;
               1,  1;
              -1,  1];

    % Matrice C per isolamento dei dof
    C = c_matrix(coords);

    idx_disp = [1, 4, 7, 10];
    idx_rotx = [2, 5, 8, 11];
    idx_roty = [3, 6, 9, 12];

    % Preallocazioni
    Nw  = zeros(1,4);
    Ntx = zeros(1,4);
    Nty = zeros(1,4);
    dNw_dxi  = zeros(1,4);
    dNw_deta = zeros(1,4);

    % Costruisco gli alpha per ciascun dof del nodo k
    for k = 1:4
        e_w = zeros(12,1);   e_w(idx_disp(k)) = 1;
        e_tx = zeros(12,1);  e_tx(idx_rotx(k)) = 1;
        e_ty = zeros(12,1);  e_ty(idx_roty(k)) = 1;

        alpha_w  = C \ e_w;
        alpha_tx = C \ e_tx;
        alpha_ty = C \ e_ty;

        % Assemblaggio vettore polinomiale NN(row1→w, row2→tx, row3→ty)
        NN = [
            1, -xi, -eta, -xi^2/2, -xi^3/6, -(xi^2*eta/2 + eta^4/12), -(xi^3*eta/6), -eta^2/2, -(eta^2*xi/2 + xi^4/12), -eta^3/6, -(eta^3*xi/6), -xi*eta/2;
            0,   1,    0,  xi,       xi^2/2,         xi*eta,           (xi^2*eta)/2,        0,      eta^2/2 + xi^3/3,  0,       eta^3/6,    eta/2;
            0,   0,    1,   0,             0,    (eta^3/3 + xi^2/2), xi^3/6,             eta, xi*eta, eta^2/2, (eta^2*xi)/2, xi/2
        ];

        % calcolo valori di Nw, Ntx, Nty
        val_w  = NN * alpha_w;
        val_tx = NN * alpha_tx;
        val_ty = NN * alpha_ty;

        Nw(k)  = val_w(1);
        Ntx(k) = val_tx(1);
        Nty(k) = val_ty(1);

        % --- calcolo dNN/dxi e dNN/deta della prima riga (w) ---
        dNN_dxi_row1 = [
            0, -1,  0,     -xi,        -xi^2/2,    -(xi*eta),      -(xi^2*eta/2),    0,   -(eta^2/2 + xi^3/3), 0,   -(eta^3/6),    -eta/2
        ];
        dNN_deta_row1 = [
            0,  0, -1,      0,            0,    -(xi^2/2 + eta^3/3), -(xi^3/6),   -eta,    -(eta*xi),  -eta^2/2,   -(eta^2*xi/2),  -xi/2
        ];

        % derivate parametriche di Nw
        dNw_dxi(k)  = dNN_dxi_row1  * alpha_w;
        dNw_deta(k) = dNN_deta_row1 * alpha_w;
    end
end
