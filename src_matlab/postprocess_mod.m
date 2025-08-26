function postprocess_mod(data, DEBUG, truscale, res)
% POSTPROCESS_MOD: draw each element with true isoparametric mapping
%   - bordo nero via shape‐functions lungo la frontiera
%   - interno colore interpolato da Wdisp
%
%   data     struct con .coords (Nx2), .conn (Mx4), .U (3N×1)
%   DEBUG    boolean: se true disegna la figura
%   truscale boolean: se true usa scale=1, altrimenti adatta
%   res      campionamenti per lato dell’elemento (default 21)

    if nargin<4, res = 21; end
    if ~DEBUG, return; end

    coords = data.coords;      % [Nnodes × 2]
    conn   = data.conn;        % [Nelements × 4]
    U      = data.U;           % [3*Nnodes × 1]

    % estraggo spostamento w
    dispAll = reshape(U,3,[])'; % [Nnodes × 3]
    Wdisp   = dispAll(:,1);

    % scala grafico
    dims = [range(coords(:,1)), range(coords(:,2))];
    if truscale
        scale = 1;
    elseif max(abs(Wdisp))>0
        scale = 0.1 * max(dims) / max(abs(Wdisp));
    else
        scale = 1;
    end

    figure; hold on;
    colormap(flipud(jet));
    caxis([min(Wdisp), max(Wdisp)]);

    % parametri isoparametrici
    xi_vals  = linspace(-1,1,res);
    eta_vals = linspace(-1,1,res);
    [Xi, Eta] = meshgrid(xi_vals, eta_vals);

    for e = 1:size(conn,1)
        nodes = conn(e,:);
        Xe = coords(nodes,1);
        Ye = coords(nodes,2);
        We = Wdisp(nodes);

        % 1) superficie interna colorata
        Xg = zeros(res,res);
        Yg = zeros(res,res);
        Zg = zeros(res,res);
        Cg = zeros(res,res);

        for i = 1:res
            for j = 1:res
                xi  = Xi(i,j);
                eta = Eta(i,j);
                Nw  = shape_functions_at_gauss_points(xi, eta);  % [1×4]

                % mapping isoparametrico
                Xg(i,j) = Nw * Xe;
                Yg(i,j) = Nw * Ye;

                % spostamento scalato
                wloc    = Nw * We;
                Zg(i,j) = wloc * scale;
                Cg(i,j) = wloc;
            end
        end

        % disegno l’interno senza linee
        surf(Xg, Yg, Zg, Cg, 'EdgeColor','none');

        % 2) bordo nero tramite campionamento shape‐functions
        %    (4 lati: xi=±1, eta=±1)
        %    creiamo array di punti per ciascun lato e li plottiamo
        %    come curve continue
        % lato xi = -1
        eta_b = eta_vals;
        xi_b  = -1 * ones(1,res);
        plot_iso_edge(xi_b, eta_b, Xe, Ye, We, scale);

        % lato xi = +1
        xi_b  = +1 * ones(1,res);
        plot_iso_edge(xi_b, eta_b, Xe, Ye, We, scale);

        % lato eta = -1
        xi_b  = xi_vals;
        eta_b = -1 * ones(1,res);
        plot_iso_edge(xi_b, eta_b, Xe, Ye, We, scale);

        % lato eta = +1
        eta_b = +1 * ones(1,res);
        plot_iso_edge(xi_b, eta_b, Xe, Ye, We, scale);
    end

    colorbar;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('POSTPROCESS\_MOD: Deformed (scale×%.2f)', scale));
    view(3); axis equal; grid on;
    hold off;
end

function plot_iso_edge(xi_b, eta_b, Xe, Ye, We, scale)
    % helper: campiona la frontiera (xi_b,eta_b) e ne traccia il bordo
    % Xe,Ye: vettori 4×1 di coordinate nodali
    % We: spostamenti nodali
    nb = numel(xi_b);
    Xb = zeros(1,nb);
    Yb = zeros(1,nb);
    Zb = zeros(1,nb);
    for k = 1:nb
        Nw = shape_functions_at_gauss_points(xi_b(k), eta_b(k));
        Xb(k) = Nw * Xe;
        Yb(k) = Nw * Ye;
        Zb(k) = (Nw * We) * scale;
    end
    plot3(Xb, Yb, Zb, 'k-', 'LineWidth',1.5);
end
