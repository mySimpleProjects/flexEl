function data = preprocessCircle(F, R, n)
    % PREPROCESSCIRCLE: mesh quadrangolare strutturata su disco di raggio R,
    % bordo esterno incastrato, carico puntuale al centro.
    % Input:
    %   R : raggio del cerchio
    %   n : numero di suddivisioni per lato della griglia (mesh (n+1)×(n+1))

    %% 1) Genera griglia strutturata su quadrato [-R,R]×[-R,R]
    xi = linspace(-R, R, n+1);
    yi = linspace(-R, R, n+1);
    [X, Y] = meshgrid(xi, yi);
    coords = [X(:), Y(:)];

    %% 2) Connnectivity quadrangolare standard
    conn = zeros(n*n, 4);
    e = 0;
    NN = n+1;
    for j = 1:n
        for i = 1:n
            e = e + 1;
            n1 = (j-1)*NN + i;
            n2 = n1 + 1;
            n3 = n2 + NN;
            n4 = n1 + NN;
            conn(e,:) = [n1, n2, n3, n4];
        end
    end
    % **Inverti** ogni riga per passare a CCW
    conn = conn(:, [1 4 3 2]);
    data.conn = conn;

    %% 3) Proietta i nodi esterni sul cerchio
    r = sqrt(coords(:,1).^2 + coords(:,2).^2);
    outside = r > R;
    coords(outside,1) = coords(outside,1)./r(outside) * R;
    coords(outside,2) = coords(outside,2)./r(outside) * R;

    %% 4) Salva in data
    data.coords    = coords;
    data.conn      = conn;
    data.thickness = 0.2;
    data.E         = 1e6;
    data.nu        = 0.3;

    %% 5) Gauss–Legendre integrazione
    nGP = 3;
    [data.xi, data.w] = lgwt(nGP,-1,1);

    %% 6) Boundary Conditions
    nNodes = size(coords,1);
    data.BC.dirichlet.dofs   = [];
    data.BC.dirichlet.values = [];
    data.BC.force            = zeros(nNodes*3,1);

    % helper Dirichlet
    function addD(nds, type, val)
        types = {'W','betaX','betaY'};
        col   = find(strcmp(types,type));
        for nd = nds(:)'
            gd = (nd-1)*3 + col;
            data.BC.dirichlet.dofs(end+1)   = gd;
            data.BC.dirichlet.values(end+1) = val;
        end
    end
    % helper forces
    function addF(nd, type, val)
        types = {'W','betaX','betaY'};
        col   = find(strcmp(types,type));
        gd    = (nd-1)*3 + col;
        data.BC.force(gd) = data.BC.force(gd) + val;
    end

    % a) incastra tutti i nodi sul bordo (r≈R)
    tol = 1e-6;
    onBoundary = abs(sqrt(data.coords(:,1).^2 + data.coords(:,2).^2) - R) < tol;
    bnodes = find(onBoundary);
    addD(bnodes,'W',0);
    addD(bnodes,'betaX',0);
    addD(bnodes,'betaY',0);

    % b) applica carico puntuale al centro (0,0)
    dists = sqrt(data.coords(:,1).^2 + data.coords(:,2).^2);
    [~, cnode] = min(dists);
    addF(cnode,'W',F);

    %% 7) Visualizza
    plotGeometry(data.coords, data.conn, data.BC);
    title(sprintf('PREPROCESSCIRCLE: R=%.2f, mesh=%dx%d', R, n, n));
    view(3); axis equal; grid on;
end
