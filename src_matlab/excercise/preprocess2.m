function data = preprocess2()
    % PREPROCESS2: mesh 4×4 su [0,10]×[0,10], clamp bordo, carico centro,
    % genera conn oraria e poi la inverte CCW

    %% Parametri
    Lx = 10; Ly = 10;
    nx = 4;  ny = 4;

    %% Nodi
    dx = Lx/nx; dy = Ly/ny;
    [X,Y] = meshgrid(0:dx:Lx, 0:dy:Ly);
    data.coords = [X(:), Y(:)];

    %% Connectivity (inizialmente ORARIO)
    conn = zeros(nx*ny,4);
    e = 0;
    for j = 1:ny
        for i = 1:nx
            e = e + 1;
            n1 = (j-1)*(nx+1) + i;
            n2 = n1 + 1;
            n3 = n2 + (nx+1);
            n4 = n1 + (nx+1);
            conn(e,:) = [n1, n2, n3, n4];
        end
    end

    % **Inverti** ogni riga per passare a CCW
    conn = conn(:, [1 4 3 2]);
    data.conn = conn;

    %% Gauss
    nGP = 3;
    [data.xi, data.w] = lgwt(nGP,-1,1);

    %% Materiali
    data.thickness = 0.2;
    data.E         = 1e6;
    data.nu        = 0.3;

    %% BC init
    nNodes = size(data.coords,1);
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
    % helper Forze
    function addF(nd, type, val)
        types = {'W','betaX','betaY'};
        col   = find(strcmp(types,type));
        gd    = (nd-1)*3 + col;
        data.BC.force(gd) = data.BC.force(gd) + val;
    end

    %% Clamped bordo esterno
    bnodes = find( ...
      abs(data.coords(:,1))<1e-8 | abs(data.coords(:,1)-Lx)<1e-8 | ...
      abs(data.coords(:,2))<1e-8 | abs(data.coords(:,2)-Ly)<1e-8 );
    addD(bnodes,'W',0); addD(bnodes,'betaX',0); addD(bnodes,'betaY',0);

    %% Carico al centro
    center = [Lx/2, Ly/2];
    dists  = hypot(data.coords(:,1)-center(1), data.coords(:,2)-center(2));
    [~, cnode] = min(dists);
    addF(cnode,'W',-1);

    %% Visualizza
    plotGeometry(data.coords, data.conn, data.BC);
    title('PREPROCESS2: connettività invertita CCW');
    view(3); axis equal; grid on;
end
