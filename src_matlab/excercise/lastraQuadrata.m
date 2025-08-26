function data = lastraQuadrata(data)
% PREPROCESS3: mesh ne×ne su [0,L]×[0,L], clamped al bordo, carico al centro
fprintf("LASTRA QUADRATA VINCOLATA AI BORDI CON CARICO CONCENTRATO\n\n")

% INPUT DATI:
data.F = -1000;
data.L = 100;
data.ne = 20;
data.q = 0;

% per semplicità del codice
F = data.F;
L = data.L;
ne = data.ne;

%% Nodi
dx = L/ne; dy = L/ne;
[X,Y] = meshgrid(0:dx:L, 0:dy:L);
data.coords = [X(:), Y(:)];

%% Connectivity CCW
conn = zeros(ne*ne,4);
e = 0;
nnx = ne+1;
for j = 1:ne
    for i = 1:ne
        e = e + 1;
        n1 = (j-1)*nnx + i;      % bottom-left
        n2 = n1 + 1;             % bottom-right
        n3 = n2 + nnx;           % top-right
        n4 = n1 + nnx;           % top-left
        conn(e,:) = [n1 n2 n3 n4];
    end
end

% **Inverti** ogni riga per passare a CCW
conn = conn(:, [1 4 3 2]);
data.conn = conn;

%% Gauss–Legendre integrazione
nGP = 3;
[data.xi, data.w] = lgwt(nGP, -1, 1);

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

%% Clamped su bordo esterno
bnodes = find( ...
    abs(data.coords(:,1))<1e-8 | abs(data.coords(:,1)-L)<1e-8 | ...
    abs(data.coords(:,2))<1e-8 | abs(data.coords(:,2)-L)<1e-8 );
addD(bnodes,'W',0);
addD(bnodes,'betaX',0);
addD(bnodes,'betaY',0);

%% Carico puntuale al centro
center = [L/2, L/2];
dists  = hypot(data.coords(:,1)-center(1), data.coords(:,2)-center(2));
[~, cnode] = min(dists);
addF(cnode,'W',F);

%% Visualizza
plotGeometry(data.coords, data.conn, data.BC);
title(sprintf('PREPROCESS3: %dx%d Mesh, L=%.2f, Clamped, Center Load', ne, ne, L));
view(3); axis equal; grid on;
end
