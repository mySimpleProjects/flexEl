function data = lastraQuadrataSS(data, DEBUG)
% LASTRAQUADRATASS: genera mesh ne×ne su [0,L]²,
% piastra semplicemente appoggiata (solo W=0 al bordo),
% e definisce un carico distribuito uniforme q.
%
% INPUT:
%   data.L  = lato della lastra
%   data.ne = numero di elementi per lato (pari)
%
% OUTPUT in data:
%   coords, conn, xi, w, thickness, E, nu,
%   BC.dirichlet (solo W sul bordo),
%   data.q = valore della pressione uniforme (N/m²)

fprintf('LASTRA SEMPLICEMENTE APPOGGIATA CON CARICO DISTRIBUITO\n\n');

% INPUT DATI:
% Carico distribuito uniforme (negativo = verso il basso)
data.q = -0.1;  % [N/m^2]
data.F = 0;
data.L = 100;

% per semplicità del codice
F = data.F;
L = data.L;
ne = data.ne;


% Gauss–Legendre integrazione 3×3
nGP = 3;
[data.xi, data.w] = lgwt(nGP, -1, 1);

% Nodi
dx = L/ne; dy = L/ne;
[X,Y] = meshgrid(0:dx:L, 0:dy:L);
data.coords = [X(:), Y(:)];

% Connectivity CCW
nnx = ne+1;
conn = zeros(ne*ne,4);
e = 0;
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
% Inverti per CCW (anti-orario)
data.conn = conn(:,[1 4 3 2]);

% BC init
nNodes = size(data.coords,1);
data.BC.dirichlet.dofs   = [];
data.BC.dirichlet.values = [];
data.BC.force            = zeros(nNodes*3,1);

% helper Dirichlet (solo W)
    function addD(nds)
        for nd = nds(:)'
            gd = (nd-1)*3 + 1;      % W è il primo DOF
            data.BC.dirichlet.dofs(end+1)   = gd;
            data.BC.dirichlet.values(end+1) = 0;
        end
    end

% Trova i nodi di bordo e fissa W=0
tol = 1e-8;
border = find( ...
    abs(data.coords(:,1))<tol    | abs(data.coords(:,1)-L)<tol | ...
    abs(data.coords(:,2))<tol    | abs(data.coords(:,2)-L)<tol );
addD(border);

% Nota: le forze distribuite (data.q) verranno integrate nel solver

% Visualizza mesh e vincoli
if DEBUG
    plotGeometry(data.coords, data.conn, data.BC);
    hold on;
    plotDistributedLoad(data.coords, data.conn, data.q)
    hold off;

    title(sprintf('PREPROCESS SS: %dx%d Mesh, L=%.2f, q=%.1f N/m^2', ne, ne, L, data.q));
    view(3); axis equal; grid on;
end
end
