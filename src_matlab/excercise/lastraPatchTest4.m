function data = lastraPatchTest4(data, DEBUG)
% LASTRAPATCHTEST4: patch-test KSBQP con polinomio di 4° ordine
%   impone W, betaX, betaY ai nodi di bordo secondo
%   w(x,y)=a0 + a1 x + a2 y + a3 x^2 + a4 x y + a5 y^2
%                + a6 x^3 + a7 x^2 y + a8 x y^2 + a9 y^3
%                + a10 x^4 + a11 x^3 y + a12 x^2 y^2 + a13 x y^3 + a14 y^4
%
% INPUT:
%   data.ne   = n. elementi per lato (pari)
%   DEBUG     = true|false
%
% OUTPUT in data:
%   data.L, data.poly, data.q,
%   coords, conn, xi, w, thickness, E, nu,
%   BC.dirichlet.dofs/.values, BC.force

fprintf('PATCH-TEST grado 4: impongo W, βx, βy sui bordi\n');

% dati di esempio
data.L    = 100;
data.poly = [ ...
    0.001,  0.01,  -0.02,   0.005,  0.003, -0.001, ...
    0.0008, 0.0006, -0.0004, 0.0002, -0.0001, ...
    0.00015, -0.00015, 0.00005, -0.00005]./100;
data.q = 0;  % nessun carico

L  = data.L;
ne = data.ne;
a  = data.poly;

% materiale
data.thickness = 0.2;
data.E         = 1e6;
data.nu        = 0.3;

% Gauss 3×3
nGP = 3;
[data.xi, data.w] = lgwt(nGP, -1, 1);

% nodi
dx = L/ne; dy = L/ne;
[X,Y] = meshgrid(0:dx:L, 0:dy:L);
data.coords = [X(:), Y(:)];

% connessione CCW
nnx = ne+1;
conn = zeros(ne*ne,4);
e = 0;
for j = 1:ne
    for i = 1:ne
        e = e + 1;
        n1 = (j-1)*nnx + i;
        n2 = n1 + 1;
        n3 = n2 + nnx;
        n4 = n1 + nnx;
        conn(e,:) = [n1 n2 n3 n4];
    end
end
data.conn = conn(:, [1 4 3 2]);

% inizializza BC
nNodes = size(data.coords,1);
data.BC.dirichlet.dofs   = [];
data.BC.dirichlet.values = [];
data.BC.force            = zeros(nNodes*3,1);

% helper per Dirichlet
    function addD(nds)
        for nd = nds(:)'
            x = data.coords(nd,1);
            y = data.coords(nd,2);

            % w di 4° ordine
            wv = a(1) ...
                + a(2)*x + a(3)*y ...
                + a(4)*x^2 + a(5)*x*y + a(6)*y^2 ...
                + a(7)*x^3 + a(8)*x^2*y + a(9)*x*y^2 + a(10)*y^3 ...
                + a(11)*x^4 + a(12)*x^3*y + a(13)*x^2*y^2 + a(14)*x*y^3 + a(15)*y^4;

            % derivate
            dwdx = a(2) ...
                + 2*a(4)*x   + a(5)*y ...
                + 3*a(7)*x^2 + 2*a(8)*x*y + a(9)*y^2 ...
                + 4*a(11)*x^3 + 3*a(12)*x^2*y + 2*a(13)*x*y^2 + a(14)*y^3;

            dwdy = a(3) ...
                + a(5)*x     + 2*a(6)*y ...
                + a(8)*x^2   + 2*a(9)*x*y + 3*a(10)*y^2 ...
                + a(12)*x^3  + 2*a(13)*x^2*y + 3*a(14)*x*y^2 + 4*a(15)*y^3;

            bx = -dwdx;
            by = -dwdy;

            % impongo W, βx, βy
            data.BC.dirichlet.dofs(end+1)   = (nd-1)*3 + 1;
            data.BC.dirichlet.values(end+1) = wv;
            data.BC.dirichlet.dofs(end+1)   = (nd-1)*3 + 2;
            data.BC.dirichlet.values(end+1) = bx;
            data.BC.dirichlet.dofs(end+1)   = (nd-1)*3 + 3;
            data.BC.dirichlet.values(end+1) = by;
        end
    end

% applica Dirichlet su bordo
tol = 1e-8;
border = find( ...
    abs(data.coords(:,1))<tol    | abs(data.coords(:,1)-L)<tol | ...
    abs(data.coords(:,2))<tol    | abs(data.coords(:,2)-L)<tol );
addD(border);

% plot di debug
if nargin>1 && DEBUG
    plotGeometry(data.coords, data.conn, data.BC);
    title('PATCH-TEST grado 4 – BC imposte');
    axis equal; view(3); grid on;
end
end
