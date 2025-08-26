function data = lastraPatchTest3(data, DEBUG)
    % LASTRAPATCHTEST: patch-test per KSBQP su lastra [0,L]²
    %   impone W, betaX, betaY ai nodi di bordo secondo
    %   w(x,y)=a0+a1 x + a2 y + a3 x^2 + a4 x y + a5 y^2
    %
    % INPUT:
    %   data.L    = lato della lastra
    %   data.ne   = n. elementi per lato (pari)
    %   data.poly = [a0,a1,a2,a3,a4,a5]
    %   DEBUG (bool): se true fa il plotGeometry
    %
    % OUTPUT (in data):
    %   coords, conn, xi, w, thickness, E, nu,
    %   BC.dirichlet.dofs/.values (W,betaX,betaY sui bordi),
    %   BC.force = zeros (nessuna forza esterna)

    fprintf('PATCH-TEST: enforce boundary displacements and rotations\n');
    
    data.L = 100;
    data.poly= [0.01, 0.01, -0.02, 0.005, 0.003, -0.001, 0.001, 0.002, 0.003, -0.001]/150;
    % Carico distribuito uniforme (negativo = verso il basso)
    data.q = 0;  % [N/m^2]

    L  = data.L;
    ne = data.ne;

    a = data.poly;  % [a0,a1,a2,a3,a4,a5]

    % Materiale
    data.thickness = 0.2;
    data.E         = 1e6;
    data.nu        = 0.3;

    % Gauss 3×3
    nGP = 3;  
    [data.xi, data.w] = lgwt(nGP,-1,1);

    % 1) Nodi
    dx = L/ne; dy = L/ne;
    [X,Y] = meshgrid(0:dx:L, 0:dy:L);
    data.coords = [X(:), Y(:)];

    % 2) Connectivity CCW
    nnx = ne+1;
    conn = zeros(ne*ne,4);
    e = 0;
    for j = 1:ne
      for i = 1:ne
        e = e+1;
        n1 = (j-1)*nnx + i;
        n2 = n1 + 1;
        n3 = n2 + nnx;
        n4 = n1 + nnx;
        conn(e,:) = [n1,n2,n3,n4];
      end
    end
    data.conn = conn(:,[1 4 3 2]);  % garantisco CCW

    % 3) Inizializza BC
    nNodes = size(data.coords,1);
    data.BC.dirichlet.dofs   = [];
    data.BC.dirichlet.values = [];
    data.BC.force            = zeros(nNodes*3,1);

    % helper per Dirichlet
    function addD(nds, wv, bx, by)
      types = {'W','betaX','betaY'};
      for nd = nds(:)'
        x = data.coords(nd,1);
        y = data.coords(nd,2);
        % calcolo i valori
        valW  = a(1) + a(2)*x + a(3)*y + a(4)*x^2 + a(5)*x*y + a(6)*y^2 + a(7)*x^3 + a(8)*x^2*y + a(9)*x*y^2 + a(10)*y^3;
        valBx = -(          a(2)   + 2*a(4)*x   + a(5)*y + 3*a(7)*x^2 + 2*a(8)*x*y + a(9)*y^2);   % dw/dx
        valBy = -(          a(3)   + a(5)*x     + 2*a(6)*y + a(8)*x^2 + 2*a(9)*x*y + 3*a(10)*y^2); % dw/dy

        % se voglio specificare che fisso W e rotazioni
        if wv
          gd = (nd-1)*3 + 1; data.BC.dirichlet.dofs(end+1)=gd; data.BC.dirichlet.values(end+1)=valW;
        end
        if bx
          gd = (nd-1)*3 + 2; data.BC.dirichlet.dofs(end+1)=gd; data.BC.dirichlet.values(end+1)=valBx;
        end
        if by
          gd = (nd-1)*3 + 3; data.BC.dirichlet.dofs(end+1)=gd; data.BC.dirichlet.values(end+1)=valBy;
        end
      end
    end

    % 4) Applico Dirichlet su bordo (solo W, betaX, betaY)
    tol = 1e-8;
    border = find( ...
      abs(data.coords(:,1))<tol    | abs(data.coords(:,1)-L)<tol | ...
      abs(data.coords(:,2))<tol    | abs(data.coords(:,2)-L)<tol );
    addD(border, true, true, true);

    % 5) Visualizzo se DEBUG
    if nargin>1 && DEBUG
      plotGeometry(data.coords, data.conn, data.BC);
      title('PATCH-TEST: Spostamenti e rotazioni al bordo');
      view(3); axis equal; grid on;
    end
end
