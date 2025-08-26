function data = preprocess()
    % PREPROCESS: definisce geometria, materiale, integrazione e BC

    % Geometry: 6 nodes, 2 elements
    data.coords = [ 0,  0;
                    10, 0;
                    20, 0;
                     0,10;
                    10,10;
                    20,10];
    data.conn   = [1,2,5,4;
                   2,3,6,5];

    % Material & thickness
    data.thickness = 0.2;
    data.E         = 1e6;
    data.nu        = 0.3;

    % Gauss–Legendre integration (3×3)
    nGP = 3;
    [data.xi, data.w] = lgwt(nGP, -1, 1);

    % Initialize BC structures
    nNodes = size(data.coords,1);
    nDOF   = nNodes * 3;             % [W, betaX, betaY] per node
    data.BC.dirichlet.dofs   = [];
    data.BC.dirichlet.values = [];
    data.BC.force            = zeros(nDOF,1);

    % Helper: add Dirichlet BC
    function addD(node, type, val)
        types = {'W','betaX','betaY'};
        col   = find(strcmp(types,type));
        for nd = node
            gdof = (nd-1)*3 + col;
            data.BC.dirichlet.dofs(end+1)   = gdof;
            data.BC.dirichlet.values(end+1) = val;
        end
    end

    % Helper: add point force
    function addF(node, type, val)
        types = {'W','betaX','betaY'};
        col   = find(strcmp(types,type));
        gdof  = (node-1)*3 + col;
        data.BC.force(gdof) = data.BC.force(gdof) + val;
    end

    % Define BCs (example)
    % Fix all DOFs at nodes 1 and 4
    addD(1, 'W',     0); addD(1, 'betaX', 0); addD(1, 'betaY', 0);
    addD(4, 'W',     0); addD(4, 'betaX', 0); addD(4, 'betaY', 0);
    % Apply vertical point loads at nodes 2 and 5
    addF(2, 'W', -1);
    addF(5, 'W', -1);

    % Plot 3D geometry with BC symbols & loads
    plotGeometry(data.coords, data.conn, data.BC);
    title('PREPROCESS: Geometry & BCs');
    view(3); axis equal; grid on;
end
