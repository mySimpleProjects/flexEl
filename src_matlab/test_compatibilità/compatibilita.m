clear all
close all
clc

%% Parametri materiali (non usati qui direttamente)
thickness = 1;
E = 210000;  % MPa
nu = 0.3;

%% Coordinate locali nel dominio [-1, 1]
coords = [-1, -1;
           1, -1;
           1,  1;
          -1,  1];

% Calcolo matrice C
C = c_matrix(coords);

% Griglia su dominio [-1, 1]
xv = linspace(-1, 1, 50);
yv = linspace(-1, 1, 50);
[X, Y] = meshgrid(xv, yv);

% Preallocazione delle funzioni di forma
Z_disp = zeros([size(X), 4]);   % w
Z_rotx = zeros([size(X), 4]);   % theta_x
Z_roty = zeros([size(X), 4]);   % theta_y

% Indici associati ai 12 gdl
idx_disp = [1, 4, 7, 10];
idx_rotx = [2, 5, 8, 11];
idx_roty = [3, 6, 9, 12];

% Loop su ciascun nodo
for k = 1:4
    rhs = zeros(12,1); rhs(idx_disp(k)) = 1;     % Spostamenti
    alpha = C \ rhs;

    rhsx = zeros(12,1); rhsx(idx_rotx(k)) = 1;   % Rotazione x
    alphax = C \ rhsx;

    rhsy = zeros(12,1); rhsy(idx_roty(k)) = 1;   % Rotazione y
    alphay = C \ rhsy;

    for i = 1:numel(X)
        x = X(i); y = Y(i);

        NN = [
            1, -x, -y, -x^2/2, -x^3/6, -((x^2*y)/2 + y^4/12), -((x^3*y)/6), -y^2/2, -((y^2*x)/2 + x^4/12), -y^3/6, -((y^3*x)/6), -x*y/2;
            0,  1,  0, x, x^2/2, x*y, (x^2*y)/2, 0, y^2/2 + x^3/3, 0, y^3/6, y/2;
            0,  0,  1, 0, 0, y^3/3 + x^2/2, x^3/6, y, x*y, y^2/2, (y^2*x)/2, x/2
        ];

        val    = NN * alpha;
        valx   = NN * alphax;
        valy   = NN * alphay;

        Z_disp(i + (k-1)*numel(X)) = val(1);
        Z_rotx(i + (k-1)*numel(X)) = valx(1);
        Z_roty(i + (k-1)*numel(X)) = valy(1);
    end
end

%% Plot funzioni di forma - Spostamenti verticali (w)
figure('Name','Funzioni di forma - Spostamenti (w)');
for k = 1:4
    subplot(2,2,k);
    surf(X, Y, Z_disp(:,:,k));
    title(['Nodo ', num2str(k), ' - w']);
    xlabel('x'); ylabel('y'); zlabel('N_w');
    axis equal tight;
end

%% Plot funzioni di forma - Rotazioni attorno a x
figure('Name','Funzioni di forma - Rotazioni \theta_x');
for k = 1:4
    subplot(2,2,k);
    surf(X, Y, Z_rotx(:,:,k));
    title(['Nodo ', num2str(k), ' - \theta_x']);
    xlabel('x'); ylabel('y'); zlabel('N_{\theta x}');
    axis equal tight;
end

%% Plot funzioni di forma - Rotazioni attorno a y
figure('Name','Funzioni di forma - Rotazioni \theta_y');
for k = 1:4
    subplot(2,2,k);
    surf(X, Y, Z_roty(:,:,k));
    title(['Nodo ', num2str(k), ' - \theta_y']);
    xlabel('x'); ylabel('y'); zlabel('N_{\theta y}');
    axis equal tight;
end
