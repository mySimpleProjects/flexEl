function postprocess(data, DEBUG, truscale)
% POSTPROCESS: plot underformed (dashed), deformed colored,
% outline + nodes, con segno coerente W->Z.

coords = data.coords;
conn   = data.conn;
U      = data.U;

% Estraggo W = prima colonna
dispAll = reshape(U,3,[])';  % [W, betaX, betaY]
Wdisp   = dispAll(:,1);

% Scala adattiva
dims    = [range(coords(:,1)), range(coords(:,2))];

if truscale
    scale = 1;
else
    scale   = 1;
    if max(abs(Wdisp))>0
        scale = 0.1 * max(dims) / max(abs(Wdisp));
    end
end

% Vertici deformati
defVerts = [coords, Wdisp * scale];

if DEBUG
    figure; hold on;
    % Originale tratteggiata
    verts0 = [coords, zeros(size(coords,1),1)];
    for e=1:size(conn,1)
        patch('Vertices', verts0, 'Faces', conn(e,:), ...
            'FaceColor','none','EdgeColor','k','LineStyle','--');
    end

    % Deformata colorata
    patch('Vertices', defVerts, 'Faces', conn, ...
        'FaceVertexCData', Wdisp, 'FaceColor','interp','EdgeColor','none');
    colormap(flipud(jet)); colorbar;
    caxis([min(Wdisp), max(Wdisp)]);

    % Outline + nodi
    patch('Vertices', defVerts, 'Faces', conn, ...
        'FaceColor','none','EdgeColor','k','LineWidth',1);
    plot3(defVerts(:,1), defVerts(:,2), defVerts(:,3), ...
        'ko','MarkerFaceColor','k','MarkerSize',5);

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('POSTPROCESS: Deformed Shape (Scale ×%.2f)',scale));
    view(3); axis equal; grid on; hold off;
    end
end
