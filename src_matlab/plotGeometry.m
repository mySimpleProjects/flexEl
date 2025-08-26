function plotGeometry(coords, conn, BC)
% PLOTGEOMETRY: 3D mesh + nodi+etichette+BC+carichi

DEBUG = false;

if nargin<3, BC = []; end
figure; hold on;

% Mesh a z=0
verts3D = [coords, zeros(size(coords,1),1)];
patch('Vertices', verts3D, 'Faces', conn, ...
    'FaceColor',[0.8 0.9 1],'FaceAlpha',0.5,'EdgeColor','k');

if DEBUG
    % Nodi e label
    for i=1:size(coords,1)
        x=coords(i,1); y=coords(i,2); z=0;
        plot3(x,y,z,'ko','MarkerFaceColor','k');
        text(x,y,z,[' N' num2str(i)],'Color','k','FontSize',10, ...
            'HorizontalAlignment','left','VerticalAlignment','bottom');
    end

    % Etichette elementi
    for e=1:size(conn,1)
        cen = mean(coords(conn(e,:),:),1);
        text(cen(1),cen(2),0,['E' num2str(e)], ...
            'FontSize',12,'FontWeight','bold','Color','b', ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

% Dirichlet (cubetti rossi)
if isfield(BC,'dirichlet') && ~isempty(BC.dirichlet.dofs)
    for gd=BC.dirichlet.dofs
        nd = ceil(gd/3);
        pt = [coords(nd,:),0];
        scatter3(pt(1),pt(2),pt(3),100,'s','filled','r');
    end
end

% Carichi (frecce rosse)
if isfield(BC,'force') && any(BC.force~=0)
    aL = max(range(coords(:,1)),range(coords(:,2)))/10;
    for k=find(BC.force~=0)'
        if mod(k-1,3)+1==1  % W
            nd  = ceil(k/3);
            val = BC.force(k);
            st  = [coords(nd,:),0];
            quiver3(st(1),st(2),st(3),0,0,sign(val)*aL,'r','LineWidth',1.5,'AutoScale','off');
        end
    end
end

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Plate Geometry with BCs & Loads');
view(3); axis equal; grid on; hold off;
end
