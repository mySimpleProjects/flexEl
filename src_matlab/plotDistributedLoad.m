function plotDistributedLoad(coords, conn, q)
    % PLOTDISTRIBUTEDLOAD: disegna frecce più piccole che arrivano col vertice
    % al baricentro di ogni elemento (punta sulla superficie z=0)
    % coords : Nx2 array di nodi [x,y]
    % conn   : Mx4 array di elementi (CCW)
    % q      : carico superficiale uniforme (N/m^2)

    % 1) Calcola baricentri (z = 0)
    nElem = size(conn,1);
    cents = zeros(nElem,3);
    for e = 1:nElem
        c = mean(coords(conn(e,:),:),1);
        cents(e,:) = [c, 0];
    end

    % 2) Arrow length ridotta (5% della dimensione massima)
    span    = max(range(coords(:,1)), range(coords(:,2)));
    arrowLen = span * 0.05;

    % 3) Punto di partenza e direzione
    sgn = sign(q);
    zStart = -sgn * arrowLen;        % inverte: punta sempre verso z=0
    dZ     =  sgn * arrowLen;        

    X0 = cents(:,1);
    Y0 = cents(:,2);
    Z0 = repmat(zStart, nElem, 1);
    dX = zeros(nElem,1);
    dY = zeros(nElem,1);
    dZ = repmat(dZ, nElem, 1);

    % 4) Disegna le frecce rosse con testa in z=0
    quiver3( ...
      X0, Y0, Z0, ...   % punti di partenza sopra/sotto l’elemento
      dX, dY, dZ, ...   % vettore freccia
      'r', 'LineWidth',1.2, 'MaxHeadSize',0.3, 'AutoScale','off' );
end
