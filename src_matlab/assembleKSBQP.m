function [Kglob, Fglob, Ae] = assembleKSBQP(coords, conn, thickness, E, nu, xi, w, q)
    % ASSEMBLEKSBQP: assembla matrice di rigidezza e vettore forze per KSBQP
    % coords, conn       : geometria della mesh
    % thickness, E, nu   : proprietà del materiale
    % xi, w              : punti e pesi di Gauss
    % q                  : carico distribuito costante (N/m^2)
    %
    % OUTPUT:
    %   Kglob : matrice globale (nDOF×nDOF)
    %   Fglob : vettore globale delle forze nodali coerenti (nDOF×1)
    %   Ae    : vettore aree elementi (nElem×1)

    nNodes = size(coords,1);
    nElem  = size(conn,1);
    nDOF   = nNodes*3;

    Kglob  = zeros(nDOF);
    Fglob  = zeros(nDOF,1);
    Ae     = zeros(nElem,1);

    for e = 1:nElem
        % Estrai nodi e coordinate dell'elemento
        nodes      = conn(e,:);
        elemCoords = coords(nodes,:);

        % Calcola Ke e area elemento
        % [Ke, Ae(e)] = ksbqp_element_mod(elemCoords, thickness, E, nu, xi, w);
        [Ke, Ae(e)] = ksbqp_element(elemCoords, thickness, E, nu, xi, w);

        % Assembla rigidezza
        dofMap = reshape([(nodes-1)*3+1; (nodes-1)*3+2; (nodes-1)*3+3], 1, []);
        Kglob(dofMap, dofMap) = Kglob(dofMap, dofMap) + Ke;

        % Calcola forze nodali consistenti dovute a q
        Fe = assembleForceElement(elemCoords, xi, w, q);

        % Assembla forze
        Fglob(dofMap) = Fglob(dofMap) + Fe;
    end
end
