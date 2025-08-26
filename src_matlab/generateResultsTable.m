function T = generateResultsTable(data)
    % GENERATERESULTSTABLE: crea e stampa la tabella dei risultati
    % INPUT:
    %   data.coords        Nx2 array di nodi
    %   data.fieldFun      @(x,y)→[W_ex,bx_ex,by_ex]  (opzionale)
    %   U, reactions       vettori globali [nDOF×1]
    %
    % OUTPUT:
    %   T  = MATLAB table con colonne:
    %        Node, X, Y,
    %        W_num, betaX_num, betaY_num,
    %        R_W, R_betaX, R_betaY,
    %        (se fieldFun) W_ex, betaX_ex, betaY_ex,
    %        (se fieldFun) Err_W, Err_betaX, Err_betaY
    
    U = data.U;
    reactions = data.reactions;
    coords = data.coords;
    nNodes = size(coords,1);

    % Indici nodi e coordinate
    Node = (1:nNodes).';
    X    = coords(:,1);
    Y    = coords(:,2);

    % Estrai numerici
    Umat      = reshape(U,3,[])';  % [W,betaX,betaY] per riga
    W_num     = Umat(:,1);
    betaX_num = Umat(:,2);
    betaY_num = Umat(:,3);

    % Estrai reazioni
    Rmat     = reshape(reactions,3,[])';
    R_W      = Rmat(:,1);
    R_betaX  = Rmat(:,2);
    R_betaY  = Rmat(:,3);

    % Costruisci tabella base
    T = table(Node, X, Y, ...
              W_num, betaX_num, betaY_num, ...
              R_W, R_betaX, R_betaY);

    % Se è definita una funzione di campo esatto, aggiungo colonne
    if isfield(data,'fieldFun') && isa(data.fieldFun,'function_handle')
        W_ex     = zeros(nNodes,1);
        betaX_ex = zeros(nNodes,1);
        betaY_ex = zeros(nNodes,1);
        for i=1:nNodes
            [W_ex(i), betaX_ex(i), betaY_ex(i)] = data.fieldFun(X(i), Y(i));
        end
        Err_W     = abs(W_num     - W_ex);
        Err_betaX = abs(betaX_num - betaX_ex);
        Err_betaY = abs(betaY_num - betaY_ex);

        T.W_ex      = W_ex;
        T.betaX_ex  = betaX_ex;
        T.betaY_ex  = betaY_ex;
        T.Err_W     = Err_W;
        T.Err_betaX = Err_betaX;
        T.Err_betaY = Err_betaY;
    end

    % Stampa a schermo
    disp(T);
    fprintf("Total reaction Rz= = %.2f mm\n", sum(R_W));
    fprintf("Total reaction Mx= = %.2f mm\n", sum(R_betaX));
    fprintf("Total reaction My= = %.2f mm\n", sum(R_betaY));
end
