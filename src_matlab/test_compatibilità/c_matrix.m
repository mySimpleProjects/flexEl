%% c_matrix.m
function C = c_matrix(coords)
    % Calcola la matrice C per un elemento a 4 nodi (Kirchhoff plate element)
    % Input:
    %   coords: 4x2 matrice delle coordinate dei nodi [x1 y1; x2 y2; x3 y3; x4 y4]
    % Output:
    %   C: matrice 12x12 di interpolazione per i gradi di libertà

    C = zeros(12);  % Preallocazione
    for i = 1:4
        Pi = P_matrix(coords(i,1), coords(i,2));  % Matrice P associata al nodo i
        C((i-1)*3 + (1:3), :) = Pi;
    end
end
