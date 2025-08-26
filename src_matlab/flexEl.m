%% flexEl.m

close all; clearvars; clc;
format shortG;
warning('off');
PLOT = true;

%% Materiali (default)
data.thickness = 1;
data.E         = 210000;  % MPa
data.nu        = 0.3;

%% Parametri di simulazione
nelem = [5];            % numero di elementi (es. per lato)
w_num = zeros(size(nelem));     % preallocazione

%% Loop sulle mesh
for i = 1:length(nelem)
    fprintf("\n--- Simulazione con mesh %dx%d ---\n", nelem(i), nelem(i));
    
    data.ne = nelem(i);  % numero di elementi per lato
    % data = lastraQuadrataSS(data,PLOT);  % genera il problema
    % data = lastraPatchTest0(data, PLOT);
    % data = lastraPatchTest1(data, PLOT);
    % data = lastraPatchTest2(data, PLOT);
    data = lastraPatchTest3(data, PLOT);
    
    % Solver FEM
    [data.U, data.reactions] = solver(data);
    
    % Postprocess (opzionale)
    % postprocess(data,PLOT,true);
    postprocess_mod(data,PLOT,true,10);

    T = generateResultsTable(data);

    % Estrazione massimo spostamento
    if abs(min(data.U)) > abs(max(data.U))
        wn = min(data.U);
    else
        wn = max(data.U);
    end

    w_num(i) = wn;  % salva il risultato per la convergenza
    h(i) = sqrt(2)*data.L/data.ne;

    fprintf("Wn: %.6f mm\n", wn);
end

%% Calcolo una sola volta la soluzione analitica
wa = lastraQuadrataAnalytica(data);
fprintf("\nAnalytical Solution: Wa = %.6f mm\n", wa);

%% Plot convergenza
plot_convergenza(h, w_num, wa, 2);
