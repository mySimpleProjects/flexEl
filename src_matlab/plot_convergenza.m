function plot_convergenza(h, w_num, w_ref, ordine_atteso)
% plot_convergenza - Crea un grafico di convergenza bilogaritmico
%
% Sintassi:
%   plot_convergenza(h, w_num, w_ref, ordine_atteso)
%
% Input:
%   h              - Vettore delle lunghezze caratteristiche (es. lato medio elemento)
%   w_num          - Vettore dei valori numerici corrispondenti
%   w_ref          - Valore analitico di riferimento
%   ordine_atteso  - (facoltativo) Ordine atteso di convergenza (es. 2)
%
% Output:
%   Un grafico log-log dell'errore relativo rispetto alla lunghezza caratteristica

    % Calcola errore relativo
    errore = abs((w_num - w_ref) / w_ref);

    % Crea grafico
    figure;
    loglog(h, errore, '-o', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel('Characteristic length h [m]');
    ylabel('Relative error');
    title('Convergence plot');
    hold on;

    % Aggiungi retta di riferimento se richiesto
    if nargin == 4
        ref_h = [min(h), max(h)];
        % Calcola linea di riferimento con stesso ordine
        ref_err = errore(1) * (ref_h / h(1)).^ordine_atteso;
        loglog(ref_h, ref_err, '--k', 'LineWidth', 1.5);
        legend('Numerical error', sprintf('Expected order = %d', ordine_atteso), 'Location', 'southwest');
    else
        legend('Numerical error', 'Location', 'southwest');
    end

    hold off;
end
