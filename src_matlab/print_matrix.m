%% print_matrix.m
function print_matrix(M, name, digits)
    if nargin<3, digits=3; end
    if nargin<2, name='Matrix'; end
    [r,c] = size(M);
    fieldW = digits + 7;
    fmt = ['%', num2str(fieldW), '.', num2str(digits), 'f '];
    fprintf('\n%s (%dx%d):\n', name, r, c);
    for i = 1:r
        for j = 1:c
            fprintf(fmt, M(i,j));
        end
        fprintf('\n');
    end
end