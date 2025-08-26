%% Q_matrix.m
function Q = Q_matrix(x, y)
    Q = [0,0,0,1,   x,   y,  x*y, 0,   x^2, 0,   0,  0;
         0,0,0,0,   0,   y^2,0,   1,   x,   y,   x*y,0;
         0,0,0,0,   0,   2*x,x^2,0,   2*y, 0,   y^2,1];
end