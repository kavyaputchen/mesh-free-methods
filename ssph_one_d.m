n = 100;
x0 = 0;
xn = 1;
lambda = 1;
G = 1.04823;
L = (xn - x0)/n;
x = [0:L:1]';
h = 1.2*(L);
f = sin(x);
P = zeros(2,n);
Q = zeros(2,n);
dx = ones(n+1,1)*L;
dx(1) = 0.5*L; dx(end) = 0.5*L;
% y = [];


for i = 1:n+1
    % For the neighbours of each point in sample space
    d_ij = [];
    z_i = [];
    w_ij = [];
    P_ijT = [];
    Q_ij = [];
    dx(i) = L;
    R = zeros(2,1);
    A = zeros(2,2);
    nei = [];
    
    
    % for all neighbours of a point, calc d_ij, z_i
    for j = 1:n+1
        if i ~= j
            dij = norm(x(j) - x(i));
            xi = x(j)-x(i);
            if dij <= 2*h
                nei = [nei,j];
                P_ijT = [xi^2 xi^3;xi^3 xi^4];
                w_ij = (G/((h*sqrt(pi))^lambda))*(exp(-(dij^2)/h^2)-exp(-4));
                A = A + P_ijT*w_ij*dx(j);
                R = R + (f(j)-f(i))*[xi;xi^2]*w_ij*dx(j);
            end   
        end
    end
    Q_ij = A\R;
    Q_ijT = Q_ij';
    dfdx(i) = Q_ij(1);
    d2fdx2(i) = Q_ij(2)*2;
    
    end