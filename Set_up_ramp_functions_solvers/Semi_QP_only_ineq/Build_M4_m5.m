function [M4, m5,Qinit,actset] = Build_M4_m5(Ha_inv, Hb, Hess, c, ha, hb)
    na = size(Ha_inv, 1);  % number of active constraints
    nb = size(Hb, 1);      % number of inequality constraints
    nc = na+nb;

    % Projection matrices
    M_H = Ha_inv;
    M_H2 = M_H;
    MT_H2 = M_H2';

    % m5 construction
    m5_upper = MT_H2 * (Hess * M_H * ha + c);
    m5_lower = Hb * M_H * ha - hb;
    m5 = sparse([m5_upper; m5_lower]);

    % M4 construction
    M4_upper_left = speye(na) - MT_H2 * Hess * M_H2;
    M4_upper_right = MT_H2 * Hb';
    M4_upper = [M4_upper_left, M4_upper_right];

    M4_lower_left = -M4_upper_right';
    M4_lower_right = speye(nb);
    M4_lower = [M4_lower_left, M4_lower_right];

    M4 = sparse([M4_upper; M4_lower]);

    %Check if it is possible it initialize Q^{-1} with all first na ramps
    %being active
    A_upper_left = M_H2'*Hess*M_H2;
    if rank(full(A_upper_left)) == size(A_upper_left, 1)
        D_lower_right = M4_lower_right;
        C_lower_left = Hb*M_H2; 
        Q_inv_upper_left = inv(A_upper_left);
        Q_inv_lower_left = -D_lower_right*C_lower_left/A_upper_left;
        Qinit = [Q_inv_upper_left,zeros(na,nb);Q_inv_lower_left,D_lower_right];
        actset = sparse([ones(na,1);zeros(nb,1)]);
        disp("Starting with the first na ramps being active")
    else
        Qinit = speye(nc);
        actset = zeros(nc,1);
    end
end
