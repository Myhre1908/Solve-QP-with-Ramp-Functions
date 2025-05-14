function [M4, m5,Qinit,actset] = Build_M4_m5(M_HT, Hb, Hess, c, ha, hb, g_eq)
    na = size(ha,1);
    nb = size(hb,1);
    nc = na+nb;
    M_H = M_HT';
    
    M_H2 = M_H(:, end-na+1:end);
    MT_H2 = M_H2';
    
    bHa = [g_eq;ha];
    
    m5_upper = MT_H2*(Hess*M_H*bHa+c);
    m5_lower = Hb*M_H*bHa-hb;
    m5 = sparse([m5_upper;m5_lower]);

    M4_upper_left = speye(na)-MT_H2*Hess*M_H2;
    M4_upper_right = MT_H2*Hb';
    
    M4_upper = [M4_upper_left, M4_upper_right];
    M4_lower_left = -Hb*M_H2;
    M4_lower_right = speye(nb);
    M4_lower = [M4_lower_left,M4_lower_right];
    M4 = sparse([M4_upper;M4_lower]);

    A_upper_left = M_H2'*Hess*M_H2;
    if rank(A_upper_left) == size(A_upper_left, 1)
        %This is Q_inv when the first na ramps are active
        D_lower_right = M4_lower_right;
        C_lower_left = Hb*M_H2;
        Q_inv_upper_left = inv(A_upper_left);
        Q_inv_lower_left = -D_lower_right*C_lower_left*Q_inv_upper_left;
       
        Qinit = [Q_inv_upper_left,zeros(na,nb);Q_inv_lower_left,D_lower_right];
        actset = [ones(na,1);zeros(nb,1)];
        disp("Starting with the first na ramps being active")
        disp(" ")
    else
        actset = zeros(nc,1);
        Qinit = speye(nc);
        
    end
end