function [Qinit,actset,M4,m5,M_H,ha,na,nb,nc,H_ineq_new,h_ineq_new]= intialize_Semi_def_QP_solver(selectedFile,Hess,c, H_ineq, h_ineq,G_eq,g_eq)

    use_already_filtered_matrices = true;
    % use_already_filtered_matrices = false;

    if use_already_filtered_matrices 
        try
            [Ha, Hb, ha, hb, M_HT, selected_indices] = load_initialization_values_for_semi_QP(selectedFile);
        catch
            error("A new algorithm is needed to create the invertible matrix [M^T (H_a)^T]")
        end
    else
        disp("Finding the invertible matrix [M^T (H_a)^T] manually. This may take a while")
        %Just to check if there exist one where inf => 1e6;

        m_eq = size(G_eq, 1);
        k     = size(G_eq, 2);
        na     = k - m_eq;
        finite_mask = isfinite(h_ineq);
        H_ineq2 = H_ineq(finite_mask, :);
        h_ineq2 = h_ineq(finite_mask);
        num_removed = length(h_ineq) - length(h_ineq2);
        if length(h_ineq2) > na
            h_ineq = h_ineq2;
            H_ineq = H_ineq2;
            fprintf("Removed %d Inf or -Inf values from h_ineq\n", num_removed);
        else
            new_extreme = 1e9;
            h_ineq(h_ineq ==  Inf) =  new_extreme;
            h_ineq(h_ineq == -Inf) = -new_extreme;
        end
        [idx, Ha, Hb, hb, ha, M_HT] = find_invertible_partition_matrix(G_eq, H_ineq, h_ineq);
    end
    M_H = M_HT';

    new_extreme = 1e9;
    ha(ha ==  max( 1e6,max(ha))) =  new_extreme;
    ha(ha ==  min(-1e6,min(ha))) = -new_extreme;
    hb(hb ==  max( 1e6,max(hb))) =  new_extreme;
    hb(hb ==  min(-1e6,min(hb))) = -new_extreme;
    na = length(ha);
    nb = length(hb);
    nc = na+nb;
   
    % disp("na = " + na)
    % disp("nb = " + nb)
    % disp(" ")

    [M4,m5,Qinit,actset] = Build_M4_m5(M_HT, Hb, Hess, c, ha, hb, g_eq);

   
end