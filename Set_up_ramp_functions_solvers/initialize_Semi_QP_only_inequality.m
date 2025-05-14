function [Qinit,actset,M4,m5,Ha_inv,ha,na,nb,nc] = initialize_Semi_QP_only_inequality(SelectedFile,Q,c,H_ineq,h_ineq)
    
    use_already_filtered_matrices = true;
    % use_already_filtered_matrices = false;


    if use_already_filtered_matrices 
        try
            [Ha, Hb, ha, hb, Ha_inv, ~] = load_initialization_values_for_semi_QP_ineq(SelectedFile);
        catch
            error("A new algorithm is needed to create the invertible matrix [H_a^T]")
        end
    else
        disp("Finding the partition matrix manually")
        disp("This may take a while")
        na     = size(Q,1);
        nb     = size(H_ineq,1)-na ;
        finite_mask = isfinite(h_ineq);
        H_ineq2 = H_ineq(finite_mask, :);
        h_ineq2 = h_ineq(finite_mask);

        num_removed = length(h_ineq) - length(h_ineq2);
        if length(h_ineq2) > na %+inf
            h_ineq = h_ineq2;
            H_ineq = H_ineq2;
            fprintf("Removed %d Inf or -Inf values from h_ineq\n", num_removed);
        else
            new_max_value = 1e6;
            h_ineq(h_ineq == Inf) = new_max_value;
            h_ineq(h_ineq == -Inf) = -new_max_value;
        end
        [Ha, Hb, ha, hb, Ha_inv, idx] = select_independent_rows(H_ineq, h_ineq);
    end

    na = length(ha);
    nb = length(hb);
    nc = na+nb;

    new_extreme = 1e6;
    ha(ha ==  max( 1e6,max(ha))) =  new_extreme;
    ha(ha ==  min(-1e6,min(ha))) = -new_extreme;
    hb(hb ==  max( 1e6,max(hb))) =  new_extreme;
    hb(hb ==  min(-1e6,min(hb))) = -new_extreme;
    % disp("na = " + na)
    % disp("nb = " + nb)
    % disp(" ")

    [M4, m5,Qinit,actset] = Build_M4_m5(Ha_inv, Hb, Q, c, ha, hb);
end
