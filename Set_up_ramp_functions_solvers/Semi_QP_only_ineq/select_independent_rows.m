%This is mostly AI generated
function [Ha, Hb, ha, hb, Ha_inv,selected_indices] = select_independent_rows(H_ineq, h_ineq)
    % Selects n well-conditioned, linearly independent rows from H_ineq (m x n)
    % such that Ha^T is invertible and not ill-conditioned.
    % Returns:
    %   Ha, Hb        - Selected and remaining rows of H_ineq
    %   ha, hb        - Matching entries in h_ineq
    %   selected_indices - Indices of selected rows
    %   Ha_inv        - Stable inverse of Ha using QR

    [m, n] = size(H_ineq);
    Ha = zeros(n, n);
    selected_indices = zeros(1, n);
    count = 0;

    cond_thresh = 1e8;  % Acceptable condition number

    for i = 1:m
        candidate = full(H_ineq(i, :));
        if count == 0
            Ha(1, :) = candidate;
            selected_indices(1) = i;
            count = 1;
        else
            temp = [Ha(1:count, :); candidate];
            if sparse_rank(temp) > count
                temp_full = temp;
                % Only continue if not ill-conditioned
                if count + 1 < n || cond(temp_full) < cond_thresh
                    count = count + 1;
                    Ha(count, :) = candidate;
                    selected_indices(count) = i;
                end
            end
        end
        if count == n
            break;
        end
    end

    if count < n
        error('Could not find %d well-conditioned linearly independent rows.', n);
    end

    % Remaining rows
    all_indices = 1:m;
    remaining_indices = setdiff(all_indices, selected_indices);

    Hb = full(H_ineq(remaining_indices, :));
    hb = full(h_ineq(remaining_indices));
    ha = full(h_ineq(selected_indices));
    % Ha_inv = pinv(Ha);
    

    % Stable inverse of Ha using QR
    [Q, R] = qr(Ha);
    Ha_inv = R \ Q';
end




function r = sparse_rank(A, tol)
    if nargin < 2
        tol = 1e-10; % Default tolerance
    end
    [~, R] = qr(A, 0);
    r = sum(abs(diag(R)) > tol);
end