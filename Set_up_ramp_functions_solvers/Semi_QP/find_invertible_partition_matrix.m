function [idx, Ha, Hb, hb, ha, M_HT] = find_invertible_partition_matrix(G_eq, H_ineq, h_ineq)
    m_eq = size(G_eq, 1);
    n    = size(G_eq, 2);
    na   = n - m_eq;

    Z = null(full(G_eq));   % Nullspace of G_eq
    rZ = size(Z, 2);
    if rZ < na
        error('Nullspace has lower dimension than required.');
    end

    H_proj = H_ineq * Z;

    % Modified Gram-Schmidt for picking na linearly independent rows
    idx = zeros(na, 1);
    Q = zeros(rZ, na);
    count = 0;
    threshold = 1e-10;

    for i = 1:size(H_proj,1)
        v = H_proj(i,:)';
        for j = 1:count
            v = v - (Q(:,j)' * v) * Q(:,j);
        end
        norm_v = norm(v);
        if norm_v > threshold
            count = count + 1;
            Q(:,count) = v / norm_v;
            idx(count) = i;
            if count == na
                break;
            end
        end
    end

    if count < na
        error('Could not find enough independent inequality constraints.');
    end

    remaining = setdiff(1:size(H_ineq,1), idx);  % rows not yet used

    max_retries = 1e4;
    retry = 0;

    while retry < max_retries
        Ha = H_ineq(idx,:);
        ha = h_ineq(idx,:);
        partition_mat = [G_eq; Ha]';

        % Check conditioning
        if condest(partition_mat) < 1e10
            % Good enough
            break;
        else
            % Replace worst row in Ha with a fresh one
            retry = retry + 1;
            replace_idx = randi(na);  % Randomly pick a row to replace
            new_row = remaining(1);   % Take first unused row
            remaining(1) = [];        % Remove it from pool
            idx(replace_idx) = new_row;

            if isempty(remaining)
                error('Ran out of unused inequality rows while trying to fix singularity.');
            end
        end
    end

    if retry >= max_retries
        error('Failed to build a well-conditioned partition after retries.');
    end
    % Final assembly
    mask = true(size(H_ineq,1),1);
    mask(idx) = false;
    Hb = H_ineq(mask,:);
    hb = h_ineq(mask,:);

    I = speye(n);
    % M_HT = inv(partition_mat);
    M_HT = zeros(n);
    for i = 1:n
        M_HT(:,i) = partition_mat \ I(:,i);
    end
end
