function [iz, actset, qc] = find_best_index_for_rank_update(method,M4, y, actset, Qmat0i)


    if method == 1 %Most extreme value
        [iz, actset, qc] = find_safe_set_to_add_most_extreme_not_singular(M4, y, actset, Qmat0i);
        if isempty(iz)
            [iz, actset, qc] = find_safe_set_to_remove_most_extreme_not_singular(M4, y, actset, Qmat0i);
        end

    elseif method == 2
        [iz, actset, qc] = find_safe_set_to_remove_most_extreme_not_singular(M4, y, actset, Qmat0i);
        if isempty(iz)
            [iz, actset, qc] = find_safe_set_to_add_most_extreme_not_singular(M4, y, actset, Qmat0i);
        end

    elseif method == 3
        [iz, actset, qc] = find_set_to_remove_most_extreme(y,actset);
        if isempty(iz)
            [iz, actset, qc] = find_set_to_add_most_extreme(y,actset);
        end    
    end
end
%% This ensures that only a rank 1 update is requires
function is_singular = check_singular_update(M4, i, Q_inv, qc)
    % qc = -1 for add, +1 for remove
    u      = M4(:,i);
    xi_row = Q_inv(i,:);
    if qc < 0
        qdiv = 1 - xi_row * u;        % “add” test
    else
        qdiv = 1 + xi_row * u;        % “remove” test
    end
    is_singular = abs(qdiv) < 1e-8;
    % is_singular = abs(qdiv) < 1e-8 || abs(qdiv) > 1e6;
end
%% Most extreme value
function [iz, actset, qc] = find_safe_set_to_remove_most_extreme_not_singular(M4, y, actset, Q_inv)
% Prioritizes most negative y(i) * actset(i) for removal

    iz = [];
    qc = 1;  % removal

    % Only consider elements where actset is 1 and y_remove is negative
    y_remove = y .* actset;
    % candidates = find((y_remove < 0) & (actset == 1));
    candidates = find((y_remove < 0));% & (actset == 1));

    % Sort candidates by most negative first
    [~, sorted_idx] = sort(y_remove(candidates));  % ascending order
    candidates = candidates(sorted_idx);

    for k = 1:length(candidates)
        iz_try = candidates(k);
        if ~check_singular_update(M4, iz_try, Q_inv, qc)
            iz = iz_try;
            actset(iz) = 0;
            return;
        end
    end
    % if k ~= 1
    %     disp("Selected " + k)
    %     % [iz, actset, qc] = find_safe_set_to_remove_first_possible(M4, y, actset, Q_inv);
    % end
end

function [iz, actset, qc, first] = find_safe_set_to_add_most_extreme_not_singular(M4, y, actset, Q_inv)
% Prioritizes most positive y(i) * (1 - actset(i)) for addition
    first = false;

    iz = [];
    qc = -1;  % addition

    % Only consider elements where actset is 0 and y_add is positive
    y_add = y .* (1 - actset);
    % candidates = find((y_add > 0) & (actset == 0));
    candidates = find((y_add > 0));% & (actset == 0));

    % Sort candidates by most positive first
    [~, sorted_idx] = sort(y_add(candidates), 'descend');
    candidates = candidates(sorted_idx);

    for k = 1:length(candidates)
        iz_try = candidates(k);
        if ~check_singular_update(M4, iz_try, Q_inv, qc)
            iz = iz_try;
            actset(iz) = 1;
            return;
        end
    end
    % if k ~= 1
    %     disp("Selected " + k)
    %     % [iz, actset, qc] = find_safe_set_to_add_first_possible(M4, y, actset, Q_inv);
    % end
end



%% Just most_extreme
function [iz, actset, qc] = find_set_to_remove_most_extreme(y,actset)
    lam = y.*actset;

    [i1,iz] = min(lam);
    if i1 < 0
        actset(iz) = 0;
        qc = 1; 
    else
        iz = [];
        qc = [];
    end
end


function [iz, actset, qc] = find_set_to_add_most_extreme(y,actset)
        lam = y.*actset;
        [i2,i2z] = max(y-lam);
        if (i2 <= 0)
            i2 = [];
        end

        if(~isempty(i2))
            iz  = i2z;
            actset(iz) = 1;
            qc = -1; % Add constraint to active set
        else
            iz = [];
            qc = [];
        end
    end
