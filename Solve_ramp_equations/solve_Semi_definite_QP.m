function  [z,time,actset] = solve_Semi_definite_QP(Qmat0i,actset,M4,m5,M_H,ha,na,nc,g_eq)
    
    rank2 = false;
    qC = diag([-1,1]);


    %Methods for selecting iz:
    method = 1;

    %Method 1 - most extreme option that "doesnt" cause singularity,
    %prioritizing adding elements to the active set
    %Method 2 - most extreme option that "doesnt" cause singularity,
    %prioritizing removal of elements to the active set
    %Method 3 - original method

    y0 = Qmat0i*m5;

    tend = 1; 
    tic
    for ik = 1:tend
        feasflag = false;
        solved = false;
        ix = 0;
 
        if method == 3
            Min_qdiv = 1e-8;
        else
            Min_qdiv = 0; %Ensures no rank 2 updates
        end
 
        
        while (~solved)
            ix = ix +1;
            if (ix == 1)
                y = y0;
            elseif ~rank2
                y = y-y(iz)*vAd;
            else
                y = Qmat0i*m5;
            end

            if mod(ix,500) == 0 %To reset numerical errors
                Qmat0i_candidate = inv(speye(nc)-M4*diag(actset));
                if condest(Qmat0i_candidate) < 1e12
                    Qmat0i = Qmat0i_candidate;
                    y = Qmat0i*m5;
                end
            end   

            % Try to remove a constraint first
            actset0 = actset;
            [iz, actset, qc] = find_best_index_for_rank_update(method,M4, y, actset, Qmat0i);


            if(~isempty(iz))
                qu = M4(:,iz);
                vA = Qmat0i*qu;
                qdiv = qc+vA(iz);

                if abs(qdiv) < Min_qdiv
                    rank2 = true;
                else
                    rank2 = false;
                end

                if ~rank2 %Rank 1 update
                    vAd = (1/qdiv)*vA;
                    Qmat0i = sparse(Qmat0i-vAd * Qmat0i(iz, :));
                else %Rank 2 update
                    actix = find(actset0 == 1);
                    dy = Qmat0i * qu;
                    dy = clean(dy, 1e-6);
                    dyn = find(dy(actix) < 0);
                    
                    if isempty(dyn)
                        disp("dyn empty")
                        sum(full(actset))
                        flag = -1;
                        return 
                    end
                    
                    ratios = abs(y(actix(dyn)) ./ dy(actix(dyn)));
                    [~, sorted_indices] = sort(ratios,"ascend");

                    iz4_candidates = actix(dyn(sorted_indices));
                    iz4 = iz4_candidates(1); 
                    
                    % If iz4 is same as iz, try second-iz4 if it exists
                    if iz == iz4 && length(iz4_candidates) > 1
                        iz4 = iz4_candidates(2);
                        disp("Selected second-best because of conflict with iz")
                    elseif iz == iz4
                        disp("same (no alternative found)")
                    end
                    
                    actset(iz4) = 0;
                    qu2 = M4(:,iz4);
                    qU = [qu qu2];
                    Qmat0i = sparse(Qmat0i - Qmat0i * qU * 1/(qC + Qmat0i([iz, iz4],:)*qU) * Qmat0i([iz, iz4],:));
                end
            else
                solved = true;
                feasflag = true;
            end
        end

        
        if (feasflag == 0)
            break
        end
        
        need_to_account_for_numerical_issues = true;
        % need_to_account_for_numerical_issues = false;
        if isempty(g_eq)
            inv_Ha = M_H;
            if ~need_to_account_for_numerical_issues
                ya = y(1:na);
                ra = max(ya,0);
                z = inv_Ha*(ha-ra);
            else
                y = inv(speye(nc)-M4*diag(actset))*m5;
                % y = (speye(nc) - M4 * diag(actset)) \ m5;
                ya = y(1:na);
                ra = max(ya,0);
                z = inv_Ha*(ha-ra);
            end
        else
            if ~need_to_account_for_numerical_issues
                ya = y(1:na);
                ra = max(ya,0);
                z = M_H*[g_eq;ha-ra];
            else
                y = inv(speye(nc)-M4*diag(actset))*m5;
                % y = (speye(nc) - M4 * diag(actset)) \ m5;
                ya = y(1:na);
                ra = max(ya,0);
                z = M_H*[g_eq;ha-ra];
            end
        end

            
    end
    time = toc;
    
    % ix
    % na_actset = sum(full(actset(1:na)));
    % nb_actset = sum(full(actset(na+1:end)));
    % disp("Sum of all na indices from actset: " + na_actset)
    % disp("Sum of all nb indices from actset: " + nb_actset)
    % disp("Sum of actset " + sum_total_actset)
    % size(g_eq,1)
    % disp("Finished at iteration = " + ix )
end


function x = clean(x, tol)
    x(abs(x) < tol) = 0;
end


