function [zsave,time,actset] = f_run_new_QP(M4,m5,Qinv,actset,g,Hi,c,GQGi,HQG,QH,QG)

    Qmat0i = Qinv;
    nc = length(m5);
    y0 = m5;
   
    method = 1;
    %Methods for selecting iz:
    %Method 1 - most extreme option that "doesnt" cause singularity,
    %prioritizing adding elements to the active set
    %Method 1 - most extreme option that "doesnt" cause singularity,
    %prioritizing removal of elements to the active set
    %Method 3 - original method

    rank2 = false;
    flag = 0;
    qC = [-1 0;0 1];


    tend = 1;
    tic
    for ik = 1:tend
        solved = false;
        ix = 0;
        iz = 0;

        if method == 3
            Min_qdiv = 1e-8;
        else
            Min_qdiv = 0;
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
              
                % disp("ix = " + ix)
                % disp("size of current actset = " + sum_actset)
                Qmat0i_candidate = inv(speye(nc)-M4*diag(actset));
                if condest(Qmat0i_candidate) < 1e12
                    Qmat0i = Qmat0i_candidate;
                    y = Qmat0i*m5;
                end
            end   
            actset0 = actset;
 
            [iz, actset, qc] = find_best_index_for_rank_update(method,M4, y, actset, Qmat0i);
            if ix > 20e3
                break
            end           
            if(~isempty(iz))
                qu = M4(:,iz);
                vA = Qmat0i*qu;
                qdiv = qc+vA(iz);
                if abs(qdiv) < Min_qdiv
                    rank2 = true;
                else
                    rank2 = false;
                end
   
                if ~rank2 
                    vAd = (1/qdiv)*vA;
                    Qmat0i = sparse(Qmat0i-vAd * Qmat0i(iz, :));
                else 
                    disp("Rank2")
                    actix = find(actset0 == 1);
                    dy = Qmat0i * qu;
                    dy = clean(dy, 1e-6);
                    dyn = find(dy(actix) < 0);
                    
                    if isempty(dyn)
                        qc
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
        if (flag == 0)
            % Qmat0i = inv(speye(nc)-M4.*actset);
            % y = Qmat0i*m5;
            lam = max(y,0);
            mu = -GQGi*(HQG'*lam+g);
            v = -QG*mu-QH*lam;
            zsave = v-Hi*c;
            time =  toc;
        else
            if flag == -1
                disp("Infeasible")
            else
                disp("Infinity loop")
            end
            zsave = [];
            time = [];
        end
    end
end



