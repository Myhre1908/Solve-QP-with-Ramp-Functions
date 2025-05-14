function [M4,m5, Qinv,actset,b_eq,Hi,GQGi,HQG,QH,QG] = initialize_New_QP(Hess,c,G,g,H,h)
    
    % Hess = Q
    % G = A_eq
    % g = b_eq
    % 
    % H = H_ineq
    % h = h_ineq


    %Remove inifinity constraints
    finite_mask = isfinite(h);
    H = H(finite_mask, :);
    h = h(finite_mask);
    nc = size(H,1);
    Hi = inv(Hess);

    [GQGi,HQG,HQ,QH,QG,M4] = create_matrices(Hi,G,H,nc);


    h = h + HQ*c;
    b_eq = g + G*Hi*c;

    m5 = HQG*GQGi*b_eq -h;


    actset = zeros(nc,1);
    Qinv = speye(nc); 


end

function [GQGi,HQG,HQ,QH,QG,M4] = create_matrices(Hi,G,H,nc)
    GQG = G*Hi*G'; 
    if ~isinf(det(GQG)) && det(GQG) > 1e-6
        GQGi = inv(GQG);
    elseif det(GQG) < 1e-6
        % GQGi = (inv(full(GQG)));
        GQGi =(pinv(full(GQG)));
        % disp("rank GQ{-1}G^T:" + rank(full(GQG)))
        % disp("Size GQ{-1}G^T:" + size(GQG))
        % disp("det(GQG) = very small. Expect numerical errors")
    else
        GQGi = (pinv(full(GQG)));
        disp("det(GQG) = inf. Expect numerical errors")
    end

    HQH = H*Hi*H'; HQG = H*Hi*G';
    HQ = H*Hi; QH = HQ';
    GQ = G*Hi; QG = GQ';
    M4 = sparse(speye(nc)-HQH+HQG*GQGi*HQG');
end
