function [neg_g_invh_gt_t,neg_invh_f,neg_s_t,neg_w_t,neg_g_invh_t,invh,x0] = initialize_qpramp(Q, c, H_ineq, h_ineq)
    finite_mask = isfinite(h_ineq);
    H_ineq = H_ineq(finite_mask, :);
    h_ineq = h_ineq(finite_mask);

    x0 = ones(1,length(c));
    f0 = c * pinv(x0');
    g = H_ineq;

    invh = inv(full(Q));
    s = zeros(length(h_ineq),length(c));
    w = h_ineq + H_ineq/Q*f0*x0';


    neg_g_invh_gt = -g*invh*g';
    neg_s = -s;
    neg_w = -w;
    neg_g_invh = -g*invh;
    neg_invh_f = -invh*f0;
    
    % Transposing in the loop gives big slowdown
    neg_g_invh_gt_t = neg_g_invh_gt';
    neg_s_t = neg_s';
    neg_w_t = neg_w';
    neg_g_invh_t = neg_g_invh';
end
