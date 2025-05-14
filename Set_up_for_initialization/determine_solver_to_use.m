function solver = determine_solver_to_use(Q,A_eq)
    warning("off")
    if isempty(A_eq)
        only_inequality_constraits = true;
    else
        only_inequality_constraits = false;
    end
    
    min_eig = eigs(Q, 1, 'SA');  % Smallest algebraic eigenvalue
    if isnan(min_eig)
        e = eig(full(Q));
        min_eig = min(real(e));
    end
    detQ = det(Q);
    
% [R, p] = chol(Q);
% if p == 0
%     fprintf('Q is positive definite.\n');
% else
%     fprintf('Q is not strictly positive definite (failed at pivot %d).\n', p);
% end

    if min_eig > 0 && detQ > 1e-4
        positive_definite = true;
        positive_semidefinite = false;
    elseif min_eig >= -1e-4 && (detQ < 1e-4 || isnan(detQ))
        positive_semidefinite = true;
        positive_definite = false;
    else
        disp('H is not positive semidefinite');
        solver = [];
        return
    end
    
    if only_inequality_constraits && positive_definite
        solver = "qpramp";
    elseif ~only_inequality_constraits && positive_definite
        if isfinite(det(A_eq*inv(Q)*A_eq')) 
            solver = "new_QP";
        else
            solver = "SemiQP";
        end
    elseif ~only_inequality_constraits && positive_semidefinite
        solver = "SemiQP";
    elseif only_inequality_constraits && ~positive_definite
        solver = "SemiQP_ineq";
    end
    warning("on")
end