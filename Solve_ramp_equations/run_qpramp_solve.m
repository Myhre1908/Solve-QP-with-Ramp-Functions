function [z, elapsedTime] = run_qpramp_solve(neg_g_invh_gt_t, neg_s_t, neg_w_t, neg_g_invh_t, x0, invh, ...
    c,Use_threading_to_avoid_infinity_loops)
    
    clear qpramp_solve
    if Use_threading_to_avoid_infinity_loops
        if isempty(gcp('nocreate'))
            parpool;
        end
        % Submit the job and start timing the solve
        solveStart = tic;
        f = parfeval(@qpramp_solve, 1, neg_g_invh_gt_t, neg_s_t, neg_w_t, neg_g_invh_t, x0);
        
        % Timeout setup
        timeout = 10;
        tStart = tic;
    
        % Wait loop
        while ~strcmp(f.State, 'finished') && toc(tStart) < timeout
            pause(0.1);
        end
    
        % Handle result or timeout
        if strcmp(f.State, 'finished')
            try
                gamma = fetchOutputs(f);
                z = gamma' - invh * c;
                elapsedTime = toc(solveStart);  % Time from just before solve
            catch 
                warning("Function finished with error:\n%s", f.Error.message);
                z = [];
                elapsedTime = NaN;
            end
        else
            cancel(f);
            warning("Function execution timed out after %.1f seconds.", timeout);
            z = [];
            elapsedTime = NaN;
        end
    else
        disp("Matlab might crash unless you turn on the threading code")
        try
            tic
            gamma = qpramp_solve(neg_g_invh_gt_t, neg_s_t, neg_w_t, neg_g_invh_t, x0);
            z = gamma' - invh * c;
            elapsedTime = toc;
        catch err
            disp("HERE")
            warning(sprintf('Function finished with error:\n%s', err.message));
            % warning(err.identifier, 'Function finished with error:\n%s', err.message);
            z = [];
            elapsedTime = [];
        end
    end

end

