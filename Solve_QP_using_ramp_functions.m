clear all
clc
addpath("Set_up_for_initialization")

qpramp_indexes = [41,42,43,44, 99,128];

model_index = 48;
[Q, c, H_ineq, h_ineq, A_eq, b_eq, c0, solution, SelectedFile, Main_folder] = get_model_data_from_git(model_index);
[z,time,~] = solve_Ramp_QP(Main_folder,SelectedFile, Q, c, H_ineq, h_ineq, A_eq, b_eq);
[z1,time1] = solve_with_quadprog(Q, c, H_ineq, h_ineq, A_eq, b_eq);


Ramp_solution     = 1/2*z' *Q*z  +c'*z  + c0;
Quadprog_solution = 1/2*z1'*Q*z1 +c'*z1 + c0;
Actual_solution   = solution;

disp("Ramp solver solution: " + Ramp_solution)
disp("quadprog solution   : " + Quadprog_solution)
disp("Actual solution     : " + solution)
disp(" ")
disp("Ramp solver time [ms] :" + round(time *1e3,1))
disp("quadprog    time [ms] :" + round(time1*1e3,1))
%%

function [Ramp_solution,time,actset] = solve_Ramp_QP(Main_folder,SelectedFile, Q, c, H_ineq, h_ineq, A_eq, b_eq)
    Give_time_estimate_for_finding_solution(SelectedFile,Main_folder)
    solver = determine_solver_to_use(Q, A_eq);
    
    %Might need to override the selected solver for some models. 
    %This can be done by simply uncomment the codes below
    if solver == "qpramp"
        solver = "SemiQP_ineq"
    end

    disp("Selected file  : " + SelectedFile)
    disp("Selected solver: " + solver)
    disp(" ")
    Activate_right_paths(Main_folder,solver);

    if solver == "qpramp"
        [neg_g_invh_gt_t,~,neg_s_t,neg_w_t,neg_g_invh_t,invh,x0] = initialize_qpramp(Q, c, H_ineq, h_ineq);
        %The threading code requires the "Parallel Computing Toolbox"
        threading_on = false;
        % threading_on = true;
        [z, time] = run_qpramp_solve(neg_g_invh_gt_t, neg_s_t, neg_w_t, neg_g_invh_t, x0, invh, c, threading_on);
        actset = [];
    end
    
    if solver == "SemiQP" %|| solver == "New_QP"
        [Qinit,actset,M4,m5,M_H,ha,na,~,nc] = intialize_Semi_def_QP_solver(SelectedFile, ...
            Q,c, H_ineq, h_ineq,A_eq,b_eq);
        [z,time,actset] = solve_Semi_definite_QP(Qinit,actset,M4,m5,M_H,ha,na,nc,b_eq);
    end

    if solver == "SemiQP_ineq"
        [Qinit,actset,M4,m5,Ha_inv,ha,na,~,nc] = initialize_Semi_QP_only_inequality(SelectedFile,Q,c,H_ineq,h_ineq);
        [z,time,actset] = solve_Semi_definite_QP(Qinit,actset,M4,m5,Ha_inv,ha,na,nc,b_eq);
    end
    if solver == "new_QP"
        [M4,m5, Qinv,actset,b_eq,Hi,GQGi,HQG,QH,QG] = initialize_New_QP(Q,c, ...
            A_eq,b_eq, H_ineq,h_ineq);
        [z,time,actset] = f_run_new_QP(M4,m5,Qinv,actset,b_eq,Hi,c,GQGi,HQG,QH,QG);
    end
    Ramp_solution = z;
end

function [z,time] = solve_with_quadprog(Q, c, H_ineq, h_ineq, A_eq, b_eq)
    options = optimoptions('quadprog', 'Display', 'off');
    tic
    z = quadprog(Q, c, H_ineq, h_ineq, A_eq, b_eq, [], [],[], options);
    time = toc; 
end

