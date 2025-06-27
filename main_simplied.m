clear;clc;

% 设置问题参数
params = setup_problem_parameters();

% 存储结果
results = struct();

%% 1. 求解
try
    fprintf('\n1. 尝试ipopt求解器...\n');
    tic;
    [sol_ipopt, fval_ipopt, flag_ipopt, output_ipopt] = simplied_ipopt(params);
    time_ipopt = toc;

    results.ipopt.solution = sol_ipopt;
    results.ipopt.fval = fval_ipopt;
    results.ipopt.flag = flag_ipopt;
    results.ipopt.time = time_ipopt;
    results.ipopt.output = output_ipopt;

    if flag_ipopt == 1
        fprintf('求解成功，目标值: %.4f，用时: %.2f秒\n', fval_ipopt, time_ipopt);
    else
        fprintf('求解失败或未找到最优解。 YALMIP状态: %s\n', output_ipopt.info);
    end
catch ME
    fprintf('不可用或发生错误: %s\n', ME.message);
    results.ipopt = [];
end


%%
function [solution, fval, exitflag, sol] = simplied_ipopt(params)    
    x = sdpvar(params.n_o, 1, 'full'); % Community O-D flow vector 
    y = sdpvar(params.n_l, 1, 'full'); % Link-flow vector 
    z = sdpvar(params.n_r, 1, 'full'); % Route-flow vector 
    b = sdpvar(params.n_c, 1, 'full'); % Community benefit vector
    
    % 约束条件
    constraints = [];
    
    %% 1.1.3 Network Flow Constraints 
    constraints = [constraints, params.E * y == zeros(params.n_n, 1)]; % Ey=0_nn 
    constraints = [constraints, params.F * z == y];                     % Fz=y 
    constraints = [constraints, x == (params.P_scaled / params.P_ScalingFactor) * z]; % x=Pz
    
    %% 1.1.4 Capacity Constraints 
    constraints = [constraints, params.J * z <= (1 - params.ep) * params.c_v]; % Jz <= (1-eps)cv 
    constraints = [constraints, y <= (1 - params.ep) * params.c_l];           % y <= (1-eps)cl 
    constraints = [constraints, params.K * y <= (1 - params.ep) * params.c_w]; % Ky <= (1-eps)cw 

    %% Non-negativity and bounds
    constraints = [constraints, y >= 0, z >= 0, x >= 0]; % y>=0, z>=0, x>=0 (implied)
    constraints = [constraints, b >= 0, b <= 1]; % b in range [0,1]
    
    %% 1.4 Energy Consumption (32i)
    constraints = [constraints, params.p' * z <= params.max_energy * sum(z)];
    
    %% 1.5.1 Community Benefit (26, 32g)
    W_abs = abs(params.W); 
    W_abs_row_sums = sum(W_abs, 2); 
    zero_sum_indices = find(W_abs_row_sums == 0);
    if ~isempty(zero_sum_indices)   
        W_abs_row_sums(zero_sum_indices) = 1000; 
    end
    D_inv_W_abs = sparse(diag(1 ./ W_abs_row_sums)) * W_abs; 
    constraints = [constraints, b == D_inv_W_abs * (x ./ params.e)]; 

    
    %% 1.5.4 Social Welfare Function (Objective) - 简化为 LP 目标
    objective = sum(b);
    
    %% 求解设置
    options = sdpsettings('solver', 'ipopt', 'verbose', 2, 'debug', 1);
    %options.ipopt.hessian_approximation = 'bfgs'; 
    options.ipopt.tol = 1e-6; 
    options.ipopt.print_level = 10;
    options.ipopt.output_file = 'ipopt_detailed_log.txt'; 
    options.ipopt.max_cpu_time = 7200;
    %options.ipopt.linear_solver = 'mumps';
    
    % YALMIP默认最小化,所以最大化目标时加负号
    sol = optimize(constraints, -objective, options);
    
    %% 结果提取
    if sol.problem == 0 || sol.problem == 4 
        solution.x = value(x);
        solution.y = value(y);
        solution.z = value(z);
        solution.b = value(b);
        % ... 其他被注释的变量不再返回 ...
        
        fval = value(objective);
        exitflag = 1;
    else
        solution = [];
        fval = -inf;
        exitflag = 0;
        warning('Optimization failed or did not find an optimal solution. YALMIP status: %s', sol.info);
    end

end