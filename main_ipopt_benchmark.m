clear;clc;

params = setup_problem_parameters(); 

%%
import casadi.*

opti = casadi.Opti();

% Define the variable
x = opti.variable(params.n_o,1); % Community O-D flow vector 
y = opti.variable(params.n_l,1); % Link-flow vector 
z = opti.variable(params.n_r,1); % Route-flow vector 
b = opti.variable(params.n_c,1); % Community benefit vector
n = opti.variable(params.n_c, 1);
n_prime = opti.variable(params.n_c, 1);
t = opti.variable(params.n_c, 1);
%log_arg = opti.variable(params.n_c, 1); %
u = opti.variable(params.n_c, 1);
t_min_active = opti.variable(1, 1); 
objective = opti.variable(1, 1); 


% Network Flow Constraints 
opti.subject_to(params.E * y == zeros(params.n_n, 1));
opti.subject_to(params.F * z == y);
%opti.subject_to(x == (params.P_scaled / params.P_ScalingFactor) * z);
P_matrix_for_constraint = sparse(params.P_scaled / params.P_ScalingFactor);
P_param = opti.parameter(params.n_o, params.n_r);
opti.set_value(P_param, P_matrix_for_constraint);
opti.subject_to(x == P_param * z);

% Capacity
opti.subject_to(params.J * z <= (1 - params.ep) * params.c_v);
opti.subject_to(y <= (1 - params.ep) * params.c_l);
opti.subject_to(params.K * y <= (1 - params.ep) * params.c_w);

% Non-negativity and bounds
opti.subject_to(y >= 0);
opti.subject_to(z >= 0);
opti.subject_to(x >= 0);

opti.subject_to(b <= 1); 
opti.subject_to(b >= 0);

opti.subject_to(n >= 0);

opti.subject_to(n_prime >= 0);
opti.subject_to(n_prime <= params.Delta_n_max);

opti.subject_to(t_min_active >= 0);

% Energy
opti.subject_to(params.p' * z <= params.max_energy * sum(z));

% Community Benefit 
W_abs = abs(params.W); 
W_abs_row_sums = sum(W_abs, 2); 
zero_sum_indices = find(W_abs_row_sums == 0);
if ~isempty(zero_sum_indices)   
    W_abs_row_sums(zero_sum_indices) = 1000;  % to avoid divide-by-zero
end
% D_inv_W_abs = sparse(diag(1 ./ W_abs_row_sums)) * W_abs; 
% opti.subject_to(b == D_inv_W_abs * (x ./ params.e));
D_inv = spdiags(1 ./ W_abs_row_sums, 0, params.n_c, params.n_c);
D_inv_W_abs = D_inv * W_abs;
opti.subject_to(b == D_inv_W_abs * (x ./ params.e));

% Noise 
%opti.subject_to(log_arg == params.M' * y + 1e-6);
%opti.subject_to(log_arg >= 1e-6);
%opti.subject_to(n == 10 * log10(log_arg) - 10 * log10(params.T_ratio));
opti.subject_to(params.M' * y + 1e-6 == params.T_ratio * 10.^(n / 10)); % <-- Use this line
opti.subject_to(n_prime >= n - params.a);

% Cost-Benefit Analysis
opti.subject_to(t == b .* (1 - n_prime / params.Delta_n_max));

% SWF
opti.subject_to(u <= t_min_active);
for i = 1:params.n_c
    if ismember(i, params.active_communities_idx) 
        opti.subject_to(t_min_active <= t(i));
        opti.subject_to(u(i) <= t(i) - params.Delta);
    end
end

% objective
opti.subject_to(objective == params.Delta + (1 / params.n_c_active) * sum(u));


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial-guess

opti.set_initial(x, 0.9 * ones(params.n_o, 1));
opti.set_initial(y, 0.1 * ones(params.n_l, 1));
opti.set_initial(z, 1 * ones(params.n_r, 1));
opti.set_initial(b, 0.5 * ones(params.n_c, 1));
opti.set_initial(n, 90 * ones(params.n_c, 1));
opti.set_initial(n_prime, 40 * ones(params.n_c, 1));
opti.set_initial(t, 0.05 * ones(params.n_c, 1));
opti.set_initial(u, 0.02 * ones(params.n_c, 1));
opti.set_initial(t_min_active, 0.02);
opti.set_initial(objective, 0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%


% Objective Function

opti.minimize(-objective);

% Solver options
opts = struct();
opts.verbose = false;
opts.ipopt.print_level = 5;
opts.print_time = false;

opts.ipopt.warm_start_init_point = 'yes';
opts.ipopt.max_iter = 1000;

opts.ipopt.tol = 1e-6;
opts.ipopt.acceptable_tol = 1e-4;
opts.ipopt.acceptable_obj_change_tol = 1e-5; 
opts.ipopt.mu_strategy = 'adaptive';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.ipopt.hessian_approximation = 'limited-memory';  % faster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


opts.ipopt.linear_solver = 'mumps';                   % default stable solver


opti.solver('ipopt', opts);
sol = opti.solve();



x = sol.value(x);
y = sol.value(y);
z = sol.value(z);
b = sol.value(b);
n = sol.value(n);
n_prime = sol.value(n_prime);
u = sol.value(u);
t_min_active = sol.value(t_min_active);
t = sol.value(t);
fval = value(objective);
flag = sol.problem;


