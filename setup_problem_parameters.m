function params = setup_problem_parameters()
    % This function loads all necessary matrices and parameters from files.
    % It assumes the files are in a subdirectory named 'input matrix_od30_r1.2'.

    fprintf('Loading problem parameters...\n');
    folder_path = 'input matrix_od30_r1.2';

    % Load matrices
    params.E = sparse(readmatrix(fullfile(folder_path, 'matrix_E.xlsx')));
    params.F = sparse(readmatrix(fullfile(folder_path, 'matrix_F.xlsx')));
    params.J = sparse(readmatrix(fullfile(folder_path, 'matrix_J.xlsx')));
    params.H = max(params.J,0);
    params.K = max(params.E,0);
    
    % Load .mat file for sparse matrix P
    temp_P = load(fullfile(folder_path, 'matrix_P.mat'));
    params.P = temp_P.P;
    params.P_ScalingFactor = 1e5; % 选择一个合适的缩放因子，例如 10^4 或 10^5
    params.P_scaled = params.P * params.P_ScalingFactor; % 创建缩放后的 P 矩阵

    params.M = readmatrix(fullfile(folder_path, 'matrix_M.xlsx'));       
    params.M_ScalingFactor = 1e7; % 新增 M 矩阵的数值缩放因子
    params.M_scaled = params.M / params.M_ScalingFactor; % 创建缩放后的 M 矩阵
    threshold_M = 1.1 / params.M_ScalingFactor;% 对 M_scaled 进行阈值化
    params.M_scaled(params.M_scaled < threshold_M) = 0; 
    params.M_scaled = sparse(params.M_scaled); % 将缩放后的 M_scaled 转换为稀疏
    params.M = sparse(params.M); % 如果M是结构稀疏的，这里仍然需要 sparse()

    params.W = sparse(readmatrix(fullfile(folder_path, 'matrix_W_transposed.xlsx'))');
    
    % Get dimensions from loaded matrices
    [params.n_n, params.n_l] = size(params.E);
    [~, params.n_r] = size(params.F);
    [params.n_v, ~] = size(params.J);
    [params.n_o, ~] = size(params.P);
    [~, params.n_c] = size(params.M_scaled);% 使用 M_scaled 的维度

    % Load vectors and scalar parameters
    params.c_v = readmatrix(fullfile(folder_path, 'vector_cv.xlsx'));
    params.c_l = readmatrix(fullfile(folder_path, 'vector_cl.xlsx'));
    params.c_w = readmatrix(fullfile(folder_path, 'vector_cw.xlsx'));
    params.a = readmatrix(fullfile(folder_path, 'vector_a.xlsx'));
    params.e = readmatrix(fullfile(folder_path, 'vector_e.xlsx'));
    params.p = readmatrix(fullfile(folder_path, 'vector_p.xlsx'));
    
    % Design Parameters
    params.Delta_n_max = 50;
    params.T_ratio = 35.56;
    params.ep = 0.01;
    params.max_energy = 20;
    params.Delta = 0; % Example threshold for SWF, adjust as needed

    fprintf('Parameters loaded successfully.\n');
    
    % 向量参数
    fprintf('\n向量参数统计:\n');
    vectors_to_check = {'c_v', 'c_l', 'c_w', 'a', 'e', 'p'};
    for i = 1:length(vectors_to_check)
        vec_name = vectors_to_check{i};
        vec_data = params.(vec_name); % 动态获取参数数据

        % 确保数据是数值类型且非空
        if ~isnumeric(vec_data) || isempty(vec_data)
            fprintf('  %s: 非数值或为空。\n', vec_name);
            continue;
        end

        vec_min = min(vec_data);
        vec_max = max(vec_data);
        vec_mean = mean(vec_data);
        fprintf('  params.%s: min=%g, max=%g, mean=%g\n', vec_name, vec_min, vec_max, vec_mean);
    end


    % 矩阵参数 (非零元素统计)
    fprintf('\n矩阵参数统计 (非零元素):\n');
    matrices_to_check = {'E', 'F', 'J', 'K', 'P_scaled', 'M', 'W'}; % 使用P_scaled
    for i = 1:length(matrices_to_check)
        mat_name = matrices_to_check{i};
        mat_data = params.(mat_name);

        if ~isnumeric(mat_data) || isempty(mat_data)
            fprintf('  %s: 非数值或为空。\n', mat_name);
            continue;
        end

        % 对于稀疏矩阵，只考虑非零元素
        if issparse(mat_data)
            non_zero_elements = nonzeros(mat_data);
            if isempty(non_zero_elements)
                fprintf('  params.%s: (全零稀疏矩阵)\n', mat_name);
                continue;
            end
            mat_min = min(non_zero_elements);
            mat_max = max(non_zero_elements);
            mat_mean = mean(non_zero_elements);
            fprintf('  params.%s (sparse, non-zeros): min=%g, max=%g, mean=%g, nnz=%d\n', mat_name, mat_min, mat_max, mat_mean, nnz(mat_data));

            % 对于大型矩阵，计算条件数可能很慢或耗内存。可以跳过或仅对小矩阵计算
            % 如果矩阵不是太离谱的大，可以尝试计算条件数，但这只适用于数值范围检查
            % 如果是巨大的稀疏矩阵，cond()会很慢，不要轻易尝试
            % if all(size(mat_data) < 1000) % 仅对小矩阵计算条件数
            %     mat_cond = cond(full(mat_data)); % 对稀疏矩阵，先转为全矩阵再计算条件数
            %     fprintf('    条件数: %g\n', mat_cond);
            % end

        else % 全矩阵 (P_scaled, M)
            mat_min = min(mat_data(:)); % 使用 (:) 将矩阵展平为向量
            mat_max = max(mat_data(:));
            mat_mean = mean(mat_data(:));
            fprintf('  params.%s (full): min=%g, max=%g, mean=%g\n', mat_name, mat_min, mat_max, mat_mean);

            % 对于全矩阵，可以尝试计算条件数
            % 如果矩阵太大 (P_scaled: 4677x4050)，cond() 会非常非常慢且耗内存，可能崩溃。
            % if all(size(mat_data) < 1000) % 仅对小矩阵计算条件数
            %    mat_cond = cond(mat_data);
            %    fprintf('    条件数: %g\n', mat_cond);
            % end
        end
    end

fprintf('--------------------------------------\n');
end