function [u_null_expr, v_null_expr] = identify_nullclines_general(u, v, du, dv, options)
    % 通用动力系统nullcline识别工具
    % 新增选项:
    %   u_extra_basis - u-nullcline额外基函数 (函数句柄元胞数组, 默认{})
    %   v_extra_basis - v-nullcline额外基函数 (函数句柄元胞数组, 默认{})
    
    % 初始化输出变量 (防止未赋值错误)
    u_null_expr = '0';
    v_null_expr = '0';
    
    % 设置默认选项
    if nargin < 5
        options = struct();
    end
    if ~isfield(options, 'deriv_threshold'), options.deriv_threshold = 0.05; end
    if ~isfield(options, 'u_poly_order'), options.u_poly_order = 3; end
    if ~isfield(options, 'v_poly_order'), options.v_poly_order = 3; end
    if ~isfield(options, 'lambda'), options.lambda = 0.1; end
    if ~isfield(options, 'plot_results'), options.plot_results = true; end
    % 新增额外基函数选项
    if ~isfield(options, 'u_extra_basis'), options.u_extra_basis = {}; end
    if ~isfield(options, 'v_extra_basis'), options.v_extra_basis = {}; end
    
    % 验证输入数据
    validate_inputs(u, v, du, dv);
    
    % 识别nullcline点并存储导数
    [u_null_points, v_null_points, u_null_derivs, v_null_derivs] = ...
        find_nullcline_points(u, v, du, dv, options.deriv_threshold);
    
    % 检查是否有足够的点进行拟合
    if isempty(u_null_points) || size(u_null_points, 1) < 5
        warning('u-nullcline点不足 (%d)，使用默认表达式 v=0', size(u_null_points, 1));
        u_null_expr = '0';
    else
        % 使用SINDy拟合u-nullcline表达式 (v = f(u))
        u_null_expr = fit_nullcline_with_sindy(...
            u_null_points(:,1), u_null_points(:,2), 'u', ...
            options.u_poly_order, options.lambda, options.u_extra_basis);
    end
    
    % 检查是否有足够的点进行拟合
    if isempty(v_null_points) || size(v_null_points, 1) < 5
        warning('v-nullcline点不足 (%d)，使用默认表达式 u=0', size(v_null_points, 1));
        v_null_expr = '0';
    else
        % 使用SINDy拟合v-nullcline表达式 (u = g(v))
        v_null_expr = fit_nullcline_with_sindy(...
            v_null_points(:,1), v_null_points(:,2), 'v', ...
            options.v_poly_order, options.lambda, options.v_extra_basis);
    end
    
    % 打印结果
    fprintf('\n===== Nullcline 表达式 =====\n');
    fprintf('u-nullcline: v = %s\n', u_null_expr);
    fprintf('v-nullcline: u = %s\n', v_null_expr);
    
    % 可视化结果
    if options.plot_results
        visualize_results(u, v, u_null_points, v_null_points, ...
                          u_null_derivs, v_null_derivs, u_null_expr, v_null_expr);
    end
    
    % 嵌套函数定义
    function validate_inputs(u, v, du, dv)
        % 验证输入数据维度
        if numel(u) ~= numel(v) || numel(u) ~= numel(du) || numel(u) ~= numel(dv)
            error('所有输入向量(u, v, du, dv)必须有相同的长度');
        end
        
        % 转换为列向量
        u = u(:);
        v = v(:);
        du = du(:);
        dv = dv(:);
        
        % 检查NaN或Inf值
        if any(isnan(u)) || any(isinf(u)) || ...
           any(isnan(v)) || any(isinf(v)) || ...
           any(isnan(du)) || any(isinf(du)) || ...
           any(isnan(dv)) || any(isinf(dv))
            error('输入数据包含NaN或Inf值');
        end
    end

    function [u_null_points, v_null_points, u_null_derivs, v_null_derivs] = ...
        find_nullcline_points(u, v, du, dv, deriv_threshold)
        % 识别nullcline点
        
        % 识别u-nullcline点 (du/dt ≈ 0)
        u_null_idx = abs(du) < deriv_threshold;
        u_null_points = [u(u_null_idx), v(u_null_idx)];
        u_null_derivs = [du(u_null_idx), dv(u_null_idx)]; % 存储导数
        
        % 识别v-nullcline点 (dv/dt ≈ 0)
        v_null_idx = abs(dv) < deriv_threshold;
        v_null_points = [u(v_null_idx), v(v_null_idx)];
        v_null_derivs = [du(v_null_idx), dv(v_null_idx)]; % 存储导数
        
        % 显示识别到的点数
        fprintf('识别到 %d 个 u-nullcline 点 (|du/dt| < %.3f)\n', size(u_null_points, 1), deriv_threshold);
        fprintf('识别到 %d 个 v-nullcline 点 (|dv/dt| < %.3f)\n', size(v_null_points, 1), deriv_threshold);
        
        % 检查是否有足够的点
        if size(u_null_points, 1) < 5
            warning('u-nullcline点数量较少 (%d)，拟合结果可能不准确', size(u_null_points, 1));
        end
        if size(v_null_points, 1) < 5
            warning('v-nullcline点数量较少 (%d)，拟合结果可能不准确', size(v_null_points, 1));
        end
    end

    function nullcline_expr = fit_nullcline_with_sindy(x, y, nullcline_type, poly_order, lambda, extra_basis)
        % 根据nullcline类型确定输入输出
        if strcmp(nullcline_type, 'u')
            input = x;
            output = y;
            var_name = 'u';
        else
            input = y;
            output = x;
            var_name = 'v';
        end
        
        % 使用改进的buildTheta函数构建特征矩阵
        [Theta, term_descriptions] = buildTheta(input, poly_order, extra_basis, var_name);
        
        % 应用SINDy (带正则化的线性回归)
        xi = (Theta' * Theta + lambda * eye(size(Theta, 2))) \ (Theta' * output);
        
        % 创建符号表达式
        nullcline_expr = create_symbolic_expression(xi, term_descriptions);
    end

    function [Theta, term_descriptions] = buildTheta(input, poly_order, extra_basis, var_name)
        % 构建包含多项式项、额外基函数及其乘积的特征矩阵
        n = length(input);
        
        % === 步骤1: 创建基础特征集 ===
        base_features = {};
        base_descriptions = {};
        
        % 添加原始输入变量
        base_features{end+1} = input;
        base_descriptions{end+1} = var_name;
        
        % 添加额外基函数
        if nargin >= 3 && ~isempty(extra_basis)
            for i = 1:length(extra_basis)
                f = extra_basis{i};
                base_features{end+1} = f(input);
                
                % 获取函数名
                func_name = func2str(f);
                % 简化函数名显示
                if contains(func_name, '@')
                    func_name = ['f', num2str(i)];
                end
                base_descriptions{end+1} = [func_name, '(', var_name, ')'];
            end
        end
        
        num_base_features = length(base_features);
        
        % === 步骤2: 初始化Theta和项描述 ===
        Theta = ones(n, 1);  % 常数项
        term_descriptions = {'1'};  % 常数项描述
        
        % === 步骤3: 添加所有可能的乘积组合 ===
        for k = 1:poly_order
            % 修复: 当k大于基础特征数量时的处理
            if k > num_base_features
                % 当阶数大于基础特征数量时，只生成重复组合
                combinations = nchoosek_with_replacement(1:num_base_features, k);
            else
                % 生成不重复组合
                combinations = nchoosek(1:num_base_features, k);
                
                % 添加允许重复的组合
                repeated_combs = nchoosek_with_replacement(1:num_base_features, k);
                combinations = unique([combinations; repeated_combs], 'rows');
            end
            
            for j = 1:size(combinations, 1)
                term = ones(n, 1);
                term_desc_parts = {};
                
                % 计算该项的乘积
                for idx = 1:k
                    feature_idx = combinations(j, idx);
                    term = term .* base_features{feature_idx};
                    term_desc_parts{end+1} = base_descriptions{feature_idx};
                end
                
                % 添加到Theta矩阵
                Theta = [Theta, term];
                
                % 构建完整的项描述
                [unique_terms, ~, ic] = unique(term_desc_parts);
                term_desc = '';
                for m = 1:length(unique_terms)
                    count = sum(ic == m);
                    term_part = unique_terms{m};
                    
                    % 处理指数表示
                    if count > 1
                        term_part = [term_part, '^', num2str(count)];
                    end
                    
                    % 添加到完整描述
                    if m == 1
                        term_desc = term_part;
                    else
                        term_desc = [term_desc, '*', term_part];
                    end
                end
                
                term_descriptions{end+1} = term_desc;
            end
        end
        
        % 显示特征矩阵信息
        fprintf('构建的特征矩阵包含 %d 个特征项\n', size(Theta, 2));
    end

    function combinations = nchoosek_with_replacement(set, k)
        % 生成允许重复的组合（非递减顺序）
        n = length(set);
        combinations = [];
        
        % 处理特殊情况
        if n == 0 || k == 0
            return;
        end
        
        % 使用更健壮的迭代方法
        current = ones(1, k);
        combinations = zeros(0, k);
        
        while true
            % 添加当前组合
            combinations = [combinations; set(current)];
            
            % 找到下一个组合
            i = k;
            while i >= 1
                if current(i) < n
                    current(i) = current(i) + 1;
                    current(i+1:k) = current(i);
                    break;
                end
                i = i - 1;
            end
            
            if i < 1
                break;  % 所有组合已生成
            end
        end
    end

    function expr_str = create_symbolic_expression(xi, term_descriptions)
        % 根据系数和特征描述构建数学表达式
        expr_str = '';
        first_term = true;  % 标记是否是第一项
        
        % 按系数绝对值排序，从大到小
        [~, sorted_idx] = sort(abs(xi), 'descend');
        
        % 添加最重要的项（系数显著的非零项）
        for i = 1:length(sorted_idx)
            idx = sorted_idx(i);
            coeff = xi(idx);
            term_desc = term_descriptions{idx};
            
            % 忽略接近零的系数
            if abs(coeff) < 1e-4
                continue;
            end
            
            % === 处理项显示 ===
            if strcmp(term_desc, '1')
                % 常数项处理
                if abs(coeff - 1) < 1e-4
                    term_str = '1';
                elseif abs(coeff + 1) < 1e-4
                    term_str = '-1';
                else
                    term_str = sprintf('%.4f', coeff);
                end
            else
                % 非常数项处理
                if abs(abs(coeff) - 1) < 1e-4
                    % 系数为±1
                    if coeff > 0
                        term_str = term_desc;
                    else
                        term_str = ['-', term_desc];
                    end
                else
                    % 其他系数
                    term_str = sprintf('%.4f*%s', abs(coeff), term_desc);
                    if coeff < 0
                        term_str = ['-', term_str];
                    end
                end
            end
            
            % === 添加到表达式 ===
            if first_term
                % 第一项直接添加
                expr_str = term_str;
                first_term = false;
            else
                % 后续项添加符号
                if term_str(1) == '-'
                    expr_str = [expr_str, term_str];  % 负项直接连接
                else
                    expr_str = [expr_str, ' + ', term_str];  % 正项添加加号
                end
            end
        end
        
        % 处理全零情况
        if isempty(expr_str)
            expr_str = '0';
        end
        
        % 简化表达式：移除多余的"1*"和" + -"
        expr_str = strrep(expr_str, '1*', '');
        expr_str = strrep(expr_str, '+ -', '- ');
    end

    function visualize_results(u, v, u_null, v_null, u_derivs, v_derivs, u_expr, v_expr)
        % 创建新图形
        figure('Position', [100, 100, 1200, 500]);
        
        % 相空间图
        subplot(1, 2, 1);
        plot(u, v, 'b-', 'LineWidth', 1.5);
        hold on;
        
        % 绘制nullcline点
        if ~isempty(u_null)
            scatter(u_null(:,1), u_null(:,2), 30, 'r', 'filled');
        end
        if ~isempty(v_null)
            scatter(v_null(:,1), v_null(:,2), 30, 'g', 'filled');
        end
        
        % 生成数值范围
        u_vals = linspace(min(u), max(u), 200);
        v_vals = linspace(min(v), max(v), 200);
        
        % 初始化标志
        has_u_fit = false;
        has_v_fit = false;
        
        % 修复并求值u-nullcline (v = f(u))
        if ~isempty(u_expr)
            try
                % 修复表达式：添加缺失的乘法运算符
                u_expr_fixed = regexprep(u_expr, '(\d)([a-zA-Z])', '$1*$2');
                u_expr_fixed = regexprep(u_expr_fixed, '([a-zA-Z])(\d)', '$1*$2');
                u_expr_fixed = regexprep(u_expr_fixed, '([a-zA-Z])([a-zA-Z])', '$1*$2');
                
                % 创建函数句柄
                u_func = str2func(['@(u) ' u_expr_fixed]);
                
                % 安全评估表达式
                v_fit = arrayfun(@(u_val) safe_eval(u_func, u_val), u_vals);
                
                % 过滤无效值
                valid_idx = isfinite(v_fit) & (v_fit >= min(v)) & (v_fit <= max(v));
                if any(valid_idx)
                    plot(u_vals(valid_idx), v_fit(valid_idx), 'm--', 'LineWidth', 2);
                    has_u_fit = true;
                end
            catch ME
                warning('无法求值u-nullcline: %s\n修复后表达式: %s\n错误: %s', u_expr, u_expr_fixed, ME.message);
            end
        end
        
        % 修复并求值v-nullcline (u = g(v))
        if ~isempty(v_expr)
            try
                % 修复表达式：添加缺失的乘法运算符
                v_expr_fixed = regexprep(v_expr, '(\d)([a-zA-Z])', '$1*$2');
                v_expr_fixed = regexprep(v_expr_fixed, '([a-zA-Z])(\d)', '$1*$2');
                v_expr_fixed = regexprep(v_expr_fixed, '([a-zA-Z])([a-zA-Z])', '$1*$2');
                
                % 创建函数句柄
                v_func = str2func(['@(v) ' v_expr_fixed]);
                
                % 安全评估表达式
                u_fit = arrayfun(@(v_val) safe_eval(v_func, v_val), v_vals);
                
                % 过滤无效值
                valid_idx = isfinite(u_fit) & (u_fit >= min(u)) & (u_fit <= max(u));
                if any(valid_idx)
                    plot(u_fit(valid_idx), v_vals(valid_idx), 'c--', 'LineWidth', 2);
                    has_v_fit = true;
                end
            catch ME
                warning('无法求值v-nullcline: %s\n修复后表达式: %s\n错误: %s', v_expr, v_expr_fixed, ME.message);
            end
        end
        
        title('相空间与Nullcline拟合');
        xlabel('状态变量 1');
        ylabel('状态变量 2');
        legend_entries = {'轨迹'};
        if ~isempty(u_null), legend_entries{end+1} = 'u-nullcline点'; end
        if ~isempty(v_null), legend_entries{end+1} = 'v-nullcline点'; end
        if has_u_fit, legend_entries{end+1} = '拟合u-nullcline'; end
        if has_v_fit, legend_entries{end+1} = '拟合v-nullcline'; end
        legend(legend_entries, 'Location', 'best');
        grid on;
        
        % 导数分布图
        subplot(1, 2, 2);
        % u-nullcline点的|du/dt|
        if ~isempty(u_null)
            scatter(u_null(:,1), u_null(:,2), 30, log10(abs(u_derivs(:,1)) + 1e-10), 'filled');
            hold on;
        end
        % v-nullcline点的|dv/dt|
        if ~isempty(v_null)
            scatter(v_null(:,1), v_null(:,2), 30, log10(abs(v_derivs(:,2)) + 1e-10), 'filled');
        end
        colorbar;
        title('导数大小分布 (log_{10}尺度)');
        xlabel('状态变量 1');
        ylabel('状态变量 2');
        legend_entries = {};
        if ~isempty(u_null), legend_entries{end+1} = 'u-nullcline点 (log|du/dt|)'; end
        if ~isempty(v_null), legend_entries{end+1} = 'v-nullcline点 (log|dv/dt|)'; end
        legend(legend_entries, 'Location', 'best');
        grid on;
        
        % 添加表达式到图形
        annotation('textbox', [0.15, 0.02, 0.7, 0.05], 'String', ...
            sprintf('u-nullcline: v = %s   |   v-nullcline: u = %s', u_expr, v_expr), ...
            'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
            'FontSize', 10, 'EdgeColor', 'none');
    end

    function result = safe_eval(func, x)
        % 安全评估函数，避免错误
        try
            result = func(x);
        catch
            result = NaN;
        end
    end


end
