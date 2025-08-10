%% SCR动力系统建模与分析框架
clear all; close all; clc;

%% ========== SCR时间序列数据生成 ==========
% 加载原始数据
data = [
    2537, 1.843, 228, 2527.4, 0.37978, 0;
    2765, 1.8148, 222, 2757.6, 0.26879, -0.01055;
    817, 0.64362, 264, 804.02, 1.5888, -0.0058432;
    2784, 2.1118, 268, 2735.1, 1.7579, 0.024825;
    2081, 1.5056, 200, 2078.6, 0.11479, 0.051704;
    743, 0.62357, 234, 738.86, 0.55674, -0.0039582;
    1196, 0.90989, 229, 1191, 0.41398, -0.032749;
    1867, 1.2766, 251, 1850.8, 0.86606, -0.014275;
    2894, 2.0552, 258, 2860.4, 1.1593, -0.072172;
    2913, 1.9604, 263, 2874.8, 1.312, 0.024152;
    894, 0.71464, 210, 892.24, 0.19721, 0.16401;
    2927, 1.9342, 255, 2898.6, 0.9687, -0.0046954;
    2893, 1.9182, 252, 2867.8, 0.87136, 0.016995;
    1713, 1.1182, 198, 1711.4, 0.094606, 0.014769;
    2501, 1.4563, 193, 2499.3, 0.06708, -0.0025569;
    854, 0.74164, 235, 848.88, 0.59906, -0.00072129;
    1554, 1.1203, 286, 1506.8, 3.04, 0.033281;
    2790, 1.8775, 217, 2783.8, 0.228, 0.072381;
    2481, 1.8737, 244, 2462.7, 0.73665, 0.017215;
    2899, 1.8581, 204, 2895.5, 0.12161, 0
];

% 提取变量
NOx_in = data(:,1);   % 入口NOx浓度 (ppm)
urea = data(:,2);     % 尿素喷射量 (kg/h)
T_scr = data(:,3);    % SCR温度 (°C)
NOx_out = data(:,4);  % 出口NOx浓度 (ppm)
eff = data(:,5);      % 脱硝效率 (%)
dEff_dT = data(:,6);  % 效率变化率 (%/°C)

% 时间序列转换
dt = 0.5; % 采样时间间隔（分钟）
t = (0:dt:(length(NOx_in)-1)*dt)'; % 时间向量

% 添加动态特性
T_scr_ts = zeros(size(T_scr));   % SCR温度时间序列
eff_ts = zeros(size(eff));       % 效率时间序列

for i = 1:length(T_scr)
    if i == 1
        T_scr_ts(i) = T_scr(i);
        eff_ts(i) = eff(i);
    else
        % 温度动态（热惯性）
        tau_temp = 3.0; % 温度响应时间常数（分钟）
        T_scr_ts(i) = T_scr_ts(i-1) + (T_scr(i) - T_scr_ts(i-1)) * (dt/tau_temp);
        
        % 效率动态（化学反应动力学）
        k_react = 0.05 * exp(-8000/(8.314*(T_scr_ts(i)+273.15))); % 反应速率常数
        tau_eff = 1/(k_react + 0.01); % 效率响应时间常数
        eff_ts(i) = eff_ts(i-1) + (eff(i) - eff_ts(i-1)) * (dt/tau_eff);
    end
end

%% ========== 关键状态变量提取 ==========
% 选择核心状态变量: SCR温度(T_scr)和脱硝效率(eff)
X = [T_scr_ts, eff_ts];  % 状态矩阵 [温度, 效率]

% 计算状态导数
dT_dt = gradient(X(:,1), dt);      % 温度变化率
deff_dt = gradient(X(:,2), dt);    % 效率变化率
dXdt = [dT_dt, deff_dt];           % 导数矩阵

%% ========== SINDy模型识别 ==========
% 设置偏移范围 - 针对SCR系统的物理特性优化
temp_offsets = linspace(-10, 10, 11);    % 温度偏移范围 ±10°C
eff_offsets = linspace(-5, 5, 11);       % 效率偏移范围 ±5%
offsets = [temp_offsets; eff_offsets]';  % 组合偏移矩阵

% 运行SINDy模型识别
results = sindy_model_identification_scr(X, dt, dXdt, ...
    'Offsets', offsets, ...
    'PolyOrder', 3, ...
    'Threshold', 0.1, ...
    'PlotResults', true, ...
    'URange', [min(X(:,1))-5, max(X(:,1))+5], ...
    'VRange', [min(X(:,2))-5, max(X(:,2))+5]);

%% ========== SCR系统不动点分析 ==========
% 提取最优模型
if ~isempty(results.top_models)
    top_model = results.top_models(1);
    Xi = top_model.Xi;
    poly_order = top_model.PolyOrder;
    
    % 定义符号变量
    syms u v
    vars = [u, v];
    
    % 构建符号方程
    Theta_sym = poolData_sym(vars, 2, poly_order);
    
    % 检查维度匹配
    if length(Theta_sym) == size(Xi, 1)
        du_dt = Theta_sym .* Xi(:,1);
        dv_dt = Theta_sym .* Xi(:,2);
        
        % 求解不动点
        fixed_points = vpasolve([du_dt == 0, dv_dt == 0], [u, v], ...
            [min(X(:,1)), max(X(:,1)); min(X(:,2)), max(X(:,2))]);
        
        % 显示不动点
        fprintf('\n===== SCR系统不动点分析 =====\n');
        if ~isempty(fixed_points.u)
            for i = 1:length(fixed_points.u)
                T_fp = double(fixed_points.u(i));
                eff_fp = double(fixed_points.v(i));
                
                % 计算雅可比矩阵
                J = jacobian([du_dt, dv_dt], [u, v]);
                J_fp = double(subs(J, [u, v], [T_fp, eff_fp]));
                
                % 计算特征值
                eigvals = eig(J_fp);
                
                % 稳定性判断
                if all(real(eigvals) < 0)
                    stability = '稳定';
                elseif all(real(eigvals) > 0)
                    stability = '不稳定';
                else
                    stability = '鞍点';
                end
                
                fprintf('不动点 %d: T = %.2f°C, η = %.2f%%, 类型: %s\n', ...
                    i, T_fp, eff_fp, stability);
                fprintf('特征值: %.4f %+.4fi, %.4f %+.4fi\n', ...
                    real(eigvals(1)), imag(eigvals(1)), real(eigvals(2)), imag(eigvals(2)));
            end
        else
            fprintf('未找到不动点\n');
        end
    else
        fprintf('\n===== 维度不匹配，跳过不动点分析 =====\n');
        fprintf('符号库项数: %d, 系数向量长度: %d\n', length(Theta_sym), size(Xi,1));
        fixed_points = [];
    end
else
    fprintf('\n===== 无有效模型，跳过不动点分析 =====\n');
    fixed_points = [];
end

%% ========== SCR系统相空间分析 ==========
if ~isempty(results.top_models)
    figure('Position', [100, 100, 1200, 900]);
    
    % 创建相空间网格
    u_min = min(X(:,1)) - 5; u_max = max(X(:,1)) + 5;
    v_min = min(X(:,2)) - 5; v_max = max(X(:,2)) + 5;
    [U, V] = meshgrid(linspace(u_min, u_max, 30), linspace(v_min, v_max, 30));
    
    % 计算向量场
    DU = zeros(size(U));
    DV = zeros(size(V));
    for i = 1:size(U,1)
        for j = 1:size(U,2)
            Theta_point = poolData([U(i,j), V(i,j)], 2, poly_order);
            % 确保维度匹配
            if size(Theta_point, 2) == size(Xi, 1)
                dXdt_point = Theta_point * Xi;
                DU(i,j) = dXdt_point(1);
                DV(i,j) = dXdt_point(2);
            else
                DU(i,j) = 0;
                DV(i,j) = 0;
            end
        end
    end
    
    % 归一化向量长度
    magnitude = sqrt(DU.^2 + DV.^2);
    DU_norm = DU ./ (magnitude + 1e-3);
    DV_norm = DV ./ (magnitude + 1e-3);
    
    % 绘制相空间轨迹
    subplot(2,2,1);
    streamslice(U, V, DU_norm, DV_norm, 1.5);
    hold on;
    plot(X(:,1), X(:,2), 'k-', 'LineWidth', 1.5); % 实际轨迹
    plot(X(1,1), X(1,2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); % 起点
    plot(X(end,1), X(end,2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % 终点
    
    % 标记不动点
    if ~isempty(fixed_points) && ~isempty(fixed_points.u)
        for i = 1:length(fixed_points.u)
            T_fp = double(fixed_points.u(i));
            eff_fp = double(fixed_points.v(i));
            plot(T_fp, eff_fp, 'ko', 'MarkerSize', 12, 'LineWidth', 2);
        end
    end
    
    title('SCR系统相空间');
    xlabel('温度 (°C)');
    ylabel('脱硝效率 (%)');
    grid on;
    xlim([u_min, u_max]);
    ylim([v_min, v_max]);
    
    % 绘制向量场
    subplot(2,2,2);
    quiver(U, V, DU, DV, 1.5, 'k');
    hold on;
    contour(U, V, magnitude, 20, 'LineWidth', 1.5);
    colorbar;
    title('向量场与变化率强度');
    xlabel('温度 (°C)');
    ylabel('脱硝效率 (%)');
    grid on;
    
    % 绘制温度-效率关系
    subplot(2,2,3);
    plot(X(:,1), X(:,2), 'b-', 'LineWidth', 1.5);
    hold on;
    scatter(X(:,1), X(:,2), 30, 1:length(X(:,1)), 'filled');
    colorbar;
    title('温度-效率轨迹');
    xlabel('温度 (°C)');
    ylabel('脱硝效率 (%)');
    grid on;
    
    % 标记关键操作点
    [eff_max, idx_max] = max(X(:,2));
    plot(X(idx_max,1), eff_max, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(X(idx_max,1)+1, eff_max, sprintf('最大效率: %.1f%%', eff_max));
    
    % 绘制时间演化
    subplot(2,2,4);
    yyaxis left;
    plot(t, X(:,1), 'b-', 'LineWidth', 1.5);
    ylabel('温度 (°C)');
    yyaxis right;
    plot(t, X(:,2), 'r-', 'LineWidth', 1.5);
    ylabel('脱硝效率 (%)');
    title('时间演化');
    xlabel('时间 (分钟)');
    grid on;
    legend('温度', '效率');
end

%% ========== 辅助函数 ==========
function Theta_sym = poolData_sym(vars, n_vars, poly_order)
    % 符号多项式库
    u = vars(1);
    v = vars(2);
    
    % 常数项
    Theta_sym = [1];
    
    % 线性项
    Theta_sym = [Theta_sym, u, v];
    
    % 二阶项
    if poly_order >= 2
        Theta_sym = [Theta_sym, u^2, u*v, v^2];
    end
    
    % 三阶项
    if poly_order >= 3
        Theta_sym = [Theta_sym, u^3, u^2*v, u*v^2, v^3];
    end
    
    % 四阶项
    if poly_order >= 4
        Theta_sym = [Theta_sym, u^4, u^3*v, u^2*v^2, u*v^3, v^4];
    end
    
    Theta_sym = Theta_sym.';
end

%% ========== SINDy模型识别函数 ==========
function output = sindy_model_identification_scr(X, dt, dXdt, varargin)
    % 参数解析
    p = inputParser;
    addRequired(p, 'X', @(x) size(x,2)==2);
    addRequired(p, 'dt', @isscalar);
    addRequired(p, 'dXdt', @ismatrix);
    addParameter(p, 'Offsets', [-10:1:10; -5:0.5:5]', @ismatrix);
    addParameter(p, 'Threshold', 0.1, @isscalar);
    addParameter(p, 'PolyOrder', 3, @isscalar);
    addParameter(p, 'PlotResults', false, @islogical);
    addParameter(p, 'URange', [], @(x) isempty(x) || (isvector(x) && length(x)==2));
    addParameter(p, 'VRange', [], @(x) isempty(x) || (isvector(x) && length(x)==2));
    parse(p, X, dt, dXdt, varargin{:});
    
    u_data = X(:,1);
    v_data = X(:,2);
    du_offsets = p.Results.Offsets(:,1)';
    dv_offsets = p.Results.Offsets(:,2)';
    xi_thres = p.Results.Threshold;
    poly_order = p.Results.PolyOrder;
    plot_results = p.Results.PlotResults;
    u_range = p.Results.URange;
    v_range = p.Results.VRange;
    
    % 自动确定绘图范围
    if isempty(u_range)
        u_range = [min(u_data)-5, max(u_data)+5];
    end
    if isempty(v_range)
        v_range = [min(v_data)-5, max(v_data)+5];
    end
    
    %% 1. 计算导数 (使用输入值)
    du_data = dXdt(:,1);
    dv_data = dXdt(:,2);
    
    %% 2. 设置偏移范围
    num_models = length(du_offsets)*length(dv_offsets);
    
    % 预存储结果
    results = struct('du_offset', cell(num_models,1), 'dv_offset', cell(num_models,1), ...
                     'k', cell(num_models,1), 'R2', cell(num_models,1), ...
                     'Xi', cell(num_models,1), 'PolyOrder', poly_order);
    
    % 存储所有偏移量和K值用于绘图
    all_offsets = zeros(num_models, 2);
    all_K = zeros(num_models, 1);
    
    %% 3. 主循环：应用不同偏移
    model_idx = 1;
    for i = 1:length(du_offsets)
        for j = 1:length(dv_offsets)
            % 应用偏移
            du_offset = du_offsets(i);
            dv_offset = dv_offsets(j);
            u_shifted = u_data + du_offset;
            v_shifted = v_data + dv_offset;
            
            % 构建多项式库
            Theta = poolData([u_shifted, v_shifted], 2, poly_order);
            
            % 应用 SINDy 
            Xi = sparsifyDynamics(Theta, dXdt, xi_thres, 2); 
    
            % 计算模型复杂度 k (非零项数)
            k = nnz(Xi);
            
            % 计算 R² 分数
            pred = Theta * Xi;
            R2_u = 1 - sum((du_data - pred(:,1)).^2) / sum((du_data - mean(du_data)).^2);
            R2_v = 1 - sum((dv_data - pred(:,2)).^2) / sum((dv_data - mean(dv_data)).^2);
            R2 = mean([R2_u, R2_v]);
            
            % 存储结果
            results(model_idx).du_offset = du_offset;
            results(model_idx).dv_offset = dv_offset;
            results(model_idx).k = k;
            results(model_idx).R2 = R2;
            results(model_idx).Xi = Xi;
            results(model_idx).PolyOrder = poly_order;
            
            % 存储绘图数据
            all_offsets(model_idx,:) = [du_offset, dv_offset];
            all_K(model_idx) = k;
            
            model_idx = model_idx + 1;
        end
    end
    
    %% 4. 模型评估与选择
    k_simp = 1; 
    k_max = size(Theta, 2); % 基于库大小
    
    k_values = [results.k];
    R2_values = [results.R2];
    
    % 归一化指标
    k_norm = (k_values - k_simp) / (k_max - k_simp);
    R2_norm = 1 - R2_values;           % 转换为误差
    
    % 计算综合指标 score = k_norm + R2_norm
    scores = k_norm + R2_norm;
    
    % 将得分添加到结果结构体
    for i = 1:length(results)
        results(i).score = scores(i);
    end
    
    [sorted_scores, sorted_idx] = sort(scores);
    top_models = results(sorted_idx(1:min(10, numel(results)))); % 选择前10个最佳模型
    
    % 为top_models添加offset字段用于绘图
    for i = 1:length(top_models)
        top_models(i).offset = [top_models(i).du_offset, top_models(i).dv_offset];
        top_models(i).K = top_models(i).k;
    end
    
    %% 5. 打印最优模型信息
    fprintf('\n===== 前10个最优模型 =====\n');
    fprintf('排名 | ΔT(°C) | Δη(%%) | 复杂度k | R²得分 | 综合得分\n');
    fprintf('----|--------|------|---------|--------|--------\n');
    
    for i = 1:length(top_models)
        model = top_models(i);
        
        fprintf('%2d  | %6.1f | %5.1f | %7d | %6.4f | %7.4f\n', ...
                i, model.du_offset, model.dv_offset, model.k, model.R2, model.score);
    end
    
    %% 6. 打印最优模型方程
    fprintf('\n===== 最优模型方程 =====\n');
    fprintf('排名 | 方程\n');
    fprintf('----|-------------------------------------------------\n');
    
    term_names = generate_term_names(poly_order);
    
    for i = 1:length(top_models)
        model = top_models(i);
        Xi = model.Xi;
        
        fprintf('%2d  | dT/dt = ', i);
        print_equation(Xi(:,1), term_names);
        
        fprintf('\n    | dη/dt = ');
        print_equation(Xi(:,2), term_names);
        fprintf('\n----|-------------------------------------------------\n');
    end
    
    %% 7. 可视化结果
    if plot_results && ~isempty(top_models)
        plot_top_models_combined(top_models, all_offsets, all_K, u_range, v_range, u_data, v_data);
    end
    
    %% 8. 创建输出结构体
    output = struct();
    output.all_models = results;
    output.top_models = top_models;
    output.all_offsets = all_offsets;
    output.all_K = all_K;
end

function term_names = generate_term_names(poly_order)
    % 生成项名称
    term_names = {'const', 'T', 'η'};
    
    if poly_order >= 2
        term_names = [term_names, {'T^2', 'Tη', 'η^2'}];
    end
    
    if poly_order >= 3
        term_names = [term_names, {'T^3', 'T^2η', 'Tη^2', 'η^3'}];
    end
    
    if poly_order >= 4
        term_names = [term_names, {'T^4', 'T^3η', 'T^2η^2', 'Tη^3', 'η^4'}];
    end
end

function print_equation(coeffs, term_names)
    % 打印方程
    printed = false;
    for j = 1:length(coeffs)
        if abs(coeffs(j)) > 1e-3
            if printed
                if coeffs(j) > 0
                    fprintf(' + ');
                else
                    fprintf(' - ');
                end
                fprintf('%5.3f*%s', abs(coeffs(j)), term_names{j});
            else
                printed = true;
                if coeffs(j) < 0
                    fprintf('-');
                end
                fprintf('%5.3f*%s', abs(coeffs(j)), term_names{j});
            end
        end
    end
    if ~printed
        fprintf('0');
    end
end

function plot_top_models_combined(top_models, all_offsets, all_K, u_range, v_range, u_data, v_data)
    % 融合后的可视化函数，结合了原始轨迹和nullcline分析
    
    %% 1. 创建综合可视化图形
    figure('Position', [100, 100, 1400, 600]);
    
    %% 2. K值分布与最优模型位置 (左)
    subplot(1,2,1);
    
    % 准备K值网格数据
    du_vals = unique(all_offsets(:,1));
    dv_vals = unique(all_offsets(:,2));
    [DU, DV] = meshgrid(du_vals, dv_vals);
    K_grid = reshape(all_K, length(dv_vals), length(du_vals));
    
    % 计算统一的坐标范围
    all_vals = [du_vals; dv_vals; u_range(:); v_range(:)];
    common_min = min(all_vals);
    common_max = max(all_vals);
    common_range = [common_min, common_max];
    
    % 检查K_grid是否为常数
    if all(K_grid(:) == K_grid(1))
        % 如果K值恒定，绘制热图而不是等高线
        imagesc(du_vals, dv_vals, K_grid);
        axis xy;
        colorbar;
        title('K值分布 (常数)', 'FontSize', 12);
    else
        % 绘制K值等高线图
        contourf(DU, DV, K_grid, 20, 'LineStyle', 'none');
        colormap(parula);
        colorbar;
        title('K值分布图 (模型复杂度)', 'FontSize', 12);
    end
    
    xlabel('\Delta T', 'FontSize', 10);
    ylabel('\Delta \eta', 'FontSize', 10);
    hold on;
    
    % 标记前10个最优模型位置
    top_offsets = cell2mat({top_models.offset}');
    scatter(top_offsets(:,1), top_offsets(:,2), 100, 'r', 'filled');
    
    % 添加不动点标记 (0,0)
    plot(0, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    
    % 添加图例
    legend('K值分布', '最优模型位置', '真实不动点', 'Location', 'best');
    grid on;
    
    % 应用统一坐标范围
    xlim(common_range);
    ylim(common_range);
    
    %% 3. 改进的Nullcline可视化 (右)
    subplot(1,2,2);
    
    % 创建相空间网格 (使用统一范围)
    [U, V] = meshgrid(linspace(common_min, common_max, 100),...
                      linspace(common_min, common_max, 100));
    
    % 使用分面显示方法
    hold on;
    
    % 获取10种鲜艳的区分度高的颜色
    colors = lines(10);
    
    % 存储所有绘图对象的句柄
    h_labels = gobjects(length(top_models), 2);
    
    % 绘制nullcline
    for i = 1:min(10, length(top_models))
        model = top_models(i);
        Xi = model.Xi;
        poly_order = model.PolyOrder;
        
        % 计算nullcline函数值
        F = Xi(1,1) + Xi(2,1)*U + Xi(3,1)*V;
        G = Xi(1,2) + Xi(2,2)*U + Xi(3,2)*V;
        
        % 添加高阶项
        if poly_order >= 2 && size(Xi,1) >= 6
            F = F + Xi(4,1)*U.^2 + Xi(5,1)*U.*V + Xi(6,1)*V.^2;
            G = G + Xi(4,2)*U.^2 + Xi(5,2)*U.*V + Xi(6,2)*V.^2;
        end
        if poly_order >= 3 && size(Xi,1) >= 10
            F = F + Xi(7,1)*U.^3 + Xi(8,1)*U.^2.*V + Xi(9,1)*U.*V.^2 + Xi(10,1)*V.^3;
            G = G + Xi(7,2)*U.^3 + Xi(8,2)*U.^2.*V + Xi(9,2)*U.*V.^2 + Xi(10,2)*V.^3;
        end
        
        % 绘制nullcline
        if any(F(:) ~= F(1)) % 检查是否非常数
            [Cu, ~] = contour(U, V, F, [0 0], 'LineWidth', 1.5, 'LineColor', colors(i,:));
        end
        if any(G(:) ~= G(1)) % 检查是否非常数
            [Cv, ~] = contour(U, V, G, [0 0], 'LineWidth', 1.5, 'LineColor', colors(i,:), 'LineStyle', '--');
        end
        
        % 添加模型编号标签
        if exist('Cu', 'var') && ~isempty(Cu)
            % 找到曲线上的一个点放置标签
            k = 1;
            while k < size(Cu, 2)
                n_points = Cu(2, k);
                if n_points > 10
                    mid_idx = k + floor(n_points/2);
                    x_pos = Cu(1, mid_idx);
                    y_pos = Cu(2, mid_idx);
                    
                    % 放置u-nullcline标签
                    h_labels(i, 1) = text(x_pos, y_pos, sprintf('M%d: dT/dt=0', i),...
                                          'BackgroundColor', 'white',...
                                          'FontSize', 9,...
                                          'Color', colors(i,:),...
                                          'HorizontalAlignment', 'center');
                    break;
                end
                k = k + n_points + 1;
            end
        end
        
        if exist('Cv', 'var') && ~isempty(Cv)
            % 找到曲线上的一个点放置标签
            k = 1;
            while k < size(Cv, 2)
                n_points = Cv(2, k);
                if n_points > 10
                    mid_idx = k + floor(n_points/2);
                    x_pos = Cv(1, mid_idx);
                    y_pos = Cv(2, mid_idx);
                    
                    % 放置v-nullcline标签
                    h_labels(i, 2) = text(x_pos, y_pos, sprintf('M%d: dη/dt=0', i),...
                                          'BackgroundColor', 'white',...
                                          'FontSize', 9,...
                                          'Color', colors(i,:),...
                                          'HorizontalAlignment', 'center');
                    break;
                end
                k = k + n_points + 1;
            end
        end
    end
    
    % 添加不动点标记 (0,0)
    plot(0, 0, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
    
    % 设置图形属性
    title('前10个最优模型的Nullcline', 'FontSize', 12);
    xlabel('T (°C)', 'FontSize', 10);
    ylabel('\eta (%)', 'FontSize', 10);
    grid on;
    
    % 应用统一坐标范围
    xlim(common_range);
    ylim(common_range);
    
    % 添加图例说明
    text(0.05, 0.95, '实线: dT/dt=0 (温度nullcline)', 'Units', 'normalized', 'FontSize', 9);
    text(0.05, 0.90, '虚线: dη/dt=0 (效率nullcline)', 'Units', 'normalized', 'FontSize', 9);
    
    %% 4. 添加模型图例
    % 创建图例项
    legend_items = gobjects(min(10, length(top_models)), 1);
    legend_labels = cell(min(10, length(top_models)), 1);
    
    for i = 1:min(10, length(top_models))
        % 创建不可见的线条用于图例
        legend_items(i) = plot(NaN, NaN, 'LineWidth', 2, 'Color', colors(i,:));
        legend_labels{i} = sprintf('模型 %d', i);
    end
    
    % 添加图例
    legend(legend_items, legend_labels, 'Location', 'southoutside', 'NumColumns', 5);
    
    %% 添加标题说明
    sgtitle('SCR系统动力系统模型识别与可视化综合分析', 'FontSize', 16, 'FontWeight', 'bold');
end

function Theta = poolData(yin, nVars, polyorder)
    % 创建多项式特征库
    n = size(yin,1);
    
    % 常数项
    Theta = ones(n,1);
    
    % 线性项
    Theta = [Theta, yin];
    
    % 二阶项
    if polyorder >= 2
        for i = 1:nVars
            for j = i:nVars
                Theta = [Theta, yin(:,i).*yin(:,j)];
            end
        end
    end
    
    % 三阶项
    if polyorder >= 3
        for i = 1:nVars
            for j = i:nVars
                for k = j:nVars
                    Theta = [Theta, yin(:,i).*yin(:,j).*yin(:,k)];
                end
            end
        end
    end
    
    % 四阶项
    if polyorder >= 4
        for i = 1:nVars
            for j = i:nVars
                for k = j:nVars
                    for l = k:nVars
                        Theta = [Theta, yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l)];
                    end
                end
            end
        end
    end
end

function Xi = sparsifyDynamics(Theta, dXdt, lambda, n)
    % 稀疏识别动力学
    Xi = zeros(size(Theta,2), n);
    for i = 1:n
        % 使用最小二乘法作为备用
        Xi(:,i) = pinv(Theta) * dXdt(:,i);
        
        % 应用Lasso正则化
        if lambda > 0
            [B, ~] = lasso(Theta, dXdt(:,i), 'Lambda', lambda);
            Xi(:,i) = B;
        end
    end
end