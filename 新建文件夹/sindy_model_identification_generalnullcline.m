function results = sindy_model_identification_generalnullcline(X, dt, dXdt, varargin)
% SINDY_MODEL_IDENTIFICATION 从动力系统数据识别模型并可视化结果
%   results = sindy_model_identification(X, dt)
%   results = sindy_model_identification(X, dt, 'Parameter', Value, ...)
%
% 输入参数:
%   X - 状态变量矩阵 [N×2]，每列是一个状态变量
%   dt - 时间步长
%   dXdt - 导数矩阵
%
% 可选参数:
%   'Offsets' - 偏移范围矩阵 [M×2]，默认值: [-2:0.1:2; -2:0.1:2]
%   'Threshold' - SINDy稀疏阈值，默认值: 0.25
%   'PolyOrder' - 多项式阶数，默认值: 3
%   'PlotResults' - 是否绘制结果，默认值: false
%   'URange' - u的绘图范围 [u_min, u_max]，默认自动确定
%   'VRange' - v的绘图范围 [v_min, v_max]，默认自动确定
%
% 输出:
%   results - 结构体数组，包含所有模型结果
%   top_models - 前10个最优模型的结构体数组

%% 参数解析
p = inputParser;
addRequired(p, 'X', @(x) size(x,2)==2);
addRequired(p, 'dt', @isscalar);
addParameter(p, 'Offsets', [-2:0.1:2; -2:0.1:2]', @ismatrix);
addParameter(p, 'Threshold', 0.08, @isscalar);
addParameter(p, 'PolyOrder', 3, @isscalar);
addParameter(p, 'PlotResults', false, @islogical);
addParameter(p, 'URange', [], @(x) isempty(x) || (isvector(x) && length(x)==2));
addParameter(p, 'VRange', [], @(x) isempty(x) || (isvector(x) && length(x)==2));
parse(p, X, dt, varargin{:});

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
    u_range = [min(u_data)-1, max(u_data)+1];
end
if isempty(v_range)
    v_range = [min(v_data)-1, max(v_data)+1];
end

%% 1. 计算导数 (中心差分)
%du_data = gradient(u_data, dt);
du_data = dXdt(:,1);
dv_data = dXdt(:,2);
%dv_data = gradient(v_data, dt);

%% 2. 设置偏移范围
num_models = length(du_offsets)*length(dv_offsets);

% 预存储结果
results = struct('du_offset', cell(num_models,1), 'dv_offset', cell(num_models,1), ...
                 'k', cell(num_models,1), 'R2', cell(num_models,1), ...
                 'Xi', cell(num_models,1));

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
        Theta = poolData([u_shifted, v_shifted], 2, 3);
        
        
        % 应用 SINDy 
        dXdt = [du_data, dv_data]; % 导数矩阵
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
fprintf('排名 | Δu   | Δv   | 复杂度k | R²得分 | 综合得分\n');
fprintf('----|------|------|---------|--------|--------\n');

for i = 1:length(top_models)
    model = top_models(i);
    
    fprintf('%2d  | %4.1f | %4.1f | %7d | %6.4f | %7.4f\n', ...
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
    
    fprintf('%2d  | du/dt = ', i);
    print_equation(Xi(:,1), term_names);
    
    fprintf('\n    | dv/dt = ');
    print_equation(Xi(:,2), term_names);
    fprintf('\n----|-------------------------------------------------\n');
end

%% 7. 可视化结果（可选）
if plot_results
    % 调用融合后的可视化函数
    plot_top_models_combined(top_models, all_offsets, all_K, u_range, v_range, u_data, v_data);
end

end

%% 辅助函数

function term_names = generate_term_names(poly_order)
    % 生成项名称
    term_names = {'const', 'u', 'v'};
    
    if poly_order >= 2
        term_names = [term_names, {'u^2', 'uv', 'v^2'}];
    end
    
    if poly_order >= 3
        term_names = [term_names, {'u^3', 'u^2v', 'uv^2', 'v^3'}];
    end
    
    if poly_order >= 4
        term_names = [term_names, {'u^4', 'u^3v', 'u^2v^2', 'uv^3', 'v^4'}];
    end
    %if poly_order >=5
     %   term_names = [term_names, {'u^5', 'u^4v', 'u^3v^2', 'u^2v^3', 'uv^4', 'v^5'}];
    %end
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
    
    % 绘制K值等高线图
    contourf(DU, DV, K_grid, 20, 'LineStyle', 'none');
    colormap(parula);
    colorbar;
    title('K值分布图 (模型复杂度)', 'FontSize', 12);
    xlabel('\Delta u', 'FontSize', 10);
    ylabel('\Delta v', 'FontSize', 10);
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
    h_u = gobjects(length(top_models), 1);
    h_v = gobjects(length(top_models), 1);
    h_labels = gobjects(length(top_models), 2);
    
    % 绘制nullcline
    for i = 1:min(10, length(top_models))
        model = top_models(i);
        Xi = model.Xi;
        
        % 计算nullcline函数值
        F = Xi(1,1) + Xi(2,1)*U + Xi(3,1)*V;
        G = Xi(1,2) + Xi(2,2)*U + Xi(3,2)*V;
        
        % 添加高阶项
        if size(Xi,1) >= 4
            F = F + Xi(4,1)*U.^2 + Xi(5,1)*U.*V + Xi(6,1)*V.^2;
            G = G + Xi(4,2)*U.^2 + Xi(5,2)*U.*V + Xi(6,2)*V.^2;
        end
        if size(Xi,1) >= 7
            F = F + Xi(7,1)*U.^3 + Xi(8,1)*U.^2.*V + Xi(9,1)*U.*V.^2 + Xi(10,1)*V.^3;
            G = G + Xi(7,2)*U.^3 + Xi(8,2)*U.^2.*V + Xi(9,2)*U.*V.^2 + Xi(10,2)*V.^3;
        end
        
        % 绘制nullcline并存储句柄
        [Cu, h_u(i)] = contour(U, V, F, [0 0], 'LineWidth', 1.5, 'LineColor', colors(i,:));
        hold on;
        [Cv, h_v(i)] = contour(U, V, G, [0 0], 'LineWidth', 1.5, 'LineColor', colors(i,:), 'LineStyle', '--');
        
        % 添加模型编号标签
        if ~isempty(Cu)
            % 找到曲线上的一个点放置标签
            k = 1;
            while k < size(Cu, 2)
                n_points = Cu(2, k);
                if n_points > 10
                    mid_idx = k + floor(n_points/2);
                    x_pos = Cu(1, mid_idx);
                    y_pos = Cu(2, mid_idx);
                    
                    % 放置u-nullcline标签
                    h_labels(i, 1) = text(x_pos, y_pos, sprintf('M%d: du/dt=0', i),...
                                          'BackgroundColor', 'white',...
                                          'FontSize', 9,...
                                          'Color', colors(i,:),...
                                          'HorizontalAlignment', 'center');
                    break;
                end
                k = k + n_points + 1;
            end
        end
        
        if ~isempty(Cv)
            % 找到曲线上的一个点放置标签
            k = 1;
            while k < size(Cv, 2)
                n_points = Cv(2, k);
                if n_points > 10
                    mid_idx = k + floor(n_points/2);
                    x_pos = Cv(1, mid_idx);
                    y_pos = Cv(2, mid_idx);
                    
                    % 放置v-nullcline标签
                    h_labels(i, 2) = text(x_pos, y_pos, sprintf('M%d: dv/dt=0', i),...
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
    xlabel('u', 'FontSize', 10);
    ylabel('v', 'FontSize', 10);
    grid on;
    
    % 应用统一坐标范围
    xlim(common_range);
    ylim(common_range);
    
    % 添加图例说明
    text(0.05, 0.95, '实线: du/dt=0 (u-nullcline)', 'Units', 'normalized', 'FontSize', 9);
    text(0.05, 0.90, '虚线: dv/dt=0 (v-nullcline)', 'Units', 'normalized', 'FontSize', 9);
    
    %% 4. 添加模型图例 (修复透明度问题)
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
    sgtitle('动力系统模型识别与可视化综合分析', 'FontSize', 16, 'FontWeight', 'bold');
end