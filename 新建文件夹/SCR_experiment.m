%% SCR系统脱硝性能数据生成与零线重建分析
clear all; close all; clc;

%% 模型参数设置
k1 = 0.001;     % 主反应速率常数
Ea = 80000;     % 活化能 (J/mol)
R = 8.314;      % 气体常数 (J/(mol·K))
T_ref = 500;    % 参考温度 (K)

%% 生成实验数据
% 增加数据点数量以提高精度
n_samples = 50;
NOx_in = randi([500, 3000], n_samples, 1);
urea_flow = NOx_in * 0.0007 + randn(size(NOx_in)) * 0.1;
urea_flow = max(0.5, min(2.2, urea_flow));
T = randi([180, 290], n_samples, 1) + 273.15;

% 初始化结果数组
NOx_out = zeros(n_samples, 1);
efficiency = zeros(n_samples, 1);
dEff_dT = zeros(n_samples, 1);

%% 计算脱硝效率
for i = 1:n_samples
    k_T = k1 * exp(-Ea/R * (1/T(i) - 1/T_ref));
    efficiency(i) = 100 * (1 - exp(-k_T * urea_flow(i) * T(i) / (NOx_in(i) * 0.1)));
    efficiency(i) = max(0, min(100, efficiency(i)));
    NOx_out(i) = NOx_in(i) * (1 - efficiency(i)/100);
    
    % 计算导数
    if i > 1 && i < n_samples
        dEff_dT(i) = (efficiency(i+1) - efficiency(i-1)) / (T(i+1) - T(i-1));
    end
end

% 转换为摄氏度
T_C = T - 273.15;

%% 零线重建方法在SCR系统中的应用
% 定义状态变量
x = T_C;  % 温度 (°C)
y = NOx_in; % 入口NOx浓度 (ppm)
z = efficiency; % 脱硝效率 (%)

% 定义变换向量范围
delta_x_range = linspace(-10, 10, 15); % 温度扰动范围 ±10°C
delta_y_range = linspace(-200, 200, 15); % NOx扰动范围 ±200ppm

% 初始化K值矩阵
K_matrix = zeros(length(delta_x_range), length(delta_y_range));
R2_matrix = zeros(size(K_matrix));

%% 计算K值分布 (PDF方法核心)
for i = 1:length(delta_x_range)
    for j = 1:length(delta_y_range)
        % 施加变换向量
        x_shifted = x + delta_x_range(i);
        y_shifted = y + delta_y_range(j);
        
        % 使用SINDy识别模型 (SCR专用)
        [~, K, R2] = sindy_scr_identification(x_shifted, y_shifted, z);
        
        % 存储结果
        K_matrix(i, j) = K;
        R2_matrix(i, j) = R2;
    end
end

%% 识别零线位置 (K值下降点)
% 计算K值阈值 (PDF方法的关键观察)
K_threshold = prctile(K_matrix(:), 30); % 取K值最低的30%作为阈值

% 找出低K值区域 (潜在零线位置)
[nullcline_x, nullcline_y] = find(K_matrix < K_threshold);

% 转换为实际坐标
nullcline_points = [delta_x_range(nullcline_x)', delta_y_range(nullcline_y)'];

% 计算绝对位置
absolute_nullcline = [nullcline_points(:,1) + mean(x), ...
                      nullcline_points(:,2) + mean(y)];

%% 找出真实导数为零的点
threshold = 0.05;
zero_dT_idx = find(abs(dEff_dT) < threshold);

%% 可视化结果
figure('Position', [100, 100, 1400, 1000]);

% 子图1: K值分布与零线位置
subplot(2, 2, 1);
contourf(delta_y_range, delta_x_range, K_matrix, 20, 'LineColor', 'none');
colorbar;
hold on;
scatter(nullcline_points(:, 2), nullcline_points(:, 1), 50, 'r', 'filled');
title('K值分布与零线位置 (PDF方法)');
xlabel('ΔNOx (ppm)');
ylabel('ΔT (°C)');
legend('K值分布', '零线候选点');

% 子图2: 真实效率变化率与零线
subplot(2, 2, 2);
scatter(y, x, 100, dEff_dT, 'filled');
colorbar;
hold on;

% 标记真实导数为零的点
scatter(y(zero_dT_idx), x(zero_dT_idx), 150, 'k', 'd', 'LineWidth', 2);

% 标记K值识别出的零线点
scatter(absolute_nullcline(:, 2), absolute_nullcline(:, 1), ...
        100, 'm', 's', 'LineWidth', 1.5);

title('真实效率变化率与零线对比');
xlabel('入口NOx浓度 (ppm)');
ylabel('温度 (°C)');
legend('效率变化率', '真实零导点', 'K值识别零线');
grid on;

% 子图3: 三维K值分布
subplot(2, 2, 3);
surf(delta_y_range, delta_x_range, K_matrix);
shading interp;
colormap(jet);
xlabel('ΔNOx (ppm)');
ylabel('ΔT (°C)');
zlabel('K值');
title('三维K值分布');
view(-30, 30);

% 子图4: K值分布与效率关系
subplot(2, 2, 4);
scatter3(y, x, z, 100, K_matrix_interp(x, y, delta_x_range, delta_y_range, K_matrix), 'filled');
colorbar;
xlabel('入口NOx浓度 (ppm)');
ylabel('温度 (°C)');
zlabel('脱硝效率 (%)');
title('K值分布与脱硝效率关系');
grid on;
view(-30, 30);

%% 识别最佳工况点 (K值最小点)
[~, min_idx] = min(K_matrix(:));
[min_i, min_j] = ind2sub(size(K_matrix), min_idx);

optimal_dT = delta_x_range(min_i);
optimal_dNOx = delta_y_range(min_j);

% 计算最佳工况点
optimal_T = mean(x) + optimal_dT;
optimal_NOx = mean(y) + optimal_dNOx;

fprintf('最佳工况点识别结果:\n');
fprintf('温度: %.2f°C\n', optimal_T);
fprintf('入口NOx浓度: %.2f ppm\n', optimal_NOx);
fprintf('对应K值: %d\n', K_matrix(min_i, min_j));

%% 与真实导数为零点比较
if ~isempty(zero_dT_idx)
    real_zero_points = [y(zero_dT_idx), x(zero_dT_idx)];
    k_zero_points = absolute_nullcline;
    
    % 计算距离矩阵
    distances = pdist2(real_zero_points, k_zero_points);
    min_distances = min(distances, [], 2);
    
    fprintf('\n零线识别精度评估:\n');
    fprintf('平均最小距离: %.2f\n', mean(min_distances));
    fprintf('最大最小距离: %.2f\n', max(min_distances));
    fprintf('匹配点比例 (距离<10): %.2f%%\n', mean(min_distances < 10) * 100);
else
    fprintf('\n警告: 未找到真实导数为零的点\n');
end

%% SCR专用函数实现

% SINDy模型识别函数
function [model, K, R2] = sindy_scr_identification(x, y, z)
    % SCR专用基函数库
    Theta = [ones(size(x)), x, y, x.*y, x.^2, y.^2, ...
             exp(-80000./(8.314*(y+273)))]; % Arrhenius项
    
    % 目标变量
    target = z;
    
    % 稀疏回归 (使用Lasso)
    [B, FitInfo] = lasso(Theta, target, 'Alpha', 0.7, 'CV', 5);
    
    % 选择最佳lambda
    idx = FitInfo.IndexMinMSE;
    Xi = B(:, idx);
    
    % 添加截距项
    Xi = [FitInfo.Intercept(idx); Xi];
    
    % 计算非零项数量 (排除截距)
    K = nnz(abs(Xi(2:end)) > 0.01);
    
    % 计算R²
    pred = [ones(size(Theta,1),1), Theta] * Xi;
    R2 = 1 - sum((target - pred).^2) / sum((target - mean(target)).^2);
    
    % 构建模型对象
    model.Xi = Xi;
    model.predict = @(x,y) [1, x, y, x*y, x^2, y^2, ...
                            exp(-80000/(8.314*(y+273)))] * Xi;
end

% K值插值函数
function K_values = K_matrix_interp(x, y, delta_x_range, delta_y_range, K_matrix)
    % 将点映射到变换空间
    x_trans = x - mean(x);
    y_trans = y - mean(y);
    
    % 插值获取K值
    K_values = interp2(delta_y_range, delta_x_range, K_matrix, ...
                      y_trans, x_trans, 'linear', mean(K_matrix(:)));
end