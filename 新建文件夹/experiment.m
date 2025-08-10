% 生成FHN模型数据
%a = 0; b = 0.5; c = 1; d = 1; eps = 1;
%dt = 0.01;
%tspan = 0:dt:40;
%h = 1;
%u0 = [0.1; 0.1];

%f = @(t,u) [-u(1)^3 + c*u(1)^2 + d*u(1) - u(2); 
 %           eps*(u(1) - b*u(2) + a)];
%[~, u_orig] = ode45(f, tspan, u0);
%X = [u_orig(:,1), u_orig(:,2)];
%du = gradient(u_orig(:,1), dt);
%dv = gradient(u_orig(:,2), dt);
% 运行模型识别
%results = sindy_model_identification_generalnullcline(X, dt, 'PlotResults', true);
%options = struct('deriv_threshold', 0.03, 'u_poly_order', 3);
%[u_expr, v_expr] = identify_nullclines_general(u_orig(:,1), u_orig(:,2) , du, dv, options);

% 参数设置
%a = 1.0; b = 3.0; 
%Du = 0.1; Dv = 0.1; % 扩散系数（若无需空间扩散则设为0）

% 动力学方程
% g=@(t, y) [
%    y(2); % 离散Laplacian简化
%    h*(1-y(1)^2)*y(2)-y(1)
%];

% 数值模拟
%initial_condition=[[0;1];[20;1];[0.5;0];[10;0.5];[12;0];[20;0];[30;0];[35;0];[26;0];[32;0.5];[36;0.2];[0;6];[0;8]];
%[t, Y] = ode15s(g, tspan, [0; 1]); % 初始浓度[0.5, 0.5]
%du = gradient(Y(:,1), dt);
%dv = gradient(Y(:,2), dt);
%[t1, Y1] = ode15s(g, tspan, [6; 1]);
%[t2, Y2] = ode15s(g, tspan, [0.5; 0]);
%[t3, Y3] = ode15s(g, tspan, [5; 0.5]);
%[t4, Y4] = ode15s(g, tspan, [6; 0]);
%[t5, Y5] = ode15s(g, tspan, [4; 0]);
%[t6, Y6] = ode15s(g, tspan, [3; 0]);
%[t7, Y7] = ode15s(g, tspan, [3.5; 0]);
%[t8, Y8] = ode15s(g, tspan, [5.5; 0]);
%[t9, Y9] = ode15s(g, tspan, [7; 0.5]);
%[t10, Y10] = ode15s(g, tspan, [6.2; 0.2]);
%[t11, Y11] = ode15s(g, tspan, [0; 6]);
%[t12, Y12] = ode15s(g, tspan, [0; 8]);
%[t13, Y13] = ode15s(g, tspan, [1; 1]);
%[t14, Y14] = ode15s(g, tspan, [0.1; 1]);
%[t15, Y15] = ode15s(g, tspan, [0; 4]);
%[stacked_states, stacked_derivatives] = prepareSINDyData({Y,Y1, Y2,Y3,Y4,Y5}, {t, t1, t2, t3,t4,t5});
% 运行模型识别

R = 8.314;          % 气体常数 (J/mol·K)
Ea_range = [75000, 95000];  % 活化能搜索范围 (J/mol)
T_ref = 500;        % 参考温度 (K)

% 加载数据
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

% 温度单位转换 (°C → K)
T_K = T_scr + 273.15;



%results = sindy_model_identification_generalnullcline(X, dt, X_der, 'PlotResults', true);
%通过观察，nullcline函数可能含有三角函数项或某些有理函数如sinx,1/x,1/x^2,1/1-x,1/1-x^2。为了拟合出正确形式的nullcline函数，我们手动在Theta中添加这些基函数。

%options = struct('deriv_threshold', 0.03, 'u_poly_order', 3);
%options.u_extra_basis = {
%    @(u) sin(u), @(u) 1./u, @(u) 1./(u-1), @(u) 1./(u.^2), @(u) 1./(u.^2-1)
%};

%[u_expr, v_expr] = identify_nullclines_general(X(:,1), X(:,2) , X_der(:,1), X_der(:,2), options);

%param = [10; 28; 8/3]; % Lorenz system parameters (chaotic)
%n = 3; % number of states
%x0 = [-8; 8; 27];  % Initial condition
%dt = 0.001; % time step
%tspan = dt:dt:20; 
%options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
%[t,y]=ode45(@(t,y) lorenz(t,y,param),tspan,x0,options);
%[t, y] = ode45(@lorenz_system, [0 100], [1; 1; 1]);

%u = y(:,1);  % x分量
%v = y(:,2);  % y分量
%du = gradient(u, t(2)-t(1));
%dv = gradient(v, t(2)-t(1));

% 识别nullcline
%options = struct('deriv_threshold', 0.1, 'v_poly_order', 2);
%[u_expr, v_expr] = identify_nullclines_general(u, v, du, dv, optioduns);
% 示例1: Lorenz系统
%sigma = 10; beta = 8/3; rho = 28;
%dt = 0.01; T = 100; tspan = 0:dt:T;

%lorenz = @(t, x) [sigma*(x(2)-x(1)); 
 %                  x(1)*(rho-x(3))-x(2); 
  %                 x(1)*x(2)-beta*x(3)];

%[~, X] = ode45(lorenz, tspan, [1; 1; 1]);
%dX = zeros(size(X));
%for i = 1:3
 %   dX(:, i) = gradient(X(:, i), dt);
%end

%options = struct('dim1', 1, 'dim2', 2, 'slice_dim', 3, ...
 %                'deriv_threshold', 0.1, 'poly_order', 2, ...
  %               'num_slice_candidates', 40, 'min_points', 10, ...
  %               'max_segments', 2);

%[best_u_expr, best_v_expr, best_slice] = identify_nullclines_general_highdimensions(X, dX, options);
%% 零线重建方法在SCR系统中的应用
%clear all; close all; clc;

% 使用提供的SCR数据生成代码
% ... [这里插入您提供的SCR数据生成代码] ...
% 假设已生成以下变量: NOx_in, urea_flow, T, efficiency, dEff_dT

%% 1. 准备数据
% 将温度转换为开尔文
%T_K = T;

% 定义状态变量 (二维系统)
x = T_K - 273.15;  % 温度 (°C)
y = NOx_in;        % 入口NOx浓度 (ppm)

% 计算状态导数
dt = 1; % 假设时间间隔为1
dx = gradient(x, dt);
dy = gradient(y, dt);

%% 2. 定义变换向量范围
delta_x_range = linspace(-10, 10, 15); % 温度扰动范围 ±10°C
delta_y_range = linspace(-200, 200, 15); % NOx扰动范围 ±200ppm

% 初始化K值矩阵
K_matrix = zeros(length(delta_x_range), length(delta_y_range));
R2_matrix = zeros(size(K_matrix));
delta_matrix = zeros(size(K_matrix));

%% 3. 计算K值分布 (PDF方法核心)
for i = 1:length(delta_x_range)
    for j = 1:length(delta_y_range)
        % 施加变换向量 (PDF方法)
        x_shifted = x + delta_x_range(i);
        y_shifted = y + delta_y_range(j);
        
        % 使用SINDy识别模型 (SCR专用)
        [~, K, R2] = sindy_scr_identification(x_shifted, y_shifted, eff);
        
        % 计算δ值 (吸引子距离)
        delta_val = compute_attractor_distance([x; eff], [x_shifted; eff]);
        
        % 存储结果
        K_matrix(i, j) = K;
        R2_matrix(i, j) = R2;
        delta_matrix(i, j) = delta_val;
    end
end

%% 4. 识别零线位置 (K值下降点)
% 计算K值阈值 (PDF方法中的关键观察)
K_threshold = prctile(K_matrix(:), 30); % 取K值最低的30%作为阈值

% 找出低K值区域 (潜在零线位置)
[nullcline_x, nullcline_y] = find(K_matrix < K_threshold);

% 转换为实际坐标
nullcline_points = [delta_x_range(nullcline_x)', delta_y_range(nullcline_y)'];

%% 5. 可视化结果
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
zero_dT_idx = find(abs(dEff_dT) < 0.05);
scatter(y(zero_dT_idx), x(zero_dT_idx), 150, 'k', 'd', 'LineWidth', 2);

% 标记K值识别出的零线点
scatter(nullcline_points(:, 2) + mean(y), nullcline_points(:, 1) + mean(x), ...
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
scatter3(y, x, efficiency, 100, K_matrix_grid, 'filled');
colorbar;
xlabel('入口NOx浓度 (ppm)');
ylabel('温度 (°C)');
zlabel('脱硝效率 (%)');
title('K值分布与脱硝效率关系');
grid on;
view(-30, 30);

%% 6. 识别最佳工况点 (K值最小点)
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

%% 7. 与真实导数为零点比较
real_zero_points = [y(zero_dT_idx), x(zero_dT_idx)];
k_zero_points = [nullcline_points(:, 2) + mean(y), nullcline_points(:, 1) + mean(x)];

% 计算距离矩阵
distances = pdist2(real_zero_points, k_zero_points);
min_distances = min(distances, [], 2);

fprintf('\n零线识别精度评估:\n');
fprintf('平均最小距离: %.2f\n', mean(min_distances));
fprintf('最大最小距离: %.2f\n', max(min_distances));
fprintf('匹配点比例 (距离<10): %.2f%%\n', mean(min_distances < 10) * 100);

%% SCR专用函数实现

% SINDy模型识别函数
function [model, K, R2] = sindy_scr_identification(x, y, efficiency)
    % SCR专用基函数库
    Theta = [ones(size(x)), x, y, ...
             x.*y, x.^2, y.^2, ...
             exp(-80000./(8.314*(y+273))), ... % Arrhenius项
             log(x+1), log(y+273)];             % 对数项
    
    % 目标变量
    target = efficiency;
    
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
                            exp(-80000/(8.314*(y+273))), ...
                            log(x+1), log(y+273)] * Xi;
end

% 吸引子距离计算 (δ值)
function delta = compute_attractor_distance(orig_data, perturbed_data)
    % 计算主成分分析（PCA）投影
    [coeff, score_orig] = pca(orig_data');
    
    % 投影扰动数据
    score_pert = coeff(:,1:2)' * perturbed_data;
    
    % 计算Hausdorff距离
    D = pdist2(score_orig(:,1:2), score_pert');
    min_dist_orig = min(D, [], 2);
    min_dist_pert = min(D, [], 1);
    
    delta = max(prctile(min_dist_orig, 95), prctile(min_dist_pert, 95));
end

% 网格化K值矩阵 (用于3D可视化)
function K_grid = compute_K_grid(x, y, efficiency, delta_x_range, delta_y_range)
    [X_grid, Y_grid] = meshgrid(linspace(min(x), max(x), 20), ...
                               linspace(min(y), max(y), 20));
    K_grid = zeros(size(X_grid));
    
    for i = 1:size(X_grid, 1)
        for j = 1:size(X_grid, 2)
            % 找到最近的delta点
            [~, dx_idx] = min(abs(delta_x_range - (X_grid(i,j) - mean(x))));
            [~, dy_idx] = min(abs(delta_y_range - (Y_grid(i,j) - mean(y))));
            
            K_grid(i,j) = K_matrix(dx_idx, dy_idx);
        end
    end
end