%% SCR系统时间序列数据生成
clear all; close all; clc;

% 加载原始数据（来自之前的实验）
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

%% 时间序列转换（添加时间维度）
dt = 0.5; % 采样时间间隔（分钟）
t = (0:dt:(length(NOx_in)-1)*dt)'; % 时间向量

% 添加动态特性（时间延迟和瞬态响应）
% 1. 入口NOx动态（考虑发动机负载变化）
NOx_in_ts = zeros(size(NOx_in));
for i = 1:length(NOx_in)
    if i == 1
        NOx_in_ts(i) = NOx_in(i);
    else
        % 带一阶惯性环节的动态
        tau_NOx = 2.0; % NOx变化时间常数（分钟）
        NOx_in_ts(i) = NOx_in_ts(i-1) + (NOx_in(i) - NOx_in_ts(i-1)) * (dt/tau_NOx);
    end
end

% 2. 尿素喷射动态（考虑尿素水解延迟）
urea_ts = zeros(size(urea));
for i = 1:length(urea)
    if i == 1
        urea_ts(i) = urea(i);
    else
        % 带一阶惯性环节的动态
        tau_urea = 1.5; % 尿素响应时间常数（分钟）
        urea_ts(i) = urea_ts(i-1) + (urea(i) - urea_ts(i-1)) * (dt/tau_urea);
    end
end

% 3. SCR温度动态（考虑热惯性）
T_scr_ts = zeros(size(T_scr));
for i = 1:length(T_scr)
    if i == 1
        T_scr_ts(i) = T_scr(i);
    else
        % 带一阶惯性环节的动态
        tau_temp = 3.0; % 温度响应时间常数（分钟）
        T_scr_ts(i) = T_scr_ts(i-1) + (T_scr(i) - T_scr_ts(i-1)) * (dt/tau_temp);
    end
end

% 4. 出口NOx动态（考虑化学反应动力学）
NOx_out_ts = zeros(size(NOx_out));
for i = 1:length(NOx_out)
    if i == 1
        NOx_out_ts(i) = NOx_out(i);
    else
        % 基于温度的反应速率动态
        k_react = 0.05 * exp(-8000/(8.314*(T_scr_ts(i)+273.15))); % 反应速率常数
        tau_react = 1/(k_react + 0.01); % 反应时间常数
        
        % 一阶动态模型
        NOx_out_ts(i) = NOx_out_ts(i-1) + (NOx_out(i) - NOx_out_ts(i-1)) * (dt/tau_react);
    end
end

% 5. 效率动态（基于出口NOx计算）
eff_ts = 100 * (1 - NOx_out_ts ./ NOx_in_ts);

%% 构建SCR动力系统状态空间模型
% 状态变量: [入口NOx浓度; 催化剂温度; 表面氨覆盖率; 出口NOx浓度]
% 输入变量: [发动机负载; 环境温度; 尿素喷射量]

% 系统参数（基于物理模型）
R = 8.314; % 气体常数
Ea = 80000; % 活化能 (J/mol)
k_ads = 0.02; % 氨吸附速率常数
k_des = 0.005; % 氨解吸速率常数
k_react = 0.001; % 主反应速率常数

% 状态空间方程离散化
A = zeros(4,4); % 状态转移矩阵
B = zeros(4,3); % 输入矩阵
C = [0 0 0 1];  % 输出矩阵（出口NOx）
D = 0;          % 直接传输矩阵

% 时间序列状态估计
X = zeros(4, length(t)); % 状态矩阵
Y = zeros(1, length(t)); % 输出矩阵

% 初始状态
X(:,1) = [NOx_in_ts(1); T_scr_ts(1); 0.5; NOx_out_ts(1)];

% 输入向量（发动机负载、环境温度、尿素喷射量）
% 简化假设：发动机负载与入口NOx成正比，环境温度=20°C
U = [NOx_in_ts/1000, 20*ones(size(t)), urea_ts]';

for k = 1:length(t)-1
    % 当前状态
    NOx_in_k = X(1,k);
    T_k = X(2,k);
    theta_NH3_k = X(3,k);
    NOx_out_k = X(4,k);
    
    % 输入
    engine_load = U(1,k);
    T_amb = U(2,k);
    urea_flow = U(3,k);
    
    % 系统动力学方程 (连续时间)
    dNOx_in_dt = 0.1*(engine_load*1000 - NOx_in_k); % 入口NOx动态
    dT_dt = 0.05*(T_scr_ts(k) - T_k) + 0.01*(T_amb - T_k); % 温度动态
    dtheta_NH3_dt = k_ads*urea_flow*(1 - theta_NH3_k) - k_des*theta_NH3_k; % 氨覆盖率动态
    dNOx_out_dt = k_react*exp(-Ea/(R*(T_k+273.15)))*theta_NH3_k*NOx_in_k - 0.2*NOx_out_k; % 出口NOx动态
    
    % 离散化（前向欧拉法）
    X(1,k+1) = X(1,k) + dt * dNOx_in_dt;
    X(2,k+1) = X(2,k) + dt * dT_dt;
    X(3,k+1) = X(3,k) + dt * dtheta_NH3_dt;
    X(4,k+1) = X(4,k) + dt * dNOx_out_dt;
    
    % 输出
    Y(k) = C * X(:,k);
end

% 最终输出
Y(end) = C * X(:,end);

%% 时间序列可视化
figure('Position', [100, 100, 1400, 1000]);

% 入口NOx浓度
subplot(3,2,1);
plot(t, NOx_in_ts, 'b', 'LineWidth', 2);
hold on;
plot(t, X(1,:), 'r--', 'LineWidth', 1.5);
title('入口NOx浓度动态');
xlabel('时间 (分钟)');
ylabel('NOx (ppm)');
legend('测量值', '模型预测');
grid on;

% SCR温度
subplot(3,2,2);
plot(t, T_scr_ts, 'b', 'LineWidth', 2);
hold on;
plot(t, X(2,:), 'r--', 'LineWidth', 1.5);
title('SCR温度动态');
xlabel('时间 (分钟)');
ylabel('温度 (°C)');
legend('测量值', '模型预测');
grid on;

% 尿素喷射量
subplot(3,2,3);
plot(t, urea_ts, 'g', 'LineWidth', 2);
title('尿素喷射量');
xlabel('时间 (分钟)');
ylabel('尿素 (kg/h)');
grid on;

% 氨覆盖率（模型状态）
subplot(3,2,4);
plot(t, X(3,:), 'm', 'LineWidth', 2);
title('催化剂表面氨覆盖率');
xlabel('时间 (分钟)');
ylabel('覆盖率 (\theta_{NH_3})');
ylim([0 1]);
grid on;

% 出口NOx浓度
subplot(3,2,5);
plot(t, NOx_out_ts, 'b', 'LineWidth', 2);
hold on;
plot(t, Y, 'r--', 'LineWidth', 1.5);
title('出口NOx浓度动态');
xlabel('时间 (分钟)');
ylabel('NOx (ppm)');
legend('测量值', '模型预测');
grid on;

% 脱硝效率
subplot(3,2,6);
eff_model = 100*(1 - Y'./NOx_in_ts);
plot(t, eff_ts, 'b', 'LineWidth', 2);
hold on;
plot(t, eff_model, 'r--', 'LineWidth', 1.5);
title('脱硝效率动态');
xlabel('时间 (分钟)');
ylabel('效率 (%)');
legend('测量值', '模型预测');
grid on;

%% 系统辨识 - 使用SINDy方法识别动力学方程
% 状态变量: x = [NOx_in; T_scr; theta_NH3; NOx_out]
% 输入变量: u = [engine_load; T_amb; urea_flow]

% 构建特征库
n = length(t)-1;
X_current = X(:,1:end-1)'; % 当前状态
X_next = X(:,2:end)';     % 下一时刻状态
dX_dt = (X_next - X_current)/dt; % 状态导数

% 输入矩阵
U_current = U(:,1:end-1)';

% 特征库函数
feature_library = @(x,u) [
    ones(size(x,1),1), ...                  % 常数项
    x, ...                                  % 状态变量
    u, ...                                  % 输入变量
    x(:,1).*x(:,4), ...                     % NOx_in * NOx_out
    x(:,2).*x(:,3), ...                     % T_scr * theta_NH3
    x(:,1).*u(:,3), ...                     % NOx_in * urea
    exp(-Ea./(R*(x(:,2)+273.15))), ...       % Arrhenius项
    x(:,3).*x(:,1).*exp(-Ea./(R*(x(:,2)+273.15))), ... % 反应速率项
    x(:,1).^2, x(:,2).^2, x(:,3).^2, x(:,4).^2, ... % 二次项
    x(:,1).*x(:,2), x(:,1).*x(:,3), x(:,1).*x(:,4), ... % 交叉项
    x(:,2).*x(:,3), x(:,2).*x(:,4), x(:,3).*x(:,4), ...
    sin(2*pi*x(:,2)/100), cos(2*pi*x(:,2)/100) ... % 周期性温度影响
];

% 构建大特征矩阵
Theta = feature_library(X_current, U_current);

% 使用稀疏回归识别动力学方程
lambda = 0.1; % 稀疏参数
Xi = zeros(size(Theta,2),4); % 系数矩阵

for i = 1:4
    % 使用Lasso回归
    [B, FitInfo] = lasso(Theta, dX_dt(:,i), 'Lambda', lambda);
    Xi(:,i) = B;
    
    % 计算拟合优度
    y_pred = Theta * B;
    R2(i) = 1 - sum((dX_dt(:,i) - y_pred).^2)/sum((dX_dt(:,i) - mean(dX_dt(:,i))).^2);
end

% 显示识别结果
state_names = {'dNOx_in/dt', 'dT_scr/dt', 'dtheta_NH3/dt', 'dNOx_out/dt'};
feature_names = {
    '1', 'NOx_in', 'T_scr', 'theta_NH3', 'NOx_out', ...
    'Engine_load', 'T_amb', 'Urea_flow', ...
    'NOx_in*NOx_out', 'T_scr*theta_NH3', 'NOx_in*Urea', ...
    'exp(-Ea/RT)', 'theta_NH3*NOx_in*exp(-Ea/RT)', ...
    'NOx_in^2', 'T_scr^2', 'theta_NH3^2', 'NOx_out^2', ...
    'NOx_in*T_scr', 'NOx_in*theta_NH3', 'NOx_in*NOx_out', ...
    'T_scr*theta_NH3', 'T_scr*NOx_out', 'theta_NH3*NOx_out', ...
    'sin(2πT_scr/100)', 'cos(2πT_scr/100)'
};

fprintf('识别出的SCR系统动力学方程:\n');
for i = 1:4
    fprintf('\n%s = ', state_names{i});
    nonzero_idx = find(abs(Xi(:,i)) > 0.01);
    for j = 1:length(nonzero_idx)
        idx = nonzero_idx(j);
        if Xi(idx,i) > 0 && j > 1
            fprintf('+');
        end
        fprintf('%.4f*%s ', Xi(idx,i), feature_names{idx});
    end
    fprintf('\t(R²=%.4f)', R2(i));
end

%% 相空间分析
figure('Position', [100, 100, 1200, 900]);

% 3D相图: NOx_in vs T_scr vs 效率
subplot(2,2,1);
scatter3(X(1,:), X(2,:), eff_model, 40, t, 'filled');
colorbar;
xlabel('入口NOx (ppm)');
ylabel('SCR温度 (°C)');
zlabel('脱硝效率 (%)');
title('系统相空间轨迹');
grid on;

% 氨覆盖率 vs 温度 vs 反应速率
subplot(2,2,2);
reaction_rate = k_react * exp(-Ea./(R*(X(2,:)+273.15))) .* X(3,:) .* X(1,:);
scatter3(X(3,:), X(2,:), reaction_rate, 40, t, 'filled');
colorbar;
xlabel('氨覆盖率 (\theta_{NH_3})');
ylabel('SCR温度 (°C)');
zlabel('反应速率');
title('催化反应动力学相图');
grid on;

% 效率变化率 vs 温度
subplot(2,2,3);
dEff_dt = gradient(eff_model, t);
scatter(X(2,:), dEff_dt, 40, X(1,:), 'filled');
colorbar;
xlabel('SCR温度 (°C)');
ylabel('dEfficiency/dt (%/min)');
title('效率变化率相图');
grid on;
hold on;
plot([min(X(2,:)), max(X(2,:))], [0,0], 'r--');

% 状态变量自相关分析
subplot(2,2,4);
autocorr_NOx = autocorr(X(1,:), 'NumLags', 20);
autocorr_T = autocorr(X(2,:), 'NumLags', 20);
autocorr_theta = autocorr(X(3,:), 'NumLags', 20);
plot(0:20, autocorr_NOx, 'b-o', 'LineWidth', 1.5);
hold on;
plot(0:20, autocorr_T, 'r-s', 'LineWidth', 1.5);
plot(0:20, autocorr_theta, 'g-^', 'LineWidth', 1.5);
title('状态变量自相关函数');
xlabel('时间滞后 (采样点)');
ylabel('自相关系数');
legend('入口NOx', 'SCR温度', '氨覆盖率');
grid on;

%% 动态模型验证 - 阶跃响应测试
% 设置阶跃输入条件
t_test = (0:dt:30)'; % 30分钟测试
n_test = length(t_test);

% 初始条件
x0 = [1000; 200; 0.5; 800]; % [NOx_in; T_scr; theta_NH3; NOx_out]

% 输入信号
u_test = zeros(3, n_test);
u_test(1,:) = 1.0; % 发动机负载 = 1.0 (对应NOx_in=1000ppm)
u_test(2,:) = 20;  % 环境温度 = 20°C
u_test(3,:) = 0.7; % 尿素喷射量 = 0.7 kg/h (前15分钟)
u_test(3, floor(n_test/2):end) = 1.4; % 15分钟后阶跃到1.4 kg/h

% 模拟系统响应
X_test = zeros(4, n_test);
X_test(:,1) = x0;

for k = 1:n_test-1
    % 使用识别出的动力学方程
    Theta_k = feature_library(X_test(:,k)', u_test(:,k)');
    dX_dt = Theta_k * Xi;
    
    % 欧拉积分
    X_test(:,k+1) = X_test(:,k) + dt * dX_dt';
end

% 计算效率
eff_test = 100 * (1 - X_test(4,:)' ./ X_test(1,:)');

% 绘制阶跃响应
figure('Position', [100, 100, 1200, 800]);

% 尿素喷射量阶跃
subplot(3,1,1);
plot(t_test, u_test(3,:), 'g', 'LineWidth', 2);
title('尿素喷射量阶跃输入');
xlabel('时间 (分钟)');
ylabel('尿素 (kg/h)');
grid on;

% 出口NOx响应
subplot(3,1,2);
plot(t_test, X_test(4,:), 'b', 'LineWidth', 2);
title('出口NOx浓度响应');
xlabel('时间 (分钟)');
ylabel('NOx (ppm)');
grid on;

% 脱硝效率响应
subplot(3,1,3);
plot(t_test, eff_test, 'r', 'LineWidth', 2);
title('脱硝效率动态响应');
xlabel('时间 (分钟)');
ylabel('效率 (%)');
grid on;

%% 保存时间序列数据
time_series_data = [t, NOx_in_ts, T_scr_ts, urea_ts, NOx_out_ts, eff_ts];
save('scr_time_series_data.mat', 'time_series_data', 'X', 'U', 't');

fprintf('\n时间序列数据已保存为 scr_time_series_data.mat\n');