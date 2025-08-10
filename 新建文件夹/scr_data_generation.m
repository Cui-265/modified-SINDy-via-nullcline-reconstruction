% SCR系统脱硝性能数据生成的MATLAB代码
clear all; close all; clc;

%% 模型参数设置
% 参考论文中的动力学参数
k1 = 0.001;     % 主反应速率常数
k2 = 0.0005;    % 副反应速率常数
Ea = 80000;     % 活化能 (J/mol)
R = 8.314;      % 气体常数 (J/(mol·K))
T_ref = 500;    % 参考温度 (K)

%% 生成无控制变量的实验数据
% 随机生成入口NOx浓度 (500-3000 ppm)
NOx_in = randi([500, 3000], 20, 1);

% 尿素喷射量与入口NOx浓度成比例 (尿素:NOx摩尔比约为1:1)
urea_flow = NOx_in * 0.0007 + randn(size(NOx_in)) * 0.1;
urea_flow = max(0.5, min(2.2, urea_flow)); % 限制在合理范围

% 随机生成SCR温度 (180-290°C，转换为K)
T = randi([180, 290], 20, 1) + 273.15;

% 初始化结果数组
NOx_out = zeros(20, 1);
efficiency = zeros(20, 1);
dEff_dT = zeros(20, 1);
dEff_dNOx = zeros(20, 1);
dEff_dUrea = zeros(20, 1);

%% 计算脱硝效率 (基于SCR反应动力学模型)
for i = 1:length(NOx_in)
    % 温度对反应速率的影响 (Arrhenius方程)
    k_T = k1 * exp(-Ea/R * (1/T(i) - 1/T_ref));
    
    % 计算平衡常数 (简化模型)
    K = 100 * exp(5000/T(i));
    
    % 计算脱硝效率 (考虑温度、入口NOx和尿素量的综合影响)
    % 参考论文中的效率公式 (如式84-85)
    efficiency(i) = 100 * (1 - exp(-k_T * urea_flow(i) * T(i) / (NOx_in(i) * 0.1)));
    
    % 限制效率在0-100%之间
    efficiency(i) = max(0, min(100, efficiency(i)));
    
    % 计算出口NOx浓度
    NOx_out(i) = NOx_in(i) * (1 - efficiency(i)/100);
    
    % 计算效率对各变量的导数 (使用中心差分法)
    if i > 1 && i < length(NOx_in)
        % 对温度的导数
        dEff_dT(i) = (efficiency(i+1) - efficiency(i-1)) / (T(i+1) - T(i-1));
        
        % 对入口NOx的导数
        dEff_dNOx(i) = (efficiency(i+1) - efficiency(i-1)) / (NOx_in(i+1) - NOx_in(i-1));
        
        % 对尿素量的导数
        dEff_dUrea(i) = (efficiency(i+1) - efficiency(i-1)) / (urea_flow(i+1) - urea_flow(i-1));
    end
end

%% 找出导数为0的点
% 定义导数接近0的阈值
threshold = 0.05;

% 寻找导数为0的点
zero_derivative_points = find(abs(dEff_dT) < threshold | abs(dEff_dNOx) < threshold | abs(dEff_dUrea) < threshold);

%% 显示结果
% 创建表格
data_table = table(NOx_in, urea_flow, T-273.15, NOx_out, efficiency, dEff_dT, 'VariableNames', {'入口NOx(ppm)', '尿素喷射量(kg/h)', 'SCR温度(°C)', '出口NOx(ppm)', '脱硝效率(%)', '效率变化率(%/°C)'});

% 显示完整数据
disp('SCR系统脱硝性能数据:');
disp(data_table);

% 显示导数为0的关键点
disp('导数为0的关键点:');
disp(data_table(zero_derivative_points, :));

%% 可视化数据
figure('Position', [100, 100, 1200, 800]);

% 绘制脱硝效率 vs 温度
subplot(2, 2, 1);
scatter(T-273.15, efficiency, 50, NOx_in, 'filled');
colorbar;
title('脱硝效率 vs 温度');
xlabel('温度 (°C)');
ylabel('脱硝效率 (%)');
grid on;

% 标记导数为0的点
hold on;
scatter((T(zero_derivative_points)-273.15), efficiency(zero_derivative_points), 100, 'r', 'filled', 'MarkerEdgeColor', 'k');
legend('数据点', '导数为0的点');

% 绘制脱硝效率 vs 入口NOx浓度
subplot(2, 2, 2);
scatter(NOx_in, efficiency, 50, T-273.15, 'filled');
colorbar;
title('脱硝效率 vs 入口NOx浓度');
xlabel('入口NOx (ppm)');
ylabel('脱硝效率 (%)');
grid on;

% 绘制脱硝效率 vs 尿素喷射量
subplot(2, 2, 3);
scatter(urea_flow, efficiency, 50, T-273.15, 'filled');
colorbar;
title('脱硝效率 vs 尿素喷射量');
xlabel('尿素喷射量 (kg/h)');
ylabel('脱硝效率 (%)');
grid on;

% 绘制效率变化率 vs 温度
subplot(2, 2, 4);
scatter(T-273.15, dEff_dT, 50, NOx_in, 'filled');
colorbar;
title('效率变化率 vs 温度');
xlabel('温度 (°C)');
ylabel('效率变化率 (%/°C)');
grid on;
hold on;
plot([min(T-273.15), max(T-273.15)], [0, 0], 'r--', 'LineWidth', 2);    