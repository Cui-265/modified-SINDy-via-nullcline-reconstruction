function [X_stack, dX_stack] = prepareSINDyData(trajectories, time_vectors, options)
% PREPARESINDYDATA 整合多条轨迹数据用于SINDy算法
%   输入：
%       trajectories: 元胞数组，每条轨迹的状态数据
%       time_vectors: 元胞数组，每条轨迹的时间向量
%       options: 可选参数结构体
%           .derivative_method: 求导方法 ('central', 'spline', 'savitzky-golay')
%           .smooth_derivatives: 是否平滑导数 (true/false)
%   输出：
%       X_stack:  堆叠后的状态数据
%       dX_stack: 堆叠后的导数数据

% 设置默认选项
if nargin < 3
    options = struct();
end
if ~isfield(options, 'derivative_method')
    options.derivative_method = 'central';
end
if ~isfield(options, 'smooth_derivatives')
    options.smooth_derivatives = true;
end

% 验证输入
validateInput(trajectories, time_vectors);

% 初始化堆叠矩阵
X_stack = [];
dX_stack = [];

% 处理每条轨迹
for i = 1:numel(trajectories)
    X = trajectories{i};
    t = time_vectors{i};
    
    % 计算数值导数
    dX = computeDerivatives(X, t, options);
    
    % 堆叠数据
    X_stack = [X_stack; X];
    dX_stack = [dX_stack; dX];
end

% 验证输出
validateOutput(X_stack, dX_stack);

fprintf('成功整合%d条轨迹数据\n', numel(trajectories));
fprintf('总数据点数: %d, 状态变量数: %d\n', size(X_stack, 1), size(X_stack, 2));
end

%% 辅助函数
function validateInput(trajectories, time_vectors)
if ~iscell(trajectories) || ~iscell(time_vectors)
    error('输入必须是元胞数组');
end
if numel(trajectories) ~= numel(time_vectors)
    error('轨迹数据和时间向量数量不一致');
end
end

function validateOutput(X_stack, dX_stack)
if size(X_stack, 1) ~= size(dX_stack, 1) || size(X_stack, 2) ~= size(dX_stack, 2)
    error('状态数据和导数数据维度不匹配');
end
end

function dX = computeDerivatives(X, t, options)
% 高级导数计算
n = size(X, 1);
if n == 1
    dX = zeros(1, size(X, 2));
    return;
end

switch options.derivative_method
    case 'spline'
        % 样条插值求导
        dX = zeros(size(X));
        for j = 1:size(X, 2)
            pp = spline(t, X(:, j));
            dpp = fnder(pp);  % 导数样条
            dX(:, j) = ppval(dpp, t);
        end
        
    case 'savitzky-golay'
        % Savitzky-Golay 滤波求导
        order = min(3, n-1);  % 多项式阶数
        framelen = min(min(11, n), floor(n/2)*2+1);  % 奇数窗口
        
        dX = zeros(size(X));
        for j = 1:size(X, 2)
            [~, g] = sgolay(order, framelen);
            dt = mean(diff(t));
            halfwin = (framelen-1)/2;
            
            for i = 1:n
                % 处理边界
                start_idx = max(1, i-halfwin);
                end_idx = min(n, i+halfwin);
                segment = X(start_idx:end_idx, j);
                
                % 调整窗口中心
                center = i - start_idx + 1;
                
                % 计算导数
                dX(i, j) = g(center, 2)' * segment / dt;
            end
        end
        
    otherwise % 'central' 或其他
        % 向量化中心差分法
        dX = centralDifference(X, t);
end

% 可选平滑
if options.smooth_derivatives && n > 5
    dX = smoothDerivatives(dX);
end
end

function dX = centralDifference(X, t)
% 向量化中心差分实现
n = size(X, 1);
dX = zeros(size(X));

% 边界处理
if n > 1
    dX(1, :) = (X(2, :) - X(1, :)) / (t(2) - t(1));
end
if n > 2
    dX(end, :) = (X(end, :) - X(end-1, :)) / (t(end) - t(end-1));
end

% 内部点（向量化）
if n > 2
    dt_forward = diff(t(2:end));   % t(i+1) - t(i)
    dt_backward = diff(t(1:end-1)); % t(i) - t(i-1)
    total_dt = dt_backward + dt_forward;
    
    % 避免除以零
    total_dt(total_dt == 0) = eps;
    
    dX(2:end-1, :) = (X(3:end, :) - X(1:end-2, :)) ./ total_dt;
end
end

function dX_smooth = smoothDerivatives(dX)
% 导数平滑处理
[n, m] = size(dX);
dX_smooth = zeros(size(dX));
window_size = min(5, floor(n/2));  % 自适应窗口大小

for j = 1:m
    % 使用移动中值滤波去除离群点
    dX_clean = movmedian(dX(:, j), window_size, 'omitnan');
    
    % 应用高斯平滑
    dX_smooth(:, j) = smoothdata(dX_clean, 'gaussian', window_size);
end
end