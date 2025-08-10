function Theta_new = add_basis_functions(Theta, X, new_basis)
% 向现有特征矩阵Theta中添加自定义基函数
% 输入:
%   Theta - 现有特征矩阵 (N×M矩阵)
%   X - 原始输入数据 (N×1向量)
%   new_basis - 要添加的新基函数，可以是：
%               1. 单个函数句柄
%               2. 函数句柄的元胞数组
% 输出:
%   Theta_new - 添加新基函数后的特征矩阵

% 验证输入
if size(Theta, 1) ~= size(X, 1)
    error('Theta的行数必须与X的长度匹配');
end

if ~isa(new_basis, 'function_handle') && ~iscell(new_basis)
    error('new_basis必须是函数句柄或函数句柄的元胞数组');
end

% 如果是单个函数句柄，转换为元胞数组
if isa(new_basis, 'function_handle')
    new_basis = {new_basis};
end

% 初始化新特征矩阵
Theta_new = Theta;

% 添加每个新基函数
for i = 1:length(new_basis)
    func = new_basis{i};
    
    % 验证是否为有效函数句柄
    if ~isa(func, 'function_handle')
        error('new_basis中的元素%d不是函数句柄', i);
    end
    
    % 应用基函数并添加到矩阵
    try
        new_col = func(X);
    catch ME
        error('基函数%d执行失败: %s', i, ME.message);
    end
    
    % 验证结果维度
    if ~isequal(size(new_col), size(X))
        error('基函数%d的输出维度与输入不匹配', i);
    end
    
    % 添加到特征矩阵
    Theta_new = [Theta_new, new_col];
end

% 显示添加信息
fprintf('已添加 %d 个新基函数\n', length(new_basis));
fprintf('特征矩阵维度从 %d×%d 扩展到 %d×%d\n', ...
        size(Theta, 1), size(Theta, 2), ...
        size(Theta_new, 1), size(Theta_new, 2));
end