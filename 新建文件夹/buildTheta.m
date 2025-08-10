function Theta = buildTheta(input, poly_order, extra_basis)
% 构建包含多项式项、额外基函数及其乘积的特征矩阵
% 输入:
%   input: 输入向量 (N×1)
%   poly_order: 多项式阶数
%   extra_basis: 额外基函数元胞数组 {f1, f2, ...}
% 输出:
%   Theta: 特征矩阵 (N×M)

    % 步骤1: 创建基础特征集
    n = length(input);
    base_features = {};
    
    % 添加原始输入
    base_features{end+1} = input;
    feature_names = {'x'};
    
    % 添加额外基函数
    if nargin >= 3 && ~isempty(extra_basis)
        for i = 1:length(extra_basis)
            f = extra_basis{i};
            base_features{end+1} = f(input);
            feature_names{end+1} = ['f', num2str(i), '(x)'];
        end
    end
    
    num_base_features = length(base_features);
    
    % 步骤2: 初始化Theta
    Theta = ones(n, 1);  % 常数项
    term_descriptions = {'1'};  % 用于记录每项的描述
    
    % 步骤3: 添加所有可能的乘积组合
    for k = 1:poly_order
        % 获取所有可能的k元组合（允许重复）
        combinations = nchoosek_with_replacement(1:num_base_features, k);
        
        for j = 1:size(combinations, 1)
            term = ones(n, 1);
            term_name_parts = {};
            
            % 计算该项的乘积
            for idx = 1:k
                feature_idx = combinations(j, idx);
                term = term .* base_features{feature_idx};
                
                % 构建该项的描述
                term_name_parts{end+1} = feature_names{feature_idx};
            end
            
            % 添加到Theta矩阵
            Theta = [Theta, term];
            
            % 构建完整的项描述
            term_desc = strjoin(term_name_parts, '*');
            term_descriptions{end+1} = term_desc;
        end
    end
    
    % 显示特征矩阵信息
    fprintf('构建的特征矩阵包含 %d 个特征项:\n', size(Theta, 2));
    for i = 1:min(20, length(term_descriptions))  % 最多显示20项
        fprintf('  Term %2d: %s\n', i, term_descriptions{i});
    end
    if length(term_descriptions) > 20
        fprintf('  ... 共 %d 项\n', length(term_descriptions));
    end
end

function combinations = nchoosek_with_replacement(set, k)
% 生成允许重复的组合
% 例如: nchoosek_with_replacement([1,2,3], 2) 返回:
%   [1,1; 1,2; 1,3; 2,2; 2,3; 3,3]

    combinations = [];
    indices = ones(1, k);  % 初始组合
    
    while true
        % 添加当前组合（非递减顺序）
        combinations = [combinations; indices];
        
        % 找到下一个组合
        i = k;
        while i >= 1 && indices(i) == max(set)
            i = i - 1;
        end
        
        if i < 1
            break;  % 所有组合已生成
        end
        
        indices(i) = indices(i) + 1;
        indices(i+1:k) = indices(i);
    end
end