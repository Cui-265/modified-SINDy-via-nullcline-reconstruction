a1 = -0.3;   % y^3 项系数 (dx/dt)
a2 = 0.05;   % y^5 项系数 (dx/dt)
b1 = -0.2;   % x^3 项系数 (dy/dt)
b2 = -0.01;  % x^5 项系数 (dy/dt)
c = -0.4;    % y^3 项系数 (dy/dt)
dt = 0.01;
tspan = 0:dt:40;
u0 = [0.1; 0.1];

% 定义ODE系统
f = @(t, u) [-u(1)^3 + u(1)^2 + u(1)-u(2);
   -0.5*u(2)^3 + 0.5*u(2)^2-(1/3)*u(2) + u(1)];
[t, Y] = ode15s(f, tspan, [0; 1]); % 初始浓度[0.5, 0.5] x'=1-4x+x^2*y, y'=3x-x^2*y x-null:y=(4x-1)/x^2. y-null:x=0 or x=3/y
du = gradient(Y(:,1), dt);
dv = gradient(Y(:,2), dt);
[t1, Y1] = ode15s(f, tspan, [-4; 0.5]);
[t2, Y2] = ode15s(f, tspan, [1; -3]);
[t3, Y3] = ode15s(f, tspan, [3; -3.5]);
%[t4, Y4] = ode15s(f, tspan, [6; 0]);
%[t5, Y5] = ode15s(f, tspan, [4; 0]);
%[t6, Y6] = ode15s(f, tspan, [3; 0]);
%[t7, Y7] = ode15s(f, tspan, [3.5; 0]);
%[t8, Y8] = ode15s(f, tspan, [5.5; 0]);
%[t9, Y9] = ode15s(f, tspan, [7; 0.5]);
%[t10, Y10] = ode15s(f, tspan, [6.2; 0.2]);
%[t11, Y11] = ode15s(f, tspan, [0; 6]);
%[t12, Y12] = ode15s(f, tspan, [0; 8]);
%[t13, Y13] = ode15s(f, tspan, [1; 1]);
%[t14, Y14] = ode15s(f, tspan, [0.1; 1]);
%[t15, Y15] = ode15s(f, tspan, [0; 4]);
[stacked_states, stacked_derivatives] = prepareSINDyData({Y, Y1,Y2,Y3}, {t,t1,t2,t3});
Theta = [ones(length(stacked_states(:,1)),1), ...  % 常数项
                 stacked_states(:,1), ...                 % u
                 stacked_states(:,2), ...                 % v
                 stacked_states(:,1).^2, ...              % u^2
                stacked_states(:,1).^3, ...              % u^3
                 stacked_states(:,2).^2, ...              % v^2
                 stacked_states(:,2).^3];                 % v^3

% 使用ode45解算

results = sindy_model_identification_generalnullcline(stacked_states, dt, stacked_derivatives, 'PlotResults', true);
%results = sindy_model_identification_generalnullcline_artificial(stacked_states, dt, stacked_derivatives, 'PlotResults', true);
%options = struct('deriv_threshold', 0.03, 'u_poly_order', 3);
%[u_expr, v_expr] = identify_nullclines_general(stacked_states(:,1), stacked_states(:,2) , stacked_derivatives(:,1), stacked_derivatives(:,2), options);

%Theta = poolData(stacked_states,2,3);
%results = sparsifyDynamics(Theta,stacked_derivatives,0.05,2);

% 系数矩阵 results (Xi) 的结构:
% 第1列: du/dt 的系数
% 第2列: dv/dt 的系数
% 行顺序: [const, u, v, u^2, u^3, v^2, v^3]

% 打印识别结果
fprintf('识别出的动力系统方程:\n\n');
fprintf('du/dt = ');
fprintf('%5.3f ', results(:,1));
fprintf('\n');

fprintf('dv/dt = ');
fprintf('%5.3f ', results(:,2));
fprintf('\n\n');

% 转换为数学表达式
term_names = {'', 'u', 'v', 'u^2', 'u^3', 'v^2', 'v^3'};
%term_names = generate_term_names(3);
fprintf('数学表达式:\n');
fprintf('du/dt = ');
printed = false;
for i = 1:size(results,1)
    coeff = results(i,1);
    if abs(coeff) > 1e-3  % 忽略接近零的系数
        if printed
            if coeff > 0
                fprintf(' + ');
            else
                fprintf(' - ');
            end
            fprintf('%5.3f%s', abs(coeff), term_names{i});
        else
            printed = true;
            if coeff < 0
                fprintf('-');
            end
            fprintf('%5.3f%s', abs(coeff), term_names{i});
        end
    end
end
if ~printed
    fprintf('0');
end

fprintf('\n\ndv/dt = ');
printed = false;
for i = 1:size(results,1)
    coeff = results(i,2);
    if abs(coeff) > 1e-3  % 忽略接近零的系数
        if printed
            if coeff > 0
                fprintf(' + ');
            else
                fprintf(' - ');
            end
            fprintf('%5.3f%s', abs(coeff), term_names{i});
        else
            printed = true;
            if coeff < 0
                fprintf('-');
            end
            fprintf('%5.3f%s', abs(coeff), term_names{i});
        end
    end
end
if ~printed
    fprintf('0');
end
fprintf('\n');

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