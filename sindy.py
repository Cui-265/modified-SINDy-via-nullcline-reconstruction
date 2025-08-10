import sys
sys.path.append("D:\\pycline_venv\\Lib\\site-packages")
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter
from scipy.linalg import lstsq

# 定义所有需要的函数

def poolData(yin, nVars, polyorder):
    """构建多项式特征库"""
    n = yin.shape[0]
    yout = np.ones((n, 1))  # 常数项
    yout = np.hstack((yout, yin))  # 一阶项
    
    # 二阶项
    if polyorder >= 2:
        yout = np.hstack((yout, yin[:, 0:1]**2))
        yout = np.hstack((yout, yin[:, 0:1] * yin[:, 1:2]))
        yout = np.hstack((yout, yin[:, 1:2]**2))
    
    # 三阶项
    if polyorder >= 3:
        yout = np.hstack((yout, yin[:, 0:1]**3))
        yout = np.hstack((yout, yin[:, 0:1]**2 * yin[:, 1:2]))
        yout = np.hstack((yout, yin[:, 0:1] * yin[:, 1:2]**2))
        yout = np.hstack((yout, yin[:, 1:2]**3))
    
    return yout

def sparsifyDynamics(Theta, dXdt, lambda_val, nVars):
    """使用阈值进行稀疏回归"""
    Xi = lstsq(Theta, dXdt)[0]
    
    for k in range(10):
        small_inds = np.abs(Xi) < lambda_val
        Xi[small_inds] = 0
        
        for j in range(nVars):
            big_inds = ~small_inds[:, j]
            if np.any(big_inds):
                Xi[big_inds, j] = lstsq(Theta[:, big_inds], dXdt[:, j])[0]
    
    return Xi

def central_difference(X, t):
    """改进的向量化中心差分实现"""
    n = X.shape[0]
    dX = np.zeros_like(X)
    
    # 边界处理
    if n > 1:
        dX[0] = (X[1] - X[0]) / (t[1] - t[0])
        dX[-1] = (X[-1] - X[-2]) / (t[-1] - t[-2])
    
    # 内部点（向量化）
    if n > 2:
        dt_forward = t[2:] - t[1:-1]
        dt_backward = t[1:-1] - t[:-2]
        total_dt = dt_backward + dt_forward
        
        total_dt[total_dt == 0] = np.finfo(float).eps
        
        term1 = (X[2:] - X[1:-1]) / dt_forward[:, np.newaxis]
        term2 = (X[1:-1] - X[:-2]) / dt_backward[:, np.newaxis]
        
        weight1 = dt_backward / total_dt
        weight2 = dt_forward / total_dt
        
        dX[1:-1] = weight1[:, np.newaxis] * term1 + weight2[:, np.newaxis] * term2
    
    return dX

def prepare_sindy_data(trajectories, time_vectors, derivative_method='spline', smooth_derivatives=True):
    """整合多条轨迹数据用于SINDy算法"""
    if len(trajectories) != len(time_vectors):
        raise ValueError("轨迹数据和时间向量数量不一致")
    
    X_stack = []
    dX_stack = []
    
    for X, t in zip(trajectories, time_vectors):
        if len(X) != len(t):
            raise ValueError(f"状态数据和时间向量长度不一致: {len(X)} != {len(t)}")
        
        if derivative_method == 'spline':
            dX = np.zeros_like(X)
            for j in range(X.shape[1]):
                cs = CubicSpline(t, X[:, j])
                dX[:, j] = cs(t, 1)
                
        elif derivative_method == 'savitzky-golay':
            dX = np.zeros_like(X)
            n = len(t)
            order = min(3, n-1)
            framelen = min(21, n)
            if framelen % 2 == 0:
                framelen = max(3, framelen - 1)
                
            for j in range(X.shape[1]):
                dX[:, j] = savgol_filter(X[:, j], framelen, order, deriv=1, delta=np.mean(np.diff(t)))
                
        else:
            dX = central_difference(X, t)
        
        if smooth_derivatives and len(t) > 5:
            window_size = min(11, len(t)//2)
            if window_size % 2 == 0:
                window_size = max(3, window_size - 1)
            dX = savgol_filter(dX, window_size, 3, axis=0)
        
        X_stack.append(X)
        dX_stack.append(dX)
    
    X_stack = np.vstack(X_stack)
    dX_stack = np.vstack(dX_stack)
    
    print(f"成功整合{len(trajectories)}条轨迹数据")
    print(f"总数据点数: {X_stack.shape[0]}, 状态变量数: {X_stack.shape[1]}")
    
    return X_stack, dX_stack

def sindy_model_identification_k_plot(X, dt, dXdt, offsets, threshold=0.08, poly_order=3):
    """从动力系统数据识别模型并计算K值分布"""
    du_offsets = offsets[:, 0]
    dv_offsets = offsets[:, 1]
    
    u_data = X[:, 0]
    v_data = X[:, 1]
    
    all_offsets = []
    all_K = []
    
    for i, du_offset in enumerate(du_offsets):
        for j, dv_offset in enumerate(dv_offsets):
            u_shifted = u_data + du_offset
            v_shifted = v_data + dv_offset
            
            Theta = poolData(np.column_stack((u_shifted, v_shifted)), 2, poly_order)
            Xi = sparsifyDynamics(Theta, dXdt, threshold, 2)
            k = np.count_nonzero(Xi)
            
            all_offsets.append([du_offset, dv_offset])
            all_K.append(k)
    
    all_offsets = np.array(all_offsets)
    all_K = np.array(all_K)
    
    du_vals = np.unique(du_offsets)
    dv_vals = np.unique(dv_offsets)
    DU, DV = np.meshgrid(du_vals, dv_vals)
    K_grid = all_K.reshape(len(du_vals), len(dv_vals)).T
    
    return DU, DV, K_grid

# 定义原始系统和近似系统
def original_system(t, z):
    u, v = z
    du = -u**3 + u**2 + u - v
    dv = -0.5*v**3 + 0.5*v**2 - (1/3)*v + u
    return [du, dv]

def approx_system1(t, z):
    u, v = z
    du = 0.058 + 0.834*u - 0.906*v + 0.801*u**2 + 0.086*u*v - 0.081*v**2 - 0.8*u**3
    dv = 1.003*u - 0.364*v + 0.467*v**2 - 0.455*v**3
    return [du, dv]

def approx_system2(t, z):
    u, v = z
    du = 0.074 + 0.882*u - 0.879*v + 0.792*u**2 - 0.812*u**3 - 0.083*v**2
    dv = 1.003*u - 0.364*v + 0.467*v**2 - 0.455*v**3
    return [du, dv]

# 主程序
if __name__ == "__main__":
    # 设置参数
    dt = 0.01
    tspan = [0, 70]
    t_eval = np.linspace(0, 70, 1000)
    z0 = [0.5, 0.5]
    
    # 求解轨迹
    sol_orig = solve_ivp(original_system, tspan, z0, t_eval=t_eval, method='RK45', rtol=1e-6, atol=1e-9)
    sol_app1 = solve_ivp(approx_system1, tspan, z0, t_eval=t_eval, method='RK45', rtol=1e-6, atol=1e-9)
    sol_app2 = solve_ivp(approx_system2, tspan, z0, t_eval=t_eval, method='RK45', rtol=1e-6, atol=1e-9)
    
    # 计算误差
    error_app1 = np.sqrt((sol_orig.y[0] - sol_app1.y[0])**2 + (sol_orig.y[1] - sol_app1.y[1])**2)
    error_app2 = np.sqrt((sol_orig.y[0] - sol_app2.y[0])**2 + (sol_orig.y[1] - sol_app2.y[1])**2)
    cumulative_error1 = np.cumsum(error_app1) * (t_eval[1]-t_eval[0])
    cumulative_error2 = np.cumsum(error_app2) * (t_eval[1]-t_eval[0])
    
    # 计算Nullcline
    u = np.linspace(-0.5, 1.5, 500)
    v = np.linspace(-0.5, 1.5, 500)
    U, V = np.meshgrid(u, v)
    
    du_orig = -U**3 + U**2 + U - V
    dv_orig = -0.5*V**3 + 0.5*V**2 - (1/3)*V + U
    du_app1 = 0.058 + 0.834*U - 0.906*V + 0.801*U**2 + 0.086*U*V - 0.081*V**2 - 0.8*U**3
    dv_app1 = 1.003*U - 0.364*V + 0.467*V**2 - 0.455*V**3
    du_app2 = 0.074 + 0.882*U - 0.879*V + 0.792*U**2 - 0.812*U**3 - 0.083*V**2
    dv_app2 = 1.003*U - 0.364*V + 0.467*V**2 - 0.455*V**3
    
    # 准备SINDy数据
    tspan_k = [0, 40]
    t_eval_k = np.arange(0, 40+dt, dt)
    u0 = [0.1, 0.1]
    u0_1 = [-4, 0.5]
    u0_2 = [1, -3]
    u0_3 = [3, -3.5]
    
    sol = solve_ivp(original_system, tspan_k, u0, t_eval=t_eval_k, method='LSODA', rtol=1e-6, atol=1e-9)
    sol1 = solve_ivp(original_system, tspan_k, u0_1, t_eval=t_eval_k, method='LSODA', rtol=1e-6, atol=1e-9)
    sol2 = solve_ivp(original_system, tspan_k, u0_2, t_eval=t_eval_k, method='LSODA', rtol=1e-6, atol=1e-9)
    sol3 = solve_ivp(original_system, tspan_k, u0_3, t_eval=t_eval_k, method='LSODA', rtol=1e-6, atol=1e-9)
    
    trajectories = [
        np.column_stack((sol.y[0], sol.y[1])),
        np.column_stack((sol1.y[0], sol1.y[1])),
        np.column_stack((sol2.y[0], sol2.y[1])),
        np.column_stack((sol3.y[0], sol3.y[1]))
    ]
    
    time_vectors = [sol.t, sol1.t, sol2.t, sol3.t]
    
    stacked_states, stacked_derivatives = prepare_sindy_data(
        trajectories, time_vectors, 
        derivative_method='spline',
        smooth_derivatives=True
    )
    
    custom_offsets = np.column_stack((
        np.arange(-2, 2.1, 0.1),
        np.arange(-2, 2.1, 0.1)
    ))
    
    # 计算K值分布
    DU_k, DV_k, K_grid = sindy_model_identification_k_plot(
        X=stacked_states,
        dt=dt,
        dXdt=stacked_derivatives,
        offsets=custom_offsets,
        threshold=0.1,
        poly_order=3
    )
    
    # 创建2x2图表布局
    plt.figure(figsize=(18, 12))
    
    # 图1: 瞬时误差随时间变化
    plt.subplot(2, 2, 1)
    plt.plot(t_eval, error_app1, 'r-', linewidth=2, label='Approx System(normal)')
    plt.plot(t_eval, error_app2, 'g-', linewidth=2, label='Approx System(optimized)')
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Instantaneous Error', fontsize=12)
    plt.title('(a) Instantaneous Trajectory Error', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlim(0, 70)
    plt.ylim(0, 1.2)
    plt.text(35, 0.22, f'Final Error (t=70):\nnormal: {error_app1[-1]:.4f}\noptimized: {error_app2[-1]:.4f}',
             fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
    
    # 图2: 累积误差随时间变化
    plt.subplot(2, 2, 2)
    plt.plot(t_eval, cumulative_error1, 'r-', linewidth=2, label='Approx System(normal)')
    plt.plot(t_eval, cumulative_error2, 'g-', linewidth=2, label='Approx System(optimized)')
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Cumulative Error', fontsize=12)
    plt.title('(b) Cumulative Trajectory Error', fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlim(0, 70)
    max_cumulative = max(cumulative_error1.max(), cumulative_error2.max())
    plt.ylim(0, max_cumulative * 1.1)
    plt.text(35, max_cumulative * 0.8, 
             f'Total Error (t=70):\nnormal: {cumulative_error1[-1]:.4f}\noptimized: {cumulative_error2[-1]:.4f}',
             fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
    
    # 图3: Nullcline比较
    plt.subplot(2, 2, 3)
    # 原系统 nullcline
    plt.contour(U, V, du_orig, levels=[0], colors='blue', linestyles='dashed', linewidths=2.5, alpha=0.9)
    plt.contour(U, V, dv_orig, levels=[0], colors='blue', linestyles='solid', linewidths=2.5, alpha=0.9)
    # 近似系统1 nullcline
    plt.contour(U, V, du_app1, levels=[0], colors='red', linestyles='dashed', linewidths=1.8, alpha=0.7)
    plt.contour(U, V, dv_app1, levels=[0], colors='red', linestyles='solid', linewidths=1.8, alpha=0.7)
    # 近似系统2 nullcline
    plt.contour(U, V, du_app2, levels=[0], colors='green', linestyles='dashed', linewidths=1.8, alpha=0.7)
    plt.contour(U, V, dv_app2, levels=[0], colors='green', linestyles='solid', linewidths=1.8, alpha=0.7)
    
    # 创建自定义图例
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='blue', ls='--', lw=2.5, label='Original: u-nullcline'),
        Line2D([0], [0], color='blue', ls='-', lw=2.5, label='Original: v-nullcline'),
        Line2D([0], [0], color='red', ls='--', lw=1.8, label='Normal: u-nullcline'),
        Line2D([0], [0], color='red', ls='-', lw=1.8, label='Normal: v-nullcline'),
        Line2D([0], [0], color='green', ls='--', lw=1.8, label='Optimized: u-nullcline'),
        Line2D([0], [0], color='green', ls='-', lw=1.8, label='Optimized: v-nullcline')
    ]
    
    plt.xlabel('u', fontsize=12)
    plt.ylabel('v', fontsize=12)
    plt.title('(c) Nullcline Comparison', fontsize=14)
    plt.legend(handles=legend_elements, loc='upper right', fontsize=9)
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.xlim(-0.5, 1.5)
    plt.ylim(-0.5, 1.5)
    
    # 图4: K值分布图
    plt.subplot(2, 2, 4)
    contour = plt.contourf(DU_k, DV_k, K_grid, 20, cmap='viridis')
    plt.colorbar(contour, label='K value(complexity)')
    plt.xlabel('Δu', fontsize=12)
    plt.ylabel('Δv', fontsize=12)
    plt.title('(d) K Chart(complexity)', fontsize=14)
    plt.plot(0, 0, 'yo', markersize=10, label='real fixed point (0,0)')
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('comprehensive_analysis.png', dpi=300)
    plt.show()