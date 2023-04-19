%% 粒子群算法求解Rosenbrock函数的最小值
function pso_rosenbrock
    % 粒子数目
    n_particles = 50;
    % 迭代次数
    n_iterations = 100;
    % 搜索空间的边界
    lb = [-2, -2];
    ub = [2, 2];
    % 粒子群算法参数设置
    c1 = 2.0; % 学习因子1
    c2 = 2.0; % 学习因子2
    w = 0.7; % 惯性权重
    vmax = 0.1 * (ub - lb); % 粒子速度的最大值
    n_dim = numel(lb); % 搜索空间的维度
    % 初始化粒子的位置和速度
    x = repmat(lb, n_particles, 1) + rand(n_particles, n_dim) .* repmat(ub - lb, n_particles, 1);
    v = rand(n_particles, n_dim) .* repmat(vmax, n_particles, 1);
    % 初始化粒子的历史最佳位置和全局最佳位置
    p = x;
    fp = rosenbrock(p');
    [fgbest, fgidx] = min(fp);
    gbest = p(fgidx, :);
    % 开始迭代
    for i = 1:n_iterations
        % 更新粒子速度和位置
        r1 = rand(n_particles, n_dim);
        r2 = rand(n_particles, n_dim);
        v = w * v + c1 * r1 .* (p - x) + c2 * r2 .* repmat(gbest, n_particles, 1) - c2 * r2 .* x;
        v = max(min(v, vmax), -vmax);
        x = x + v;
        x = max(min(x, ub), lb);
        % 更新粒子历史最佳位置和全局最佳位置
        fp = rosenbrock(p');
        fx = rosenbrock(x');
        idx = fx < fp;
        p(idx, :) = x(idx, :);
        fp(idx) = fx(idx);
        [fgbest, fgidx] = min(fp);
        gbest = p(fgidx, :);
        % 可视化
        plot_particles(x, gbest, lb, ub, rosenbrock)
        pause(0.01)
    end
    % 显示最优解和最优值
    fprintf('最优解为 (%f, %f)\n', gbest(1), gbest(2))
    fprintf('最优值为 %f\n', rosenbrock(gbest'))
end

%% Rosenbrock函数
function fval = rosenbrock(x)
    fval = sum(100 * (x(2:end) - x(1:end-1).^2).^2 + (1 - x(1:end-1)).^2, 1);
end