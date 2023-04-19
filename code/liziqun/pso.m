function pso
    % 目标函数
    fun = @rosenbrock;
    % 定义搜索空间
    lb = [-2, -2];
    ub = [2, 2];
    % 初始化粒子群
    n_particles = 30;
    pos = rand(n_particles, 2) .* (ub - lb) + lb;
    vel = rand(n_particles, 2) .* (ub - lb) / 2;
    pbest = pos;
    pbest_fitness = arrayfun(fun, pos);
    [gbest_fitness, gbest_index] = min(pbest_fitness);
    gbest = pbest(gbest_index, :);
    % 参数设置
    max_iter = 100;
    w = 0.8;
    c1 = 2;
    c2 = 2;
    % 迭代优化
    for i = 1:max_iter
        % 更新速度和位置
        vel = w * vel + c1 * rand(n_particles, 2) .* (pbest - pos) + c2 * rand(n_particles, 2) .* (gbest - pos);
        pos = pos + vel;
        % 约束处理
        pos = max(pos, lb);
        pos = min(pos, ub);
        % 更新个体最优和全局最优
        fitness = arrayfun(fun, pos);
        update_index = fitness < pbest_fitness;
        pbest(update_index, :) = pos(update_index, :);
        pbest_fitness(update_index) = fitness(update_index);
        [min_fitness, min_index] = min(pbest_fitness);
        if min_fitness < gbest_fitness
            gbest = pbest(min_index, :);
            gbest_fitness = min_fitness;
        end
        % 显示粒子群状态
        fprintf('Iteration %d, Best Fitness = %.4f\n', i, gbest_fitness);
        plot_2D_function(fun, lb, ub)
        pause(0.1)
    end
end