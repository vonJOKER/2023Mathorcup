%% 二维函数可视化
function plot_2D_function(fun, lb, ub)
    % 生成网格点
    x = linspace(lb(1), ub(1), 100);
    y = linspace(lb(2), ub(2), 100);
    [X, Y] = meshgrid(x, y);
    Z = zeros(size(X));
    for i = 1:numel(X)
        Z(i) = fun([X(i), Y(i)]);
    end
    % 绘制等高线图和粒子
    contour(X, Y, Z, 20, 'LineWidth', 1);
    hold on
    scatter(lb(1), lb(2), 'g', 'filled')
    scatter(ub(1), ub(2), 'r', 'filled')
    scatter(gbest(1), gbest(2), 'k', 'filled')
    scatter(x(:, 1), x(:, 2), 20, 'b', 'filled')
    xlim([lb(1), ub(1)])
    ylim([lb(2), ub(2)])
    hold off
end