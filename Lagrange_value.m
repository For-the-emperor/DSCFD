% MATLAB 脚本：计算半步偏移一阶 Lagrange δ 卷积核系数
clear; clc;

for M = 1:4
    fprintf('==== M = %d ====\n', M);
    % 构造节点位置 xi = r + 0.5, r = -M:(M-1)
    r = -M:(M-1);
    xi = r + 0.5;
    
    % 对每个节点 xi_j，构造 Lagrange 基函数，求导并在 xi=0 处评估
    for j = 1:length(xi)
        xj = xi(j);
        
        % 构造 Lagrange 基函数 Lj(x)
        % Lj(x) = prod_{k != j} (x - xi(k)) / (xi(j) - xi(k))
        syms x
        Lj = sym(1);
        for k = 1:length(xi)
            if k ~= j
                Lj = Lj * (x - xi(k)) / (xj - xi(k));
            end
        end
        
        % 计算一阶导数并在 x=0 处评估
        dLj = diff(Lj, x);
        w = double( subs(dLj, x, 0) );
        
        % 输出结果
        fprintf('r = %+d: xi = %+4.2f, weight = %9.6f\n', r(j), xj, w);
    end
    fprintf('\n');
end
