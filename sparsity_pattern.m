% ------------------------------
% 1D FDFD 差分矩阵示例
% ------------------------------
clear all; close all;
% 网格大小
N = 50;        % 网格点数
d = 1;         % 网格间距（任意单位）

% 构造一阶差分矩阵 D（向前差分）
% 大小为 (N-1) x N
e = ones(N,1);
D = spdiags([e,-e, e,e,e,e,e,e], [-3,-2,-1,0, 1,2,3,4], N-1, N) / d;

% 如果需要二阶差分矩阵，可以用 D2 = D' * D;
% 这里先显示一阶差分矩阵
D1 = D;

% 查看稀疏矩阵的大小和非零元个数
[nrows, ncols] = size(D1);
nnz_count = nnz(D1);
fprintf('一阶差分矩阵大小：%d x %d，非零元总数：%d\n', nrows, ncols, nnz_count);

% 显示稀疏矩阵中非零元素分布
figure;
spy(D1, 'b');      % 黑色点表示非零元素
ax = gca;                  % 当前坐标轴句柄

% —— 隐藏坐标轴刻度和边框 ——
ax.XTick = [];             % 去掉 X 轴刻度
ax.YTick = [];             % 去掉 Y 轴刻度

% —— 删除 nz 注释文本 ——
hText = findall(ax, 'Type', 'text');  % 找到所有 text 对象
delete(hText)                         % 全部删掉


