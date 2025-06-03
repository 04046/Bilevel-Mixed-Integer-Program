function  result  = main_enumeration()
%% An enumeration (Brute‐Force) method for solving the Bilevel Mixed Integer Program
% Date: 2025/06/03
% Author: XY

% Example 1 from Ref.[1] Bo Zeng and Yu An, Solving Bilevel Mixed Integer
% Program by Reformulations and Decomposition, 2014.

% Example 1:
%      min -x - 10*z
% s.t. x \in Z_+
%      z \in argmax{-z:
%      -25*x + 20*z <=30,
%      x + 2*z <=10.
%      2*x - z <=15,
%      2*x + 10*z >= 15, x \in Z_+}
% -------------------------------------------------------------------------
% 算法思路：
% 1) 枚举 x 的所有可能取值（0,1,2,…,x_max）
% 2) 对于每个给定的 x，枚举 z 的所有可能取值（0,1,2,…,z_max）
%    - 判断矩阵约束 A*[x;z] <= b 是否满足
%    - 如果可行，计算下层目标 f(z) = -z，并在同一个 x 下找出所有最优 z（使 f(z) 最大）
% 3) 对于该 x 下所有最优 z 候选集，计算上层目标 F(x,z) = -x - 10*z，选出一个使得上层最优的 z
% 4) 更新全局最优 (x*, z*) 以及对应的上层目标值

clear; clc;

%% 1. 定义搜索范围和矩阵数据
x_max = 10;   % x 的枚举上限
z_max = 10;   % z 的枚举上限

% 上层目标系数向量（并不直接用于下层枚举，只在最终计算 F(x,z) 时使用）
c_upper = [-1; -10];   % F(x,z) = c_upper' * [x;z] = -x -10*z

% 下层目标系数向量（由于 f(z) = -z，只与 z 有关）
% 为了统一写成向量，与 [x;z] 相乘时 d'*[x;z] = -z
c_lower = [0; -1]; % f(z) = c_lower' * [x;z] = 0*x + (-1)*z = -z

% 统一写成 A * [x;z] ≤ b 的形式，其中 [x;z] 是 2×1 列向量
A = [ -25,   20;   % -25*x + 20*z <= 30
    +1,    2;   %   1*x +  2*z <= 10
    +2,   -1;   %   2*x -  1*z <= 15
    -2,  -10 ]; %  -2*x - 10*z <= -15  相当于 2*x +10*z >= 15
b = [ 30; 10; 15; -15 ];

% 初始化全局最优解记录
bestObj   = +Inf;   % 记录最优上层目标值（因为上层是 minimization）
bestX     = NaN;    % 记录最优的 x
bestZ     = NaN;    % 记录对应的最优 z

%% 2. 外层循环：枚举 x = 0,1,2,…,x_max
for x = 0 : x_max
    % 将 x 固定后，下层就是一个关于 z 的最大化问题：
    %   max_{z}  f(z) = -z   s.t.   A * [x;z] <= b
    %
    % 对于当前 x，我们要找到使下层 f(z) 最大的所有 z 值（因为可能有多个最优解）。
    bestInnerObj = -Inf;  % 记录当前 x 下，下层最优目标值 max(-z)
    Z_optimalSet = [];    % 记录所有达到最优下层值的 z 候选

    %% 2.1. 下层枚举：枚举 z = 0,1,2,…,z_max
    for z = 0 : z_max
        % 将 x,z 组合成一个 2×1 向量
        v = [x; z];
        % 检查是否满足 A*v <= b
        if all( A * v <= b )
            % 如果可行，则计算下层目标
            innerObj = c_lower' * v;  % = -z
            % 更新下层最优值及最优解集合
            if innerObj > bestInnerObj
                bestInnerObj = innerObj;
                Z_optimalSet = z;
            elseif innerObj == bestInnerObj
                % 如果下层目标相同，也将该 z 加入集合
                Z_optimalSet = [Z_optimalSet, z];
            end
        end
    end

    % 如果给定 x 下，下层没有可行解，则跳过该 x
    if isempty(Z_optimalSet)
        continue;
    end

    %% 2.2. 上层决策：在 Z_optimalSet 中选一个 z，使上层 F(x,z) 最小
    localBestObj_x = +Inf;  % 记录当前 x 下，上层的局部最优值
    localBestZ     = NaN;   % 记录当前 x 下，对应使上层最优的 z

    for z_cand = Z_optimalSet
        v_cand = [x; z_cand];
        F_val  = c_upper' * v_cand;  % = -x - 10*z_cand

        if F_val < localBestObj_x
            localBestObj_x = F_val;
            localBestZ     = z_cand;
        end
    end

    %% 2.3. 更新全局最优解
    if localBestObj_x < bestObj
        bestObj = localBestObj_x;
        bestX   = x;
        bestZ   = localBestZ;
    end
end

%% 3. 输出全局最优解
if isfinite(bestObj)
    fprintf('最优解：x* = %d,  z* = %d\n', bestX, bestZ);
    fprintf('对应的最优上层目标值：F(x*,z*) = %.2f\n', bestObj);
else
    disp('在给定范围内，没有找到可行解。');
end
end
