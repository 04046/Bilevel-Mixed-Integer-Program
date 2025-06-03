function  result  = main_enumeration()
%% An enumeration (Brute�\Force) method for solving the Bilevel Mixed Integer Program
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
% �㷨˼·��
% 1) ö�� x �����п���ȡֵ��0,1,2,��,x_max��
% 2) ����ÿ�������� x��ö�� z �����п���ȡֵ��0,1,2,��,z_max��
%    - �жϾ���Լ�� A*[x;z] <= b �Ƿ�����
%    - ������У������²�Ŀ�� f(z) = -z������ͬһ�� x ���ҳ��������� z��ʹ f(z) ���
% 3) ���ڸ� x ���������� z ��ѡ���������ϲ�Ŀ�� F(x,z) = -x - 10*z��ѡ��һ��ʹ���ϲ����ŵ� z
% 4) ����ȫ������ (x*, z*) �Լ���Ӧ���ϲ�Ŀ��ֵ

clear; clc;

%% 1. ����������Χ�;�������
x_max = 10;   % x ��ö������
z_max = 10;   % z ��ö������

% �ϲ�Ŀ��ϵ������������ֱ�������²�ö�٣�ֻ�����ռ��� F(x,z) ʱʹ�ã�
c_upper = [-1; -10];   % F(x,z) = c_upper' * [x;z] = -x -10*z

% �²�Ŀ��ϵ������������ f(z) = -z��ֻ�� z �йأ�
% Ϊ��ͳһд���������� [x;z] ���ʱ d'*[x;z] = -z
c_lower = [0; -1]; % f(z) = c_lower' * [x;z] = 0*x + (-1)*z = -z

% ͳһд�� A * [x;z] �� b ����ʽ������ [x;z] �� 2��1 ������
A = [ -25,   20;   % -25*x + 20*z <= 30
    +1,    2;   %   1*x +  2*z <= 10
    +2,   -1;   %   2*x -  1*z <= 15
    -2,  -10 ]; %  -2*x - 10*z <= -15  �൱�� 2*x +10*z >= 15
b = [ 30; 10; 15; -15 ];

% ��ʼ��ȫ�����Ž��¼
bestObj   = +Inf;   % ��¼�����ϲ�Ŀ��ֵ����Ϊ�ϲ��� minimization��
bestX     = NaN;    % ��¼���ŵ� x
bestZ     = NaN;    % ��¼��Ӧ������ z

%% 2. ���ѭ����ö�� x = 0,1,2,��,x_max
for x = 0 : x_max
    % �� x �̶����²����һ������ z ��������⣺
    %   max_{z}  f(z) = -z   s.t.   A * [x;z] <= b
    %
    % ���ڵ�ǰ x������Ҫ�ҵ�ʹ�²� f(z) �������� z ֵ����Ϊ�����ж�����Ž⣩��
    bestInnerObj = -Inf;  % ��¼��ǰ x �£��²�����Ŀ��ֵ max(-z)
    Z_optimalSet = [];    % ��¼���дﵽ�����²�ֵ�� z ��ѡ

    %% 2.1. �²�ö�٣�ö�� z = 0,1,2,��,z_max
    for z = 0 : z_max
        % �� x,z ��ϳ�һ�� 2��1 ����
        v = [x; z];
        % ����Ƿ����� A*v <= b
        if all( A * v <= b )
            % ������У�������²�Ŀ��
            innerObj = c_lower' * v;  % = -z
            % �����²�����ֵ�����Ž⼯��
            if innerObj > bestInnerObj
                bestInnerObj = innerObj;
                Z_optimalSet = z;
            elseif innerObj == bestInnerObj
                % ����²�Ŀ����ͬ��Ҳ���� z ���뼯��
                Z_optimalSet = [Z_optimalSet, z];
            end
        end
    end

    % ������� x �£��²�û�п��н⣬�������� x
    if isempty(Z_optimalSet)
        continue;
    end

    %% 2.2. �ϲ���ߣ��� Z_optimalSet ��ѡһ�� z��ʹ�ϲ� F(x,z) ��С
    localBestObj_x = +Inf;  % ��¼��ǰ x �£��ϲ�ľֲ�����ֵ
    localBestZ     = NaN;   % ��¼��ǰ x �£���Ӧʹ�ϲ����ŵ� z

    for z_cand = Z_optimalSet
        v_cand = [x; z_cand];
        F_val  = c_upper' * v_cand;  % = -x - 10*z_cand

        if F_val < localBestObj_x
            localBestObj_x = F_val;
            localBestZ     = z_cand;
        end
    end

    %% 2.3. ����ȫ�����Ž�
    if localBestObj_x < bestObj
        bestObj = localBestObj_x;
        bestX   = x;
        bestZ   = localBestZ;
    end
end

%% 3. ���ȫ�����Ž�
if isfinite(bestObj)
    fprintf('���Ž⣺x* = %d,  z* = %d\n', bestX, bestZ);
    fprintf('��Ӧ�������ϲ�Ŀ��ֵ��F(x*,z*) = %.2f\n', bestObj);
else
    disp('�ڸ�����Χ�ڣ�û���ҵ����н⡣');
end
end
