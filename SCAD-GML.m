function [x,k] = scadManifold(nums,G, b_analytic, nodesMat)
%UNTITLED6 此处显示有关此函数的摘要
% 函数共有三个参数，分别是
% nums: 节点个数  --  n = 4626
% G: 是一个598 * 4626的矩阵，其中598是表面的能量   --  系统矩阵A
% b_analystic: 是一个598 * 1的矩阵，是序号从属于外边界的节点能量矩阵  --   体表的荧光信号分布b
% nodesMat: 节点矩阵，前三列为节点的xyz坐标，最后一列为节点从属的区域序号1-44
% 使用牛顿迭代法迭代求解x  --  体内的光源分布信息

% step1: Initialize the x and max_iterations
tic;
%x0 = zeros(nums, 1);       % 将x初始化为一个N*1的零矩阵
x0 = pinv(G)*b_analytic;       % 将x初始化为A的伪逆×B，有助于收敛
%x0 = resetZero(x0, nodesMat);
max_iterations = 300;     % 设置最大迭代次数为300次 范围：300-1000
L = computationL(nums, nodesMat);
[x, k] = newton(x0, max_iterations, G, b_analytic, nodesMat, nums, L);    % 调用函数输出迭代结果X，以及迭代次数k
toc;
end

% step2: 定义fun1函数用来计算f(x)
function y = fun1(G, x, b_analystic, nodesMat, nums, L)
%% 正则化参数的初始化（论文里找到）
lambda1 = 0.1;  % lambda1 范围：0.1-0.001
lambda2 = 0.01; % lambda2 范围：0.1-0.0001
a = 3;
y = 0.5*(norm(G*x-b_analystic, 2)^2);
y = y + lambda1*((norm(L*x,2))^2);
% 计算scad正则化项的值
svalue = 0;
for i = 1 : nums
    if abs(x(i)) >= 0 && abs(x(i)) < lambda2
        func = lambda2 * abs(x(i));
    elseif abs(x(i)) >= lambda2 && abs(x(i)) < a * lambda2
        func = (x(i)^2 - 2*a*lambda2*abs(x(i)) + lambda2^2) / (2*(1-a));
    else
        func = ((a+1)*(lambda2^2)) / 2;
    end
    svalue = svalue + func;
end
y = y + svalue;
end

% step3: 定义fun2函数用来计算f'(x) 即f(x)的梯度
function dy = fun2(G, x, b_analystic, nodesMat, nums, L)
lambda1 = 0.01;   % 取值：0.5-0.001
lambda2 = 0.01;  % 取值：0.1-0.001 越小越稳定取0.01
a = 3;   % 
dy = (G')*(G*x - b_analystic) + 2 * lambda1 * ((L)'*L*x);
% g矩阵n*1存储函数g的偏导数（梯度）
g = zeros(nums,1);
for i = 1 : nums
    if abs(x(i)) > 0 && abs(x(i)) < lambda2
        func = lambda2 * sign(x(i));
    elseif abs(x(i)) >= lambda2 && abs(x(i)) < a * lambda2
        func = (a*lambda2*sign(x(i)) - x(i)) / (a-1);
    else
        func = 0;
    end
    g(i) = func;
end
dy = dy + g;
end

% 计算目标函数的海森矩阵（常矩阵，和x无关）
function matrix_H = computationH(G, L, nums)
a = 3.7;
lambda1 = 0.01;
matrix_scad = zeros(nums, nums);
for i = 1:nums
    matrix_scad(i, i) = 1 / (a - 1);
end
matrix_H = G'*G + 2*lambda1*L'*L + matrix_scad;
end

% 判断序号为i和j的节点是否属于同一个器官
function flag = isSameOrgan(i, j, indexArray)
% 通过在器官的节点集合里查找序号i，j来判断序号为i、j的节点是否属于同一个器官
flag = ((~isempty(find(indexArray==i))) && (~isempty(find(indexArray==j))));
end

% 计算psk---第k个器官的区域缓和器，由组成器官的不同节点计算得出
function pvalue = pskcal(nodes)
[m, ~] = size(nodes);
psk = 0;
R = 2;   % 高斯核函数半径
for i = 1 : m
    for j = i + 1 : m
        dij = sqrt(sum((nodes(i, :) - nodes(j, :)).^ 2));   % 计算mn节点之间的欧氏距离
        psk = psk + exp((-dij*dij)/(4*R*R));   
    end
end
pvalue = psk * 2;
end

% 将表面节点归零
function xk = resetZero(x, nodesMat)
BoundaryIndex = nodesMat(:, 4);   % 提取出节点从属的边界序号 n*1
index_boundary = find((BoundaryIndex==1)|(BoundaryIndex==2)|(BoundaryIndex==3)|(BoundaryIndex==4)|(BoundaryIndex==37)|(BoundaryIndex==42));
[num,~] = size(index_boundary);
for i = 1 : num
    x(index_boundary(i)) = 0;
end
xk = x;
end

function matrix_L = computationL(nums, nodesMat)
if exist('Step_scadM\L.mat','file')    % 判断文件是否存在
    disp('读取L矩阵！');
    load('Step_scadM\L.mat');
else
    %% 计算L矩阵-动态图拉普拉斯矩阵（Manifold Learning）
    % comsol软件里获得组织节点的编号（序号）  
    % 1.最外层边界不算是器官   2.左肺、右肺属于不同的器官，是两个独立的器官
    BoundaryIndex = nodesMat(:, 4);   % 提取出节点从属的边界序号 n*1
    nodes = nodesMat(:, 1:3);     % 提取出节点坐标矩阵 n*3
    R = 2;    % R 初始化为2
    index_boundary = find((BoundaryIndex==1)|(BoundaryIndex==2)|(BoundaryIndex==3)|(BoundaryIndex==4)|(BoundaryIndex==37)|(BoundaryIndex==42));
    index_heart = find((BoundaryIndex==25)|(BoundaryIndex==26)|(BoundaryIndex==27)|(BoundaryIndex==28)|(BoundaryIndex==33)|(BoundaryIndex==34)|(BoundaryIndex==35)|(BoundaryIndex==36));
    index_liver = find((BoundaryIndex==5)|(BoundaryIndex==6)|(BoundaryIndex==7)|(BoundaryIndex==8)|(BoundaryIndex==17)|(BoundaryIndex==18)|(BoundaryIndex==19)|(BoundaryIndex==20));
    index_bone = find((BoundaryIndex==38)|(BoundaryIndex==39)|(BoundaryIndex==40)|(BoundaryIndex==41)|(BoundaryIndex==43)|(BoundaryIndex==44));
    index_lungs1 = find((BoundaryIndex==13)|(BoundaryIndex==14)|(BoundaryIndex==15)|(BoundaryIndex==16)|(BoundaryIndex==29)|(BoundaryIndex==30)|(BoundaryIndex==31)|(BoundaryIndex==32));
    index_lungs2 = find((BoundaryIndex==9)|(BoundaryIndex==10)|(BoundaryIndex==11)|(BoundaryIndex==12)|(BoundaryIndex==21)|(BoundaryIndex==22)|(BoundaryIndex==23)|(BoundaryIndex==24));
    nodes_heart = nodesMat(index_heart, 1:3);   % 提取出从属于心脏组织的节点坐标集合
    psk_heart = pskcal(nodes_heart);
    nodes_liver = nodesMat(index_liver, 1:3);   % 提取出从属于肝脏器官的节点集合
    psk_liver = pskcal(nodes_liver);
    nodes_bone = nodesMat(index_bone, 1:3);   % 提取出从属于骨骼器官的节点集合
    psk_bone = pskcal(nodes_bone);
    nodes_hung1 = nodesMat(index_lungs1, 1:3);   % 提取出从属于左肺器官的节点集合
    psk_hung1 = pskcal(nodes_hung1);
    nodes_hung2 = nodesMat(index_lungs2, 1:3);   % 提取出从属于右肺器官的节点集合
    psk_hung2 = pskcal(nodes_hung2);
    L = zeros(nums, nums);   % 初始化L矩阵为n*n的0矩阵
    for i = 1 : nums
        for j = 1 : nums
            if i == j
                L(i, j) = 1;    % 公式更正为1，不然求偏导太麻烦
            else
                % 判断两个节点是否同属于一个器官
                if  isSameOrgan(i, j, index_heart) || isSameOrgan(i, j, index_liver) || isSameOrgan(i, j, index_bone) || isSameOrgan(i, j, index_lungs1) || isSameOrgan(i, j, index_lungs2)          % i, j 节点都属于同一个器官Sk
                    dij = sqrt(sum((nodes(i, :) - nodes(j, :)).^2));   % 计算节点i与节点j之间的欧几里得距离
                    if isSameOrgan(i, j, index_heart)    % 节点ij同属于心脏组织的情况
                        L(i, j) = (-exp((-dij^2)/(4*R*R))) / psk_heart;
                    elseif isSameOrgan(i, j, index_liver)    % 节点ij同属于肝脏组织的情况
                        L(i, j) = (-exp((-dij^2)/(4*R*R))) / psk_liver;
                    elseif isSameOrgan(i, j, index_bone)    % 节点ij同属于骨骼组织的情况
                        L(i, j) = (-exp((-dij^2)/(4*R*R))) / psk_bone;
                    elseif isSameOrgan(i, j, index_lungs1)    % 节点ij同属于肺1组织的情况
                        L(i, j) = (-exp((-dij^2)/(4*R*R))) / psk_hung1;
                    else    % 节点ij同属于肺2组织的情况
                        L(i, j) = (-exp((-dij^2)/(4*R*R))) / psk_hung2;
                    end
                else
                    L(i, j) = 0;
                end
            end
        end
    end
    disp('保存L矩阵！');
    parSave('Step_scadM\L.mat','L', L);   % 保存L矩阵，避免重复计算
end
matrix_L = L;
end

% step4: 使用牛顿法迭代更新x
function [x,k]=newton(x0,max_iterations, G, b_analystic, nodesMat, nums, L) 
% x0是迭代初始矩阵
% fun1是用来求取f(x)的值的函数
% fun2是用来求取f'(x)的值的函数
% max_iterations是最大迭代次数
xk=x0;
%y=fun1(x);      % 计算f(x)的值
k=1;	        % 保存迭代总次数 
H = computationH(G, L, nums);
while max_iterations >= k
    minf = fun1(G, xk, b_analystic, nodesMat, nums, L);
    minmu = 0.01; % 下山系数，范围：0.1-0.0001
    dfx=fun2(G, xk, b_analystic, nodesMat, nums, L);
%     mu = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0002, 0.0001];
%     [m, n] = size(mu);
%     for i = 1:n
%         x1 = xk;
%         %x1 = x1 - mu(i) * (H\dfx);
%         x1 = x1 - mu(i) * dfx;
%         fx = fun1(G, x1, b_analystic, nodesMat, nums, L);
%         if fx < minf
%             minf = fx;
%             minmu = mu(i);
%         end
%     end
    xk=xk-minmu*dfx;  % mu固定为一个数0.01
    %xk = xk - minmu*(H\dfx);
    xk = resetZero(xk, nodesMat);
    s = sprintf('这是第%d次循环，梯度的模长为%.6f，步长mu为%.6f，误差为%.6f!', k, sqrt(sum(dfx.*dfx)), minmu, norm(G*xk-b_analystic, 2));
    disp(s);
    k=k+1;
end
disp('保存output矩阵！');
parSave('Step_scadM\output.mat','output', xk);   % 保存output矩阵
x = xk;
end 