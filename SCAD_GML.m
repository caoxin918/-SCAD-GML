function [x] = SCAD_GML(G, b_analytic, nodesMat)
%UNTITLED6 此处显示有关此函数的摘要 该函数主要用于求解线性逆问题AX=B
% G: 系统矩阵A
% b_analystic: 节点能量矩阵  --   体表的荧光信号分布B
% nodesMat: 节点矩阵，前三列为节点的xyz坐标，最后一列为节点从属的区域序号
% 使用牛顿迭代法迭代求解X  --  体内的光源分布信息
% step1: Initialize the x and max_iterations
nums=size(nodesMat,1);
x0 = pinv(G)*b_analytic;      % 将x初始化为A的伪逆×B，有助于收敛
max_iterations = 600;         % 设置最大迭代次数为600次
L = computationL(nums, nodesMat);
[x] = newton(x0, max_iterations, G, b_analytic, nodesMat, nums, L);    % 调用函数输出迭代结果X
end

% step2: 定义fun1函数用来计算f(x)
function y = fun1(G, x, b_analystic, nums, L)
lambda1 = 0.1;
lambda2 = 0.01;
a = 3;
y = 0.5*(norm(G*x-b_analystic, 2)^2);
y = y + lambda1*((norm(L*x,2))^2);
% 计算scad项的值
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
function dy = fun2(G, x, b_analystic, nums, L)
lambda1 = 0.01;   
lambda2 = 0.01;
a = 3;
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

% 判断序号为i和j的节点是否属于同一个器官
function flag = isSameOrgan(i, j, indexArray)
% 通过在器官的节点集合里查找序号i，j来判断序号为i、j的节点是否属于同一个器官
flag = ((~isempty(find(indexArray==i))) && (~isempty(find(indexArray==j))));
end

% 计算psk---第k个器官的区域缓和器，由组成器官的不同节点计算得出
function pvalue = pskcal(nodes)
[m, ~] = size(nodes);
psk = 0;
R = 2; 
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

% 因为图流行学习正则化项的作用就是挖掘重建源的潜在流行结构，所以需要利用到模型的结构信息来计算
% 因此这里计算动态图拉普拉斯矩阵-L矩阵需要用到模型结构信息，需要判断节点是否属于同一器官组织（详见论文），所以计算时需要根据模型信息进行修改。
% 本算法使用到了图流行学习来挖掘重建源的潜在特征，提高算法的形态相似度，因此在计算时需要结合具体的模型来对算法进行修改。
% 具体来说，在计算图拉普拉斯矩阵的过程中需要用到模型的结构信息，也就是需要判断不同的节点是否属于同一个器官，这一块需要结合具体的模型来进行特定的处理，
% 这里是通过节点从属的边界序号判断节点从属的器官（可以从comsol软件中获取），进而判断出节点是否属于同一个器官组织。
function matrix_L = computationL(nums, nodesMat)
if exist('Step_scadM\L.mat','file')    % 判断文件是否存在
    disp('读取L矩阵！');
    load('Step_scadM\L.mat');
else
    % 计算L矩阵-动态图拉普拉斯矩阵（Manifold Learning）
    BoundaryIndex = nodesMat(:, 4);   % 提取出节点从属的边界序号 n*1
    nodes = nodesMat(:, 1:3);     
    R = 2;    
    index_boundary = find((BoundaryIndex==1)|(BoundaryIndex==2)|(BoundaryIndex==3)|(BoundaryIndex==4)|(BoundaryIndex==37)|(BoundaryIndex==42));
    index_heart = find((BoundaryIndex==25)|(BoundaryIndex==26)|(BoundaryIndex==27)|(BoundaryIndex==28)|(BoundaryIndex==33)|(BoundaryIndex==34)|(BoundaryIndex==35)|(BoundaryIndex==36));
    index_liver = find((BoundaryIndex==5)|(BoundaryIndex==6)|(BoundaryIndex==7)|(BoundaryIndex==8)|(BoundaryIndex==17)|(BoundaryIndex==18)|(BoundaryIndex==19)|(BoundaryIndex==20));
    index_bone = find((BoundaryIndex==38)|(BoundaryIndex==39)|(BoundaryIndex==40)|(BoundaryIndex==41)|(BoundaryIndex==43)|(BoundaryIndex==44));
    index_lungs1 = find((BoundaryIndex==13)|(BoundaryIndex==14)|(BoundaryIndex==15)|(BoundaryIndex==16)|(BoundaryIndex==29)|(BoundaryIndex==30)|(BoundaryIndex==31)|(BoundaryIndex==32));
    index_lungs2 = find((BoundaryIndex==9)|(BoundaryIndex==10)|(BoundaryIndex==11)|(BoundaryIndex==12)|(BoundaryIndex==21)|(BoundaryIndex==22)|(BoundaryIndex==23)|(BoundaryIndex==24));
    nodes_heart = nodesMat(index_heart, 1:3);  
    psk_heart = pskcal(nodes_heart);
    nodes_liver = nodesMat(index_liver, 1:3);   
    psk_liver = pskcal(nodes_liver);
    nodes_bone = nodesMat(index_bone, 1:3);   
    psk_bone = pskcal(nodes_bone);
    nodes_hung1 = nodesMat(index_lungs1, 1:3);   
    psk_hung1 = pskcal(nodes_hung1);
    nodes_hung2 = nodesMat(index_lungs2, 1:3);  
    psk_hung2 = pskcal(nodes_hung2);
    L = zeros(nums, nums);   % 初始化L矩阵为n*n的0矩阵
    for i = 1 : nums
        for j = 1 : nums
            if i == j
                L(i, j) = 1;   
            else
                % 判断两个节点是否同属于一个器官
                if  isSameOrgan(i, j, index_heart) || isSameOrgan(i, j, index_liver) || isSameOrgan(i, j, index_bone) || isSameOrgan(i, j, index_lungs1) || isSameOrgan(i, j, index_lungs2)          % i, j 节点都属于同一个器官Sk
                    dij = sqrt(sum((nodes(i, :) - nodes(j, :)).^2));   % 计算节点i与节点j之间的欧几里得距离
                    if isSameOrgan(i, j, index_heart)    % 节点ij同属于心脏组织的情况
                        L(i, j) = (-exp((-dij^2)/(4*R*R))) / psk_heart;
                    elseif isSameOrgan(i, j, index_liver)   
                        L(i, j) = (-exp((-dij^2)/(4*R*R))) / psk_liver;
                    elseif isSameOrgan(i, j, index_bone)    
                        L(i, j) = (-exp((-dij^2)/(4*R*R))) / psk_bone;
                    elseif isSameOrgan(i, j, index_lungs1)  
                        L(i, j) = (-exp((-dij^2)/(4*R*R))) / psk_hung1;
                    else    
                        L(i, j) = (-exp((-dij^2)/(4*R*R))) / psk_hung2;
                    end
                else
                    L(i, j) = 0;
                end
            end
        end
    end
    parSave('Step_scadM\L.mat','L', L);   % 保存L矩阵，避免重复计算
end
matrix_L = L;
end

% step4: 使用牛顿法迭代更新x
function [x]=newton(x0,max_iterations, G, b_analystic, nodesMat, nums, L) 
xk=x0;
k=1;	        % 保存迭代总次数 
while max_iterations >= k
    minf = fun1(G, xk, b_analystic, nums, L);
    minmu = 0.01;
    dfx=fun2(G, xk, b_analystic, nums, L);
    xk=xk-minmu*dfx;  % mu固定为一个数0.01
    xk = resetZero(xk, nodesMat);
    k=k+1;
end
parSave('Step_scadM\output.mat','output', xk);   % 保存output矩阵
x = xk;
end 