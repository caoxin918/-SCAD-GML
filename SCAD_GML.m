function [x] = SCAD_GML(G, b_analytic, nodesMat)
%UNTITLED6 �˴���ʾ�йش˺�����ժҪ �ú�����Ҫ�����������������AX=B
% G: ϵͳ����A
% b_analystic: �ڵ���������  --   ����ӫ���źŷֲ�B
% nodesMat: �ڵ����ǰ����Ϊ�ڵ��xyz���꣬���һ��Ϊ�ڵ�������������
% ʹ��ţ�ٵ������������X  --  ���ڵĹ�Դ�ֲ���Ϣ
% step1: Initialize the x and max_iterations
nums=size(nodesMat,1);
x0 = pinv(G)*b_analytic;      % ��x��ʼ��ΪA��α���B������������
max_iterations = 600;         % ��������������Ϊ600��
L = computationL(nums, nodesMat);
[x] = newton(x0, max_iterations, G, b_analytic, nodesMat, nums, L);    % ���ú�������������X
end

% step2: ����fun1������������f(x)
function y = fun1(G, x, b_analystic, nums, L)
lambda1 = 0.1;
lambda2 = 0.01;
a = 3;
y = 0.5*(norm(G*x-b_analystic, 2)^2);
y = y + lambda1*((norm(L*x,2))^2);
% ����scad���ֵ
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

% step3: ����fun2������������f'(x) ��f(x)���ݶ�
function dy = fun2(G, x, b_analystic, nums, L)
lambda1 = 0.01;   
lambda2 = 0.01;
a = 3;
dy = (G')*(G*x - b_analystic) + 2 * lambda1 * ((L)'*L*x);
% g����n*1�洢����g��ƫ�������ݶȣ�
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

% �ж����Ϊi��j�Ľڵ��Ƿ�����ͬһ������
function flag = isSameOrgan(i, j, indexArray)
% ͨ�������ٵĽڵ㼯����������i��j���ж����Ϊi��j�Ľڵ��Ƿ�����ͬһ������
flag = ((~isempty(find(indexArray==i))) && (~isempty(find(indexArray==j))));
end

% ����psk---��k�����ٵ����򻺺�������������ٵĲ�ͬ�ڵ����ó�
function pvalue = pskcal(nodes)
[m, ~] = size(nodes);
psk = 0;
R = 2; 
for i = 1 : m
    for j = i + 1 : m
        dij = sqrt(sum((nodes(i, :) - nodes(j, :)).^ 2));   % ����mn�ڵ�֮���ŷ�Ͼ���
        psk = psk + exp((-dij*dij)/(4*R*R));   
    end
end
pvalue = psk * 2;
end

% ������ڵ����
function xk = resetZero(x, nodesMat)
BoundaryIndex = nodesMat(:, 4);   % ��ȡ���ڵ�����ı߽���� n*1
index_boundary = find((BoundaryIndex==1)|(BoundaryIndex==2)|(BoundaryIndex==3)|(BoundaryIndex==4)|(BoundaryIndex==37)|(BoundaryIndex==42));
[num,~] = size(index_boundary);
for i = 1 : num
    x(index_boundary(i)) = 0;
end
xk = x;
end

% ��Ϊͼ����ѧϰ����������þ����ھ��ؽ�Դ��Ǳ�����нṹ��������Ҫ���õ�ģ�͵Ľṹ��Ϣ������
% ���������㶯̬ͼ������˹����-L������Ҫ�õ�ģ�ͽṹ��Ϣ����Ҫ�жϽڵ��Ƿ�����ͬһ������֯��������ģ������Լ���ʱ��Ҫ����ģ����Ϣ�����޸ġ�
% ���㷨ʹ�õ���ͼ����ѧϰ���ھ��ؽ�Դ��Ǳ������������㷨����̬���ƶȣ�����ڼ���ʱ��Ҫ��Ͼ����ģ�������㷨�����޸ġ�
% ������˵���ڼ���ͼ������˹����Ĺ�������Ҫ�õ�ģ�͵Ľṹ��Ϣ��Ҳ������Ҫ�жϲ�ͬ�Ľڵ��Ƿ�����ͬһ�����٣���һ����Ҫ��Ͼ����ģ���������ض��Ĵ���
% ������ͨ���ڵ�����ı߽�����жϽڵ���������٣����Դ�comsol����л�ȡ���������жϳ��ڵ��Ƿ�����ͬһ��������֯��
function matrix_L = computationL(nums, nodesMat)
if exist('Step_scadM\L.mat','file')    % �ж��ļ��Ƿ����
    disp('��ȡL����');
    load('Step_scadM\L.mat');
else
    % ����L����-��̬ͼ������˹����Manifold Learning��
    BoundaryIndex = nodesMat(:, 4);   % ��ȡ���ڵ�����ı߽���� n*1
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
    L = zeros(nums, nums);   % ��ʼ��L����Ϊn*n��0����
    for i = 1 : nums
        for j = 1 : nums
            if i == j
                L(i, j) = 1;   
            else
                % �ж������ڵ��Ƿ�ͬ����һ������
                if  isSameOrgan(i, j, index_heart) || isSameOrgan(i, j, index_liver) || isSameOrgan(i, j, index_bone) || isSameOrgan(i, j, index_lungs1) || isSameOrgan(i, j, index_lungs2)          % i, j �ڵ㶼����ͬһ������Sk
                    dij = sqrt(sum((nodes(i, :) - nodes(j, :)).^2));   % ����ڵ�i��ڵ�j֮���ŷ����þ���
                    if isSameOrgan(i, j, index_heart)    % �ڵ�ijͬ����������֯�����
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
    parSave('Step_scadM\L.mat','L', L);   % ����L���󣬱����ظ�����
end
matrix_L = L;
end

% step4: ʹ��ţ�ٷ���������x
function [x]=newton(x0,max_iterations, G, b_analystic, nodesMat, nums, L) 
xk=x0;
k=1;	        % ��������ܴ��� 
while max_iterations >= k
    minf = fun1(G, xk, b_analystic, nums, L);
    minmu = 0.01;
    dfx=fun2(G, xk, b_analystic, nums, L);
    xk=xk-minmu*dfx;  % mu�̶�Ϊһ����0.01
    xk = resetZero(xk, nodesMat);
    k=k+1;
end
parSave('Step_scadM\output.mat','output', xk);   % ����output����
x = xk;
end 