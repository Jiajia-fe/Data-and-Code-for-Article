
% 单码头不确定性模型求解代码
% 修改思路：
% 1、初始化时，直接将船舶的靠泊码头设置为其预分配码头，且在后续的修复、交叉、变异中不变化

% 算法改进方向：
% 1、自适应：随着迭代次数的增加而不断改变交叉概率以及变异概率
% 2、全局搜索与局部搜索相结合：在每一次迭代中，将遗传算法搜索出来的结果在进行局部搜索（模拟退火、禁忌搜索。。。）
% 3、多种群并行搜索：生成两个种群同时进行搜索

% 还存在的疑惑：
% 1、不确定性在论文中还是体现的不够明显，并没有得到具体的解决
% 2、算法性能不知道怎么进行比较验证，因为数据集是自己随机生成的

% 后续工作：
% 1、指标的计算还要添加。目前只能得到一个方案的示意图，但方案的相关指标还没有得到（已完成）
% 2、然后将计划到港船舶数据集扩充，目前只有20艘船，得到的最优方案有很多，它们的目标函数值都一样（数据已修改，分别为20、30、40艘船到港时建立了10个数据集）
% 3、按照算法改进方向，将算法在进行修改，使其能有所不同（交叉改为了熔合算子，可能有点老了；交叉用了自适应的，不过可能还需要修改；后面又加了个邻域搜索的，不过还想太简单了）
% 4、文字部分进行修改，指标进行画图、算法步骤进行说明等，这部分完成后可以先试着投出去（前面部分文字已修改，后续还要对结果描述）
% 5、继续研究不确定性的融入（看相关论文），看能不能找到更合适的建模方法，若能找到，在进行论文的修改（待完成）
% 6、大论文方向：继续将问题的约束等问题扩展，添加岸桥服务范围等新的实际约束，来提高论文的充实度，其他方面的约束还要看论文（待完成）
% 7、后面进行大论文的书写（待完成）
% 8、可以试下使用GPU、并行计算，看能不能加快求解速度，目前代码感觉有点臃肿，堆砌严重，同时适应度计算时因为要跑好几个样本，速度比较慢

clear;
clc;

global cost_chidao sample g cost_work cost_diff N_m L_m d_m N_q v_m cost_tr N_0 M_0 W_0 L_0 x_0 C_0 Q_first_0 d_0 N_v L M_pre W_u W_d A A_sigma v_sigma ED C_min C_max draft best_position cost_delay popsize G_max  cost_stay;
%基础数据
jiajia_1=readmatrix('D:\lunwen\ALL\paper_1\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'B3:J3');
g=jiajia_1(1);
cost_work =jiajia_1(2);
cost_diff=jiajia_1(3);
cost_stay=jiajia_1(4);
cost_chidao=jiajia_1(9);
%{
% 其他数据 %
% 岸桥干扰系数 %
g = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'B3:B3');

% 单位时间单位岸桥作业成本 %   这些成本常量可以不用作为全局变量，只需放到 fitness 函数中即可
cost_work = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'C3:C3');

% 单位时间非生产性停泊费用等待成本
cost_stay = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'E3:E3');

% 单位集装箱单位距离偏离成本 %
cost_diff = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'D3:D3');
%迟到成本
cost_chidao = 100;  % readmatrix('E:/Paper/me/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'J3:J3');
%}
%码头数据
jiajia_2=readmatrix('D:\lunwen\ALL\paper_1\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'A6:F8');
N_m=jiajia_2(1,1);
L_m=jiajia_2(:,3);
d_m=jiajia_2(:,4);
N_q=jiajia_2(:,5);
v_m =jiajia_2(:,6);
%{
% 码头数据 %
% 码头数量 %
N_m = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'A6:A6');

% 岸线长度 %
L_m = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'C6:C8')';

% 水深条件 %
d_m = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'D6:D8')';

% 岸桥数量 %
N_q = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'E6:E8')';

% 岸桥作业效率
v_m = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'F6:F8')';
%}
% 码头间单位箱量转运成本 %
cost_tr = 2 .* readmatrix('D:\lunwen\ALL\paper_1\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'C11:E13');
%已停靠船舶
jiajia_3=readmatrix('D:\lunwen\ALL\paper_1\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'A17:H24');
N_0 =jiajia_3(1,1);
M_0 =jiajia_3(:,3);
W_0 =jiajia_3(:,4);
L_0 =jiajia_3(:,5);
x_0 =jiajia_3(:,6);
C_0=jiajia_3(:,7);
Q_first_0=jiajia_3(:,8);
%{
% 已靠泊船舶 %
% 已靠泊船舶总数量 %
N_0 = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'A17:A17');

% 所靠泊码头 %
M_0 = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'C17:C24')';

% 剩余作业量 %
W_0 = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'D17:D24')';

% 船舶长度 %
L_0 = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'E17:E24')';

% 靠泊位置 %
x_0 = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'F17:F24')';

% 岸桥数量 %
C_0 = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'G17:G24')';

% 起始岸桥编号 %
Q_first_0 = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'H17:H24')';
%}
% 已靠泊船舶离港时刻
for j = 1:N_0
    d_0(j) = round(W_0(j) / (C_0(j) * v_m(M_0(j)) * g^(C_0(j)-1)), 1);
end
jiajia_4=readmatrix('D:\lunwen\ALL\paper_1\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'A28:O47');
N_v=jiajia_4(1,1);
M_pre=jiajia_4(:,3);
W_u =jiajia_4(:,4);
W_d=jiajia_4(:,5);
L=jiajia_4(:,6);
A=jiajia_4(:,7);
ED=jiajia_4(:,8);
C_min=jiajia_4(:,9);
C_max=jiajia_4(:,10);
draft=jiajia_4(:,11);
best_position=jiajia_4(:,12);
cost_delay=jiajia_4(:,13);
A_sigma=jiajia_4(:,14);
v_sigma=jiajia_4(:,15);
%{
% 待调度船舶 %
% 待调度船舶数量 %
N_v = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'A28:A28');

% 船舶长度 %
L = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'F28:F47')';

% 预分配码头 %
M_pre = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'C28:C47')';

% 待装箱量 %
W_u = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'D28:D47')';

% 待卸箱量 %
W_d = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'E28:E47')';

% 船舶预计到港时刻 % 期望值
A = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'G28:G47')';

% 船舶预计到港时刻标准差
A_sigma = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'N28:N47')';

% 岸桥作业效率标准差
v_sigma = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'O28:O47')';

% 期望离港时刻 %
ED = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'H28:H47')';

% 最小岸桥数量 %
C_min = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'I28:I47')';

% 最大岸桥数量 %
C_max = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'J28:J47')';

% 吃水深度 %
draft = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'K28:K47')';

% 偏好位置 %
best_position = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'L28:L47')';

% 单位时间离港延迟惩罚 %
cost_delay = 2 .* readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'M28:M47')';
%}
% 样本数据 %  同码头成本系数，这些样本数据也可以直接放到 fitness 函数中
% 实际到达时间
sample{1, 1} = readmatrix('D:\lunwen\ALL\paper_1\input/data_20.xlsx', 'Sheet', 'sample_9_1', 'Range', 'B4:U23');
% 岸桥实际作业效率
sample{1, 2} = readmatrix('D:\lunwen\ALL\paper_1\input/data_20.xlsx', 'Sheet', 'sample_9_2', 'Range', 'B4:U23');
%算法参数
popsize=jiajia_1(5);
G_max=jiajia_1(6);
p_0=jiajia_1(8);
%{
% 算法参数 %
% 种群规模 %
popsize = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'F3:F3');

% 最大迭代次数 %
G_max = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'G3:G3');

% 交叉概率 %
% p_c = readmatrix('E:/Paper/me/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'H3:H3');
% p_1 = 0.8;    % 交叉概率1
% p_2 = 0.9;    % 交叉概率2
% p_3 = 0.2;    % 变异概率1
%p_4 = 0.3;    % 变异概率2
% 变异概率 %
p_0 = readmatrix('D:\lunwen\ALL\paper_2\input/data_20.xlsx', 'Sheet', 'data_9', 'Range', 'I3:I3');
%}
% 种群初始化 %   
pop_init = population_init();

% 修复不可行解 %
pop_repaired = repair(pop_init);

% 计算适应度（分别针对各个目标函数）
[fitness_value, pro_all] = fitness(pop_repaired);          % #######################  至此  ####################

% 按从小到大的顺序进行排序
Matrix_1 = sortrows([pop_repaired fitness_value' pro_all'], size(pop_repaired, 2) + 1);
pop_repaired = Matrix_1( : , 1:size(pop_repaired, 2));
fitness_value = Matrix_1(:, size(pop_repaired, 2) + 1)';
pro_all = Matrix_1(:, size(pop_repaired, 2) + 2 : size(Matrix_1, 2))';

% 对合并后的矩阵进行非支配排序并计算拥挤度
% pop_non_domin_sort = non_domination_sort_mod(pop_repaired, fitness);

g_unchange = 1;
fitness_value_max = zeros(1, G_max);
fitness_value_mean = zeros(1, G_max);
best_chromosome = zeros(G_max, 8*N_v);
pro = zeros(8, G_max);

for i = 1 : G_max
    
    % 这里的pop_repaired和fitness_value都已经按从小到大的顺序排列过
%     pop_old = pop_repaired;
    pop_temp = [pop_repaired fitness_value' pro_all'];
    % 第一步：找出适应度最高的几个染色体（n），将他们直接放进下一代种群中（与新种群中适应度最差的五个进行比较后决定是否进行）
    best_five = pop_temp(0.9*popsize + 1 : popsize, :);
    
    % 至此，得到适应度最好的10%个染色体
    
    
    % 第二步：根据适应度值计算选择概率，确定各染色体的选择概率
    proba_choose = calcu_proba(fitness_value);
    % 至此，得到各染色体被选中的概率
    
    % 第三步：按照轮盘赌选择出种群规模个（或popsize-n）个体
    % 按照竞标赛选择
    pop_select = select(pop_repaired, proba_choose);
    
    % -------------------------------------------- %
    [fitness_value, ~] = fitness(pop_select);    % 适应度值矩阵
%     fitness_value_max_temp = max(fitness_value);    % 最大适应度值
%     fitness_value_mean_temp = mean(fitness_value);    % 平均适应度值
    % -------------------------------------------- %
    
    % 至此，选择步骤已完成，下一步进行交叉
    pop_cross = pop_select;
    pop_muta = zeros(popsize, 8*N_v);
    
    % 后续这里可以只产生0.9*popsize个个体，然后和前面的0.1*popsize个个体进行合并形成新种群
    for j = 1 : popsize
        
        % 第四步：选出用于交叉的染色体，然后将交叉产生的新染色体放入pop_select中（两点交叉）
        
        % 随机选取两个染色体
        r_in = randperm(popsize, 2);
        
        % 首先判断是否需要进行交叉操作，使用熔合算子进行交叉
        
        % 对交叉和变异概率进行确定    ---------------------------- %
        % Pc = P1 * (fitness_value_max - fitness_i) / (fitness_value_max - fitness_value_mean)
        % Pc = P2;
        
        f_1 = fitness_value(r_in(1));
        f_2 = fitness_value(r_in(2));
        p_c = f_1 / (f_1 + f_2);
        for k =1 : size(pop_select, 2)
            
            r = rand();
            if r <= p_c
                pop_cross(j, k) = pop_select(r_in(1), k);
            else
                pop_cross(j, k) = pop_select(r_in(2), k);
            end
        end
        % 至此，交叉操作已完成，下一步进行变异操作
        
        % 变异和交叉不同时进行？？？？
        % 交叉后的种群在重新进行修复、适应度值计算，然后再根据适应度值确定变异概率，根据变异概率进行染色体的变异操作
        
        % 首先，进行变异操作，采用均匀变异
        % 先确定变异概率
        p_m = p_0 * sqrt(g_unchange / (0.1 * G_max));
        % 然后根据变异概率进行变异操作
        pop_muta(j, :) = mutation(pop_cross(j, :), p_m);
%         pop_muta(j+1, :) = mutation(pop_cross(j+1, :));
        % 至此，变异操作已完成，进入下一个循环
    end
    % 至此，得到经过交叉变异后产生的新种群，下一步需要对产生的新种群进行基因修复
       
    
    
    
    % 第六步：对交叉、变异后的种群进行基因修复，使得解可行
    pop_repaired = repair(pop_muta);
    % 至此，完成基因的修复作业，下一步为进行适应度的计算
    
    % 第七步：进行适应度值的计算
    [fitness_value, pro_all] = fitness(pop_repaired);
    % 至此，已得到新种群的适应度，下一步，找出适应度最小的10%个个体，将其与父代中最好的10%个个体比较，择其优者10%个
    
    % 第八步：找出适应度最小的10%个个体，将其替换为父代中最好的10%个个体
    % 首先得到根据适应度值从小到大进行排序后的种群
    Matrix_2 = sortrows([pop_repaired fitness_value' pro_all'], size(pop_repaired, 2) + 1);
    
    % 然后将前10%个替换为父代最好的10%个个体
    Matrix_2(1 : 0.1*popsize, : ) = best_five;
    
    % 然后对pop_repaired重新按适应度值排序
    Matrix_2 = sortrows(Matrix_2, size(pop_repaired, 2) + 1);
    pop_repaired = Matrix_2(:, 1 : size(pop_repaired, 2));
    fitness_value = Matrix_2(:, size(pop_repaired, 2) + 1)';
    pro_all = Matrix_2(:, size(pop_repaired, 2) + 2 : size(Matrix_2, 2))';
    
    % 第九步：进行指标的计算、记录
    % 初步指标：每一代中的最大适应度值、最好的染色体、平均适应度值
    fitness_value_max(i) = max(fitness_value); % 每一次迭代中的最大适应度值
    fitness_value_mean(i) = mean(fitness_value); % 每一次迭代中的平均适应度值
    best_chromosome(i, :) = pop_repaired(size(pop_repaired, 1), :); % 每一次迭代中的最佳染色体
    pro(:, i) = pro_all(:, size(pro_all, 2));
    
    % 确定最大适应度值保持代数
    if i >=2
        if fitness_value_max(i) == fitness_value_max(i - 1)
            g_unchange = g_unchange + 1;
        else
            g_unchange = 1;
        end
    end
    
    
    % 这里对交叉以及变异概率进行确定。我们希望随着迭代次数的增加以及最优解不变代数的增加而使得交叉概率以及变异概率也随之增加
    
    % 判断最优解是否已经连续 XXX 代没有变化，若不是，则进入下一次迭代，若是，则进行邻域搜索
    % 判断迭代次数是否已超过一半，若是，则进行下一步，若不是，则进入下一代
%     if g_unchange >= 0.1 * G_max     % 判断最优解是否已经 10 代没有变化，若是，则进行下一步，若不是，则进入下一代
    t = 1;
    T_max = 100;
    T = T_max;
    while (g_unchange >= 0.1 * G_max && t < T_max)
        temp = best_chromosome(i, :);
        % pop_improved = improve(temp);    % 进行邻域搜索 --------定义邻域为当前调度方案中改变一艘船舶的靠泊信息时产生的解
        pop_improved = temp;
        r_vessel = randi([1, N_v], 1, 1);     % 随机选择一艘船舶
%         pop_improved(r_vessel) = randi([1, 3], 1, 1);     % 修改码头信息
        pop_improved(r_vessel + N_v) = randi([0, L_m(pop_improved(r_vessel)) - L(r_vessel)], 1, 1);    % 修改靠泊位置信息
        pop_improved(r_vessel + 6*N_v) = randi([-A_sigma(r_vessel) * 20, A_sigma(r_vessel) * 20 ], 1, 1)/10;     % 修改到港时间松弛
        pop_improved(r_vessel + 7*N_v) = randi([-20 * v_sigma(r_vessel), 20 * v_sigma(r_vessel) ], 1, 1)/10;
        pop_improved(r_vessel + 2*N_v) = randi([max(A(r_vessel) * 10 + pop_improved(r_vessel + 6*N_v)*10, 0), round(ED(r_vessel)*10, 0)], 1, 1)/10;    % 修改到港时间
        pop_improved(r_vessel + 3*N_v) = randi([C_min(r_vessel), C_max(r_vessel)], 1, 1); % 修改岸桥数量
        pop_improved(r_vessel + 4*N_v) = randi([1, N_q(pop_improved(r_vessel)) - pop_improved(r_vessel + 3*N_v) + 1], 1, 1); % 修改岸桥起始编号
        pop_improved(r_vessel + 5*N_v) = round(pop_improved(r_vessel + 2*N_v) + (W_u(r_vessel) + W_d(r_vessel))/(pop_improved(r_vessel+3*N_v) * (v_m(pop_improved(r_vessel)) - pop_improved(r_vessel+7*N_v)) * g^(pop_improved(r_vessel+3*N_v) - 1)), 1); % 修改离港时刻
        pop_improved = repair(pop_improved);    % 修复
        [fitness_value_new, pro_i] = fitness(pop_improved);    % 计算搜索得到的解的适应度值    pro的值要加上
        if fitness_value_new(1, 1) > fitness_value_max(i)
            % 当搜索得到的解质量更好时，将最优解替换掉
            pop_repaired(size(pop_repaired, 1), : ) = pop_improved;
            fitness_value(1, size(fitness_value, 2)) = fitness_value_new(1, 1);
            fitness_value_max(i) = fitness_value_new(1, 1);
            fitness_value_mean(i, :) = mean(fitness_value);
            best_chromosome(i, :) = pop_improved;
            pro_all(:, size(pro_all, 2)) = pro_i;
            pro(:, i) = pro_i;
            g_unchange = 1;
            t = t + 1;
        else
            % 差的解以一定概率被接受
            p = exp(-(fitness_value_max - fitness_value_new(1, 1))/T);
            r_accept = rand();
            if r_accept < p
                pop_repaired(1, : ) = pop_improved;
                fitness_value(1, 1) = fitness_value_new(1, 1);
                pro_all(:, 1) = pro_i;
                fitness_value_mean(i, :) = mean(fitness_value);
%                 disp(size(pop_repaired));
%                 disp(size(fitness_value));
                Matrix_3 = sortrows([pop_repaired fitness_value' pro_all'], size(pop_repaired, 2) + 1);
                pop_repaired = Matrix_3( : , 1 : size(pop_repaired, 2));
                fitness_value = Matrix_3( : , size(pop_repaired, 2) + 1)';
                pro_all = Matrix_3(:, size(pop_repaired, 2) + 2 : size(Matrix_3, 2))';
            end
            T = 0.8 * T;
            t = t + 1;
        end
    end
    % 进入下一代循环
end




% 画图部分
plan = zeros(8, N_v);
% 将最后一代中产生的最优染色体作为最终方案
for j = 1 : 8
    plan(j, :) = best_chromosome(G_max, (j-1)*N_v + 1 : j * N_v);
end
plan = [1 : N_v; plan];

% 最终方案中各码头的方案
plan_singel = cell(1, N_m);
%Y_0 = cell(1, N_m); % 各码头到港时间
for j = 1 : N_v
    for k = 1 : N_m
        if plan(2, j) == k
            plan_singel{1, k} = [plan_singel{1, k} plan(:, j)];
%            Y_0{1, k} =[Y_0{1, k} A(j) + plan(8, j)];
        end
    end
end

% 第十步：进行图像的呈现

figure('Name', 'fitness_value'); % 调出画图界面
clf; % 清除图像
plot(1: G_max, fitness_value_max, 'r');
% hold on;
% plot(1: G_max, fitness_value_mean, 'b');
% 至此，得到最优适应度函数值以及平均适应度函数值随代数的变化曲线

% 待调度船舶坐标
X_1 = cell(1, N_m);
Y_1 = cell(1, N_m);
Y_2 = cell(1, N_m);
% 已停靠船舶数据
X0_1 = cell(1, N_m);
X0_delta = cell(1, N_m);
Y0_1 = cell(1, N_m);
Y0_2 = cell(1, N_m);

% 对方案进行呈现 plot（X， Y）

vessel_already_0 = find_vessel_already ();

for k = 1 : N_m
    % 码头 k
    % 待调度船舶坐标数据
    % 起始横坐标集合
    X_1{1, k} = plan_singel{1, k}(3, :);
    % 起始纵坐标集合
    Y_1{1, k} = plan_singel{1, k}(4, :);
    % 终止纵坐标集合
    Y_2{1, k} = plan_singel{1, k}(7, :);
    
    % 已停靠船舶坐标数据
    % 起始横坐标
    X0_1{1, k} = vessel_already_0{1, k}(4, :);
    X0_delta{1, k} = vessel_already_0{1, k}(5, :);
    Y0_1{1, k} = vessel_already_0{1, k}(2, :);
    Y0_2{1, k} = vessel_already_0{1, k}(3, :);
    
    figure(k+1);
    clf;
    hold on;

    % 已靠泊船舶
    for i = 1 : size(vessel_already_0{1, k}, 2)
        plot([X0_1{1, k}(i), X0_1{1, k}(i), X0_1{1, k}(i) + X0_delta{1, k}(i), X0_1{1, k}(i) + X0_delta{1, k}(i), X0_1{1, k}(i)], [Y0_1{1, k}(i), Y0_2{1, k}(i), Y0_2{1, k}(i), Y0_1{1, k}(i), Y0_1{1, k}(i)], 'g');
    end
    
    % 待调度船舶
    for i = 1 : size(plan_singel{1, k}, 2)
        plot([X_1{1, k}(i), X_1{1, k}(i), X_1{1, k}(i) + L(plan_singel{1, k}(1, i)), X_1{1, k}(i) + L(plan_singel{1, k}(1, i)), X_1{1, k}(i)], [Y_1{1, k}(i), Y_2{1, k}(i), Y_2{1, k}(i), Y_1{1, k}(i), Y_1{1, k}(i)], 'r');
    end
    
    % 待调度船舶的实际总到港时间矩阵
%    for i = 1 : size(plan_singel{1, k}, 2)
%        plot([X_1{1, k}(i), X_1{1, k}(i), X_1{1, k}(i) + L(plan_singel{1, k}(1, i)), X_1{1, k}(i) + L(plan_singel{1, k}(1, i)), X_1{1, k}(i)], [Y_0{1, k}(i), Y_2{1, k}(i), Y_2{1, k}(i), Y_0{1, k}(i), Y_0{1, k}(i)], 'b--');
%    end
    
    axis([0 L_m(k) 0 150]);
end

disp('it is over!!!');






