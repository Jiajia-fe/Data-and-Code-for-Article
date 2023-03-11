function pop_init = population_init ()
    
    global popsize N_v A C_min C_max W_u W_d v_m g A_sigma v_sigma  M_pre best_position  N_q L_m L
    
    % 首先生成 popsize × 5N_v 的零矩阵，其中 1 ~ N 为靠泊码头；
    %                                      N+1 ~ 2N 为靠泊位置；
    %                                      2N+1 ~ 3N 为靠泊时间；
    %                                      3N+1 ~ 4N 为岸桥数量；
    %                                      4N+1 ~ 5N 为岸桥起始编号；
    pop_init = zeros(popsize, 8*N_v);
    
    % 初始化种群的一半染色体采用贪婪构造的方法生成
    for i = 1 : popsize/2
       for j = 1 : N_v
            % 初始解对于算法效率来说有着很重要的影响，因此，在初始化种群时，采用贪婪构造和随机生成的方式，来生成一个质量较高的解，和随机解
            
            % 第一部分：靠泊码头，初始化时将其取为其预分配码头
            % pop_init(i, j) = randi([1, 3], 1, 1);
            pop_init(i, j) = M_pre(j);
            
            % 第二部分：靠泊位置，因为前面已将船舶分配给其预分配码头
            % 因此在靠泊位置上这里初始化时直接初始化为其最优靠泊位置
            % pop_init(i, j+N_v) = randi([0, L_m(pop_init(i, j)) - L(j)], 1, 1);
            pop_init(i, j+N_v) = best_position(j);
            
            % 第七部分：到港时间的松弛参数
            pop_init(i, j+6*N_v) = A_sigma(j);
            
            % 第八部分：岸桥作业效率松弛参数
            pop_init(i, j+7*N_v) = v_sigma(j);
            
            % 第三部分：靠泊时刻
            % 这里靠泊时刻将其初始化为其到达时刻
            pop_init(i, j+2*N_v) = A(j) + pop_init(i, j+6*N_v);
            
            % 第四部分：岸桥数量
            pop_init(i, j+3*N_v) = randi([C_min(j), C_max(j)], 1, 1);
            
            % 第五部分：起始岸桥编号
            % pop_init(i, j+4*N_v) = randi([1, N_q(pop_init(i, j)) - pop_init(i, j+3*N_v) + 1], 1, 1);
            pop_init(i, j+4*N_v) = 1;
            
            % 第六部分：离港时刻计算
            pop_init(i, j+5*N_v) = round(pop_init(i, j+2*N_v) +  (W_u(j) + W_d(j)) / (pop_init(i, j+3*N_v) * (v_m(pop_init(i, j)) - pop_init(i, j+7*N_v)) * g^(pop_init(i, j+3*N_v) - 1)), 1);
            
       end
    end
    
    % 初始种群中的一半染色体采用随机生成的方法
    for i = popsize/2+1 : popsize
       for j = 1 : N_v
            
            % 初始解对于算法效率来说有着很重要的影响，因此，在初始化种群时，采用贪婪构造和随机生成的方式，来生成一个质量较高的解，和随机解
            
            % 第一部分：靠泊码头，初始化时将其取为其预分配码头
            pop_init(i, j) = randi([1, 3], 1, 1);
            % pop_init(i, j) = M_pre(j);
            
            
            % 第二部分：靠泊位置，因为前面已将船舶分配给其预分配码头
            % 因此在靠泊位置上这里初始化时直接初始化为其最优靠泊位置
            pop_init(i, j+N_v) = randi([0, L_m(pop_init(i, j)) - L(j)], 1, 1);
            % pop_init(i, j+N_v) = best_position(j);
            
            % 第七部分：到港时间的松弛参数
            pop_init(i, j+6*N_v) = A_sigma(j);
            
            % 第八部分：岸桥作业效率松弛参数
            pop_init(i, j+7*N_v) = v_sigma(j);
            
            % 第三部分：靠泊时刻
            % 这里靠泊时刻将其初始化为其到达时刻
            pop_init(i, j+2*N_v) = A(j) + pop_init(i, j+6*N_v);
            
            % 第四部分：岸桥数量
            pop_init(i, j+3*N_v) = randi([C_min(j), C_max(j)], 1, 1);
            
            % 第五部分：起始岸桥编号
            pop_init(i, j+4*N_v) = randi([1, N_q(pop_init(i, j)) - pop_init(i, j+3*N_v) + 1], 1, 1);
            % pop_init(i, j+4*N_v) = 1;
            
            % 第六部分：离港时刻计算
            pop_init(i, j+5*N_v) = round(pop_init(i, j+2*N_v) + (W_u(j) + W_d(j)) / (pop_init(i, j+3*N_v) * (v_m(pop_init(i, j)) - pop_init(i, j+7*N_v)) * g^(pop_init(i, j+3*N_v) - 1)), 1);
            
       end
    end
    
end
