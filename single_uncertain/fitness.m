
% 计算适应度

function [fitness_value, pro] = fitness (pop_repaired)
    
    global N_q  cost_chidao  cost_work  cost_diff  M_pre  ED  best_position  cost_delay    N_v   W_u  W_d  cost_stay  sample  g   L_m L % popsize
    
    % cost_zaodao = 0.9;    % 船舶提前到港单位时间成本
    
    % 适应度函数为目标函数倒数，目标函数包括码头间转运成本、延误离港惩罚成本、作业成本、泊位偏移成本
    
    % 判断指标：
    % 1、目标函数值；
    % 2、船舶的即到即靠率；-----实际到港时间与实际靠泊时间相同的船舶所占比例
    % 3、船舶的延误离港率；-----实际离港时间超过期望离港时间的船舶所占比例
    % 4、船舶延迟到港率；-------实际到港时间超过调度靠泊时间的船舶所占比例
    % 5、岸线利用率；
    % 6、岸桥利用率；
    
    % 对每个染色体进行适应度计算
    object_value = zeros(size(pop_repaired, 1), size(sample{1, 1}, 1));  % 各染色体在各样本下的目标函数值集合初始化
    object_value_mean = zeros(1, size(pop_repaired, 1));    % 各染色体在所有样本下所得目标函数值的期望集合初始化
    object_value_std = zeros(1, size(pop_repaired, 1));    % 各染色体在所有样本下所得目标函数值的方差集合初始化
    is_change_terminal = zeros(1, N_v);    % 判断船舶的实际停靠码头是否为其预分配码头
    fitness_value = zeros(1, size(pop_repaired, 1));    % 各染色体的适应度值
    pro = zeros(8, size(pop_repaired, 1));
    % pro(1, :) = pro_same = zeros(1, popsize);    % 船舶即到即靠率
    % pro(2, :) = pro_delay_leave = zeros(1, popsize);    % 船舶延误离港率
    % pro(3 : ) = pro_delay_arrive = zeros(1, popsize);    % 船舶延误靠泊率
    % pro(4: 6, : ) = pro_quay_use = zeros(popsize, 3);    % 各码头岸线利用率
    % pro(7, :) = pro_all_quay_use = zeros(1, popsize);     % 总的岸线利用率
    % pro(8, :) = pro_crane_use = zeros(1, popsize);    % 总的岸桥利用率
    
    
    % 在计算岸线利用率时还需要将计划期开始时已靠泊船舶计算进去
    % 首先得到各码头已靠泊船舶数据（船身长度、离港时间两项）
    vessel_already = find_vessel_already();
    
    for i = 1 : size(pop_repaired, 1)
        
        num_same = 0;    % 实际到港时间与实际靠泊时间相同的船舶数量
        num_delay_leave = 0;    % 实际离港时间超过期望离港时间的船舶数量
        num_delay_arrive = 0;    % 实际到港时间超过调度靠泊时间的船舶数量
        num_early_arrive = 0;    % 实际到港时间早于调度靠泊时间的船舶数量
        num_quay_use = zeros(1, 3);    % 各码头岸线总利用长度时间 N_m
        num_crane_use = zeros(1, 3);    % 各码头岸桥总利用作业时
        
        % 计算已靠泊船舶所占用的岸线资源长度时间
        for j = 1 : 3
            for m = 1 : size(vessel_already{1, j}, 2)
                num_quay_use(1, j) = num_quay_use(1, j) + vessel_already{1, j}(3, m) * vessel_already{1, j}(5, m);
                num_crane_use(1, j) = num_crane_use(1, j) + vessel_already{1, j}(3, m) * vessel_already{1, j}(6, m);
            end
        end
        
        % 先确定待调度船舶的最终停靠码头是否为其预分配码头，得到矩阵is_change_terminal
        for j = 1 : N_v
            if pop_repaired(i, j) == M_pre(j)
                is_change_terminal(j) = 1; % 停靠码头即为预分配码头
            else
                is_change_terminal(j) = 0; % 停靠码头不为其预分配码头
            end
        end
        % 至此，得到船舶是否更改其码头信息
        
        for k = 1 : size(sample{1, 1}, 1)
            for l = 1 : size(sample{1, 1}, 2)
                % 实际离港时间的计算
                % 需要先判断调度靠泊时间与实际到港时间之间的大小，选择其中大的作为实际靠泊时间
                if pop_repaired(i, l + 2*N_v) < sample{1, 1}(k, l)
                    num_delay_arrive = num_delay_arrive + 1;    % 延迟到达的船舶将其靠泊时刻记为其到港时刻
                    sample{1, 3}(k, l) = sample{1, 1}(k, l) + round((W_u(l) + W_d(l)) / (sample{1, 2}(k, l) * pop_repaired(i, l + 3*N_v) * g ^ (pop_repaired(i, l + 3*N_v) - 1)) , 1);
                elseif pop_repaired(i, l + 2*N_v) == sample{1, 1}(k, l)
                    num_same = num_same + 1;    % 即到即靠船舶数量
                    sample{1, 3}(k, l) = sample{1, 1}(k, l) + round((W_u(l) + W_d(l)) / (sample{1, 2}(k, l) * pop_repaired(i, l + 3*N_v) * g ^ (pop_repaired(i, l + 3*N_v) - 1)) , 1);
                elseif pop_repaired(i, l + 2*N_v) > sample{1, 1}(k, l)
                    num_early_arrive = num_early_arrive + 1;     % 提前到港的船舶数量
                    sample{1, 3}(k, l) = pop_repaired(i, l + 2*N_v) + round((W_u(l) + W_d(l)) / (sample{1, 2}(k, l) * pop_repaired(i, l + 3*N_v) * g ^ (pop_repaired(i, l + 3*N_v) - 1)) , 1);
                end
                
                if sample{1, 3}(k, l) > ED(l)
                    num_delay_leave = num_delay_leave + 1;    % 延迟离港的船舶数量
                end
                
                % 进行岸线利用量的计算
                % 首先，判断各船舶的靠泊码头，根据码头计算总的岸线利用量    判断靠泊时间是否大于计划期长度，若大于，则不进行该船舶计算，若不大于，则判断离港时间是否大于计划期长度，若大于，则用计划期长度去计算在港时间，若不大于，则用实际离港时间计算在港时间
                if max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)) < 72
                    if sample{1, 3}(k, l) > 72
                        if pop_repaired(i, l) == 1
                            num_quay_use(1, 1) = num_quay_use(1, 1) + L(l) * (72 - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));
                            num_crane_use(1, 1) = num_crane_use(1, 1) + pop_repaired(i, l + 3*N_v) * (72 - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));
                        elseif pop_repaired(i, l) == 2
                            num_quay_use(1, 2) = num_quay_use(1, 2) + L(l) * (72 - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));
                            num_crane_use(1, 2) = num_crane_use(1, 2) + pop_repaired(i, l + 3*N_v) * (72 - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));
                        elseif pop_repaired(i, l) == 3
                            num_quay_use(1, 3) = num_quay_use(1, 3) + L(l) * (72 - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));
                            num_crane_use(1, 3) = num_crane_use(1, 3) + pop_repaired(i, l + 3*N_v) * (72 - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));
                        end
                    else
                        if pop_repaired(i, l) == 1
                            num_quay_use(1, 1) = num_quay_use(1, 1) + L(l) * (sample{1, 3}(k, l) - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));
                            num_crane_use(1, 1) = num_crane_use(1, 1) + pop_repaired(i, l + 3*N_v) * (sample{1, 3}(k, l) - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));
                        elseif pop_repaired(i, l) == 2
                            num_quay_use(1, 2) = num_quay_use(1, 2) + L(l) * (sample{1, 3}(k, l) - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));
                            num_crane_use(1, 2) = num_crane_use(1, 2) + pop_repaired(i, l + 3*N_v) * (sample{1, 3}(k, l) - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));
                        elseif pop_repaired(i, l) == 3
                            num_quay_use(1, 3) = num_quay_use(1, 3) + L(l) * (sample{1, 3}(k, l) - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));
                            num_crane_use(1, 3) = num_crane_use(1, 3) + pop_repaired(i, l + 3*N_v) * (sample{1, 3}(k, l) - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));
                        end
                    end
                end
                % 下一步，进行适应度值的计算                
                % 问题：在样本下，当船舶的实际到港时间超过其调度靠泊时间时，就会使该船的实际靠泊时间延迟，从而可能影响到后续船舶的靠泊，因此，在进行目标函数值计算时的指标可能是不正确的
                % 解决方案？？：首先判断在样本下，该船的调度方案是否可行，若可行，则进行成本计算，若不可行，则添加一个惩罚成本。
                
                % 目标函数值的计算：
                % 1、作业成本 = （实际离港时间 - 实际靠泊时间）* 岸桥数量 * 单位作业成本；
                % 2、转运成本 = 转运量 * 单位转运成本
                % 3、位置偏离成本 = 装卸量 * 偏移距离 * 单位偏离成本
                % 4、提前到港成本 = 提前到港时间 * 单位等待成本
                % 5、延误到港成本 = 延误到港时间 * 单位惩罚成本
                % 6、延误离港成本 = 延误离港时间 * 单位延误成本
                
                object_value(i, k) = object_value(i, k) + cost_work * pop_repaired(i, l + 3*N_v) * (sample{1, 3}(k, l) - max(pop_repaired(i, l + 2*N_v), sample{1, 1}(k, l)));    % 作业成本
%                 object_value(i, k) = object_value(i, k) + cost_tr(M_pre(l), pop_repaired(i, l)) * W_u(l);    % 转运成本
                object_value(i, k) = object_value(i, k) + cost_diff * (W_u(l) + W_d(l)) * abs(pop_repaired(i, l + N_v) - best_position(l)) * is_change_terminal(l);    % 位置偏离成本
                object_value(i, k) = object_value(i, k) + cost_stay * max(pop_repaired(i, l + 2 * N_v) - sample{1, 1}(k, l), 0);    % 提前到港成本
                object_value(i, k) = object_value(i, k) + cost_chidao * max(sample{1, 1}(k, l) - pop_repaired(i, l +2*N_v), 0);    % 延误到港成本
                object_value(i, k) = object_value(i, k) + cost_delay(l) * max(sample{1, 3}(k, l) - ED(l), 0);    % 延误离港成本

            end
        end        
        % 至此，得到所有在染色体 i 所表示的解下可行样本的目标函数值
        
        % 计算期望及方差，得到适应度函数值
        object_value_mean(i) = mean(object_value(i, :));     % 所有可行样本函数值的期望
        object_value_std(i) = std(object_value(i, :));     % 所有可行样本函数值的标准差
        fitness_value(1, i) = 100000 / (object_value_mean(i) + object_value_std(i));    % 各染色体的适应度值
        % 至此，得到染色体 i 的适应度值
        
        % 计算各染色体的评价指标：
        pro(1, i) = num_same / (size(sample{1, 1}, 1) * size(sample{1, 1}, 2));          % 1、船舶即到即靠率 = 即到即靠船舶数量 /（样本数量 × 船舶数量）
        pro(2, i) = num_delay_leave / (size(sample{1, 1}, 1) * size(sample{1, 1}, 2));       % 2、船舶离港延误率 = 延误离港船舶数量 /（样本数量 × 船舶数量）
        pro(3, i) = num_delay_arrive / (size(sample{1, 1}, 1) * size(sample{1, 1}, 2));       % 3、船舶靠泊延误率 = 靠泊延误船舶数量 /（样本数量 × 船舶数量）
        for j = 1 : 3
            pro(j+3, i) = num_quay_use(1, j) / (L_m(j) * 72 * size(sample{1, 1}, 1));      % 4、码头 i 岸线利用率 = （停靠在码头 i 的船舶长度 × 船舶靠泊时间）/（计划期长度 × 岸线长度）
        end
        % 至此，得到染色体 i 的评价指标值
        % 加上所有码头整个岸线的利用率
        pro(7, i) = sum(num_quay_use) / (size(sample{1, 1}, 1) * 72 * sum(L_m));
        
        % 加上所有岸桥的利用率
        pro(8, i) = sum(num_crane_use) / (size(sample{1, 1}, 1) * 72 * sum(N_q));
        
        % 至此，得到染色体 i 的评价指标值
        
    end
    % 至此，得到所有染色体的适应度值
end
