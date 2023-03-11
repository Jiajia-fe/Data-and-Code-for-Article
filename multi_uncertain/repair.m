
% 修复函数

function pop_repaired = repair (pop_init)
    
    global  N_v d_m draft N_m L_m N_q L C_min W_u W_d v_m g N_0 C_max A M_pre best_position  % popsize
    
    pop_repaired = pop_init;  % 修复后的种群
    
    
    for i = 1 : size(pop_repaired, 1) 
        
        % 计划期开始时各码头已停靠船舶的信息
        vessel_already = find_vessel_already(); % 已完成调度船舶集合
        
        % 靠泊码头修复
        for j = 1 : N_v
            % 对分配到的码头的水深不符合其吃水要求的船舶进行码头信息修复
            if draft(j) > d_m(pop_repaired(i, j))
                if draft(j) <= d_m(2)
                    pop_repaired(i, j) = randi([2, 3], 1, 1);
                else
                    pop_repaired(i, j) = 3;
                end
            end
            
            % 将靠泊时间调整为其到达时间
            % 相应的离港时刻也进行更新
            pop_repaired(i, j + 2*N_v) = max(A(j) + pop_repaired(i ,j + 6*N_v), 0);
            pop_repaired(i, j+5*N_v) = round(pop_repaired(i, j+2*N_v) + (W_u(j) + W_d(j)) / (pop_repaired(i, j+3*N_v) * (v_m(pop_init(i, j)) - pop_init(i, j+7*N_v)) * g^(pop_repaired(i, j+3*N_v) - 1)), 1);
            
            % 当船舶实际靠泊码头为其预分配码头时，将其靠泊位置调整为其最优靠泊位置
            if pop_repaired(i, j) == M_pre(j)
                pop_repaired(i, j + N_v) = best_position(j);
            end
            
        end
        
        % 各码头待调度船舶集合
        vessel_ready = find_vessel_ready(pop_repaired(i, :)); % 码头信息修复后，各码头待调度船舶集合

        vessel_now = cell(1, N_m);
        idle_space = cell(1, N_m);
        useful_space = cell(1, N_m);

        for k = 1 : N_m
            % 对码头 k 待调度船舶靠泊位置、时间、岸桥数量、起始编号的修复
            while size(vessel_ready{1, k}, 2) > 0 % 只要还存在未完成调度的船舶
%                 disp(size(vessel_ready{1, k}, 2));
                
                % 在港船舶集合，第一行：船舶编号；第二行：靠泊时刻；第三行：离港时刻；第四行：靠泊位置；第五行：船身长度；第六行：分配岸桥数量；第七行：岸桥起始编号
                vessel_ready{1, k} = sortrows(vessel_ready{1, k}', 2)';  % 根据船舶的靠泊时间进行排序
                vessel_already{1, k} = sortrows(vessel_already{1, k}', 4)'; % 按靠泊位置以小到大进行排序
                vessel_now{1, k} = vessel_already{1, k}; % 调度时已在港船舶集合
                
                num_already = 1;
                for j = 1 : size(vessel_already{1, k}, 2) % 已完成调度船舶数量
%                     disp("---------------------");
%                     disp( size(vessel_already{1, k}, 2) );
                    if vessel_already{1, k}(3, j) <= vessel_ready{1, k}(2, 1) % 找出离港时刻早于当前调度船舶的所有船舶，并将其从已调度船舶集合中删除
                        vessel_now{1, k}(:, num_already) = [];
                    else
                        num_already = num_already + 1;
                    end
                end
                % 到此，已得到在靠泊时刻的所有在港靠泊船舶信息，即vessel_now{1, k}，下一步是得到空闲岸线区间
                
                % 得到可用岸线集合，判断当前是否有船舶，若有船舶，则先得到空闲岸线区间，若无船舶，则全段岸线均可用
                if size(vessel_now{1, k}, 2) > 1
                    
                    idle_space{1, k} = zeros(4, size(vessel_now{1, k}, 2) + 1);
                    
                    % 空闲岸线区间集合，第一行：空闲岸线区间起点坐标；第二行：空闲岸线区间终点坐标；第三行：区间内可用岸桥数量；第四行：区间内岸桥起始编号
                    
                    for j = 2 : size(vessel_now{1, k}, 2)  % 码头 k 的空闲岸桥区间数量为size()+1
                        idle_space{1, k}(1, j) = vessel_now{1, k}(4, j-1) + vessel_now{1, k}(5, j-1);  % 起点信息
                        idle_space{1, k}(2, j) = vessel_now{1, k}(4, j); % 区间末尾坐标
                        idle_space{1, k}(3, j) = vessel_now{1, k}(7, j) - vessel_now{1, k}(7, j-1) - vessel_now{1, k}(6, j-1); % 区间内可用岸桥数量
                        idle_space{1, k}(4, j) = vessel_now{1, k}(7, j-1) + vessel_now{1, k}(6, j-1); % 区间内岸桥起始编号
                    end
                    
                    % 第一段空闲岸线区间信息
                    idle_space{1, k}(1, 1) = 0; % 第一段空闲区间起始位置为0
                    idle_space{1, k}(2, 1) = vessel_now{1, k}(4, 1); % 第一段空闲区间的终止位置为最左边船舶的靠泊位置
                    idle_space{1, k}(3, 1) = vessel_now{1, k}(7, 1) - 1; % 第一段空闲区间内可用的岸桥数量
                    idle_space{1, k}(4, 1) = min(1, idle_space{1, k}(3, 1)); % 第一段空闲区间起始岸桥编号
                    
                    % 最后一段空闲区间岸线信息
                    idle_space{1, k}(1, size(vessel_now{1, k}, 2)+1) = vessel_now{1, k}(4, size(vessel_now{1 ,k}, 2)) + vessel_now{1, k}(5, size(vessel_now{1, k}, 2));  % 最后一段空闲区间起始位置
                    idle_space{1, k}(2, size(vessel_now{1, k}, 2)+1) = L_m(k); % 最后一段空闲区间终止位置
                    idle_space{1, k}(3, size(vessel_now{1, k}, 2)+1) = N_q(k) - vessel_now{1, k}(7, size(vessel_now{1, k}, 2)) - vessel_now{1, k}(6, size(vessel_now{1, k}, 2)) + 1; % 最后一段空闲区间内可用岸桥数量
                    idle_space{1, k}(4, size(vessel_now{1, k}, 2)+1) = vessel_now{1, k}(7, size(vessel_now{1, k}, 2)) + vessel_now{1, k}(6, size(vessel_now{1, k}, 2)); % 最后一段空闲区间起始岸桥编号
                    % 至此，已得到所有空闲区间，下一步为得到能够容纳正在调度船舶的船长的区间；
                    
                    % 可用岸线区间集合，第一行：空闲岸线区间起点坐标；第二行：空闲岸线区间终点坐标；第三行：区间内可用岸桥数量；第四行：区间内岸桥起始编号
                    
                    
                    useful_space{1, k} = idle_space{1, k}; % 可用岸线区间信息
                    num_use = 1;
                    for j = 1 : size(idle_space{1, k}, 2)  % 空闲岸线段数量
                        % 根据 空闲区间长度是否能满足船身长度 以及 空闲区间内岸桥数量是否满足船舶最小岸桥数量 来确定是否为可用岸线区间
                        if idle_space{1, k}(2, j) - idle_space{1, k}(1, j) < L(vessel_ready{1, k}(1, 1)) || idle_space{1, k}(3, j) < C_min(vessel_ready{1, k}(1, 1))   % 若岸线长度小于船长或区间内岸桥数量小于该船最小岸桥数量，则删除该段
                            useful_space{1, k}(:, num_use) = [];
                        else
                            num_use = num_use + 1;
                        end
                    end
                    
                elseif size(vessel_now{1, k}, 2) == 0
                    useful_space{1, k} = [];
                    useful_space{1, k}(1, 1) = 0;
                    useful_space{1, k}(2, 1) = L_m(k);
                    useful_space{1, k}(3, 1) = N_q(k);
                    useful_space{1, k}(4, 1) = 1;
                    
                elseif size(vessel_now{1, k}, 2) == 1
                    idle_space{1, k} = zeros(4, size(vessel_now{1, k}, 2) + 1);
                    % idle_space{1, k} = [];
                    idle_space{1, k}(1, 1) = 0;
                    idle_space{1, k}(2, 1) = vessel_now{1, k}(4, 1);
                    idle_space{1, k}(3, 1) = vessel_now{1, k}(7, 1) - 1;
                    idle_space{1, k}(4, 1) = min(1, idle_space{1, k}(3, 1));
                    idle_space{1, k}(1, 2) = vessel_now{1, k}(4, 1) + vessel_now{1, k}(5, 1);
                    idle_space{1, k}(2, 2) = L_m(k);
                    idle_space{1, k}(3, 2) = N_q(k) - vessel_now{1, k}(7, 1) - vessel_now{1, k}(6, 1) + 1;
                    idle_space{1, k}(4, 2) = vessel_now{1, k}(7, 1) + vessel_now{1, k}(6, 1);
                    useful_space{1, k} = idle_space{1, k}; % 可用岸线区间信息
                    num_use = 1;
                    for j = 1 : size(idle_space{1, k}, 2)  % 空闲岸线段数量
                        % 根据 空闲区间长度是否能满足船身长度 以及 空闲区间内岸桥数量是否满足船舶最小岸桥数量 来确定是否为可用岸线区间
                        if idle_space{1, k}(2, j) - idle_space{1, k}(1, j) < L(vessel_ready{1, k}(1, 1)) || idle_space{1, k}(3, j) < C_min(vessel_ready{1, k}(1, 1))   % 若岸线长度小于船长或区间内岸桥数量小于该船最小岸桥数量，则删除该段
                            useful_space{1, k}(:, num_use) = [];
                        else
                            num_use = num_use + 1;
                        end
                    end
                    
                end
                % 至此，可用区间信息已得，下一步进行基因修复
                
                % 判断是否存在可用岸线区间，若存在，则不修改靠泊时间，同时对靠泊位置、岸桥数量及编号等信息进行更新；若不存在，则需对靠泊时间以及离港时间进行更新，然后重新修复
                if size(useful_space{1, k}, 2) > 0 % 存在可用的岸线区间时，不需要修改靠泊时间
                    r = 0; % 用于表示待调度船舶停靠的空闲岸线区间编号
                    
                    % 判断是否需要修复靠泊位置
                    is_change_position = 0;  % 用于表示是否需要修复靠泊位置
                    for j = 1 : size(useful_space{1, k}, 2)
                        if vessel_ready{1, k}(4, 1) >= useful_space{1, k}(1, j) && vessel_ready{1, k}(4, 1) + vessel_ready{1, k}(5, 1) <= useful_space{1, k}(2, j)  % 判断当前调度船舶的靠泊位置是否在可用区间内
                            is_change_position = 0; % 不需要改变靠泊位置
                            r = j;  % r为当前调度船舶所停留的岸线区间编号
                            break  % 同时退出for循环
                        else
                            is_change_position = is_change_position + 1;  % 需要改变靠泊位置
                        end
                    end
                    % 至此，我们可以得到是否需要修改调度船舶的靠泊位置
                    
                    % 这里在进行岸线区间的选择时，是否要不采用随机选择的方式，而是用最靠近算法结果的岸线区间
                    
                    % 若需要修改，则对靠泊位置进行修复，若不需要，则进入下一步
                    if is_change_position > 0  % 当需要改变靠泊位置时
                        r = randi([1, size(useful_space{1, k}, 2)], 1, 1); % 随机选择一个可用岸线区间
                        vessel_ready{1, k}(4, 1) = randi([useful_space{1, k}(1, r), useful_space{1, k}(2, r) - L(vessel_ready{1, k}(1, 1))], 1, 1); % 在选择的可用岸线区间范围内随机选择一个可行靠泊位置
                    end
                    % 至此，调度船舶的靠泊位置已得到修复
                    
                    % 这里在进行岸桥数量的确定时，是否能够直接用最大岸桥数量与可用岸桥数量的较小值
                    % 这里进行起始岸桥的确定时，是否可以直接用可用岸桥集合里编号最小的岸桥
                    
                    % 确定岸桥数量以及起始岸桥编号是否需要修复
                    if vessel_ready{1, k}(6, 1) > useful_space{1, k}(3, r)   % || vessel_ready{1, k}(6, 1) < C_min(vessel_ready{1, k}(1, 1)) || vessel_ready{1, k}(6, 1) > C_max(vessel_ready{1, k}(1, 1)) % 判断当前船舶分配到的岸桥数量是否大于该区间内可用的岸桥数量，同时判断分配到的岸桥数量是否在岸桥数量范围内
                        % 若岸桥数量大于区间内岸桥数量，则需要重新分配一个可用的可用的岸桥数量，并重新确定起始岸桥编号
                        vessel_ready{1, k}(6, 1) = randi([C_min(vessel_ready{1, k}(1, 1)),  min(useful_space{1, k}(3, r), C_max(vessel_ready{1, k}(1, 1)))], 1, 1); % 重新分配一个能够区间能够提供的岸桥数量
                        vessel_ready{1, k}(7, 1) = randi([useful_space{1, k}(4, r), useful_space{1, k}(4, r) + useful_space{1, k}(3, r) - vessel_ready{1, k}(6, 1)], 1, 1); % 然后调整起始岸桥编号
                    elseif vessel_ready{1, k}(7, 1) < useful_space{1, k}(4, r) || vessel_ready{1, k}(7, 1) + vessel_ready{1, k}(6, 1) > useful_space{1, k}(4, r) + useful_space{1, k}(3, r)  % 判断岸桥起始编号是否符合要求
                        % 区间内的岸桥能满足分配的岸桥数量，但起始岸桥编号不对时，对起始岸桥编号进行更新
                        vessel_ready{1, k}(7, 1) = randi([useful_space{1, k}(4, r), useful_space{1, k}(4, r) + useful_space{1, k}(3, r) - vessel_ready{1, k}(6, 1)], 1, 1);
                    end
                    % 至此，岸桥数量以及起始岸桥编号均已得到修复，下一步，更新离港时间
                    
                    % 更新离港时间 
                    vessel_ready{1, k}(3, 1) = round(vessel_ready{1, k}(2, 1) + (W_u(vessel_ready{1, k}(1, 1)) + W_d(vessel_ready{1, k}(1, 1))) / (vessel_ready{1, k}(6, 1) * (v_m(k) - vessel_ready{1, k}(9, 1)) * g^(vessel_ready{1, k}(6, 1) - 1)), 1); % 对当前调度船舶离港时刻进行更新
                    
                    % 更新已调度船舶集合以及待调度船舶集合
                    vessel_already{1, k} = [vessel_already{1, k}, vessel_ready{1, k}(:, 1)]; % 将已分配完成的船舶移动到已调度船舶集合中
                    vessel_ready{1, k}(:, 1) = []; % 将已分配完成的船舶从待调度船舶集合中移除
                    % 至此，当前调度船舶已完成调度，接下来进行下一个待调度船舶的调度
                    
                % 当不存在可用岸线区间时，需要将当前调度船舶的靠泊时间推迟到已靠泊船舶中最早离港时刻，然后对下一个待调度船舶进行调度
                else
                    % first_leave_time = min(vessel_now{1, k}(3, :));   % 找到最早离港时刻以及其所对应的索引
                    vessel_ready{1, k}(2, 1) = min(vessel_now{1, k}(3, :));  % 将正在调度船舶的靠泊时间设为最早离港时刻
                   
                    % 离港时刻更新
                    vessel_ready{1, k}(3, 1) = round(vessel_ready{1, k}(2, 1) + (W_u(vessel_ready{1, k}(1, 1)) + W_d(vessel_ready{1, k}(1, 1))) / (vessel_ready{1, k}(6, 1) * (v_m(k) - vessel_ready{1, k}(9, 1)) * g^(vessel_ready{1, k}(6, 1) - 1)), 1);
                    
                    % 对下一个靠泊时刻最早待调度船舶进行调度
                end
            end
            
            for j = 1: size(vessel_already{1, k}, 2)
                % 判断是否为待调度船舶，若编号小于0，则为计划期开始时已靠泊船舶
                if vessel_already{1, k}(1, j) > 0
                    % 染色体的第一段已完成更新（见第二个for循环处），第一段：靠泊码头
                    % pop_repaired(i, vessel_already{1, k}(1, j)) = k;
                    % 对染色体的第二段进行更新，第二段：靠泊位置
                    pop_repaired(i, vessel_already{1, k}(1, j) + N_v) = vessel_already{1, k}(4, j);
                    % 对染色体的第三段进行更新，第三段：靠泊时刻
                    pop_repaired(i, vessel_already{1, k}(1, j) + 2*N_v) = vessel_already{1, k}(2, j);
                    % 对染色体的第四段进行更新，第四段：岸桥数量
                    pop_repaired(i, vessel_already{1, k}(1, j) + 3*N_v) = vessel_already{1 ,k}(6, j);
                    % 对染色体的第五段进行更新，第五段：起始岸桥编号
                    pop_repaired(i, vessel_already{1, k}(1, j) + 4*N_v) = vessel_already{1, k}(7, j);
                    % 对染色体的第六段进行更新，第六段：离港时刻
                    pop_repaired(i, vessel_already{1, k}(1, j) + 5*N_v) = vessel_already{1, k}(3, j);
                    % 第七段：到港时间松弛参数
                    pop_repaired(i, vessel_already{1, k}(1, j) + 6*N_v) = vessel_already{1, k}(8, j);
                    % 第八段：岸桥工作效率松弛参数
                    pop_repaired(i, vessel_already{1, k}(1, j) + 7*N_v) = vessel_already{1, k}(9, j);
                    
                end
            end
        end
        % 至此，码头一的所有待调度船舶均已修复完成，下一步需要对种群中的信息进行更新
            
        % 对种群中停靠在码头一的船舶信息进行更新

        % 至此完成种群中停靠在码头 k 的所有船舶信息的更新
        
        % 加一个判断vessel_already总数量是否为待调度船舶与已调度船舶总数
        total_number = 0;
        for k = 1 : N_m
            total_number = total_number + size(vessel_already{1, k}, 2);
        end
        if total_number ~= N_v + N_0
            disp('-------------error!!!!------------');
        end 
    end 
end
