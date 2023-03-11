
% 输入：一条染色体，即一个行向量，输出：各码头的待分配船舶集合
function vessel_ready = find_vessel_ready (pop_repaired_i)
    
    global N_v L N_m
    
    vessel_ready = cell(1, N_m);  
    
    for m = 1 : N_m
        % 码头 m 的待调度船舶集合
        k = 0;
        for j = 1 : N_v
            if pop_repaired_i(j) == m
                
                k = k + 1;
                
                % 第一行：待分配船舶编号
                vessel_ready{1, m}(1, k) = j;
                % 第二行：待分配船舶的靠泊时间
                vessel_ready{1, m}(2, k) = pop_repaired_i(j+2*N_v);
                % 第三行：离港时刻
                vessel_ready{1, m}(3, k) = pop_repaired_i(j+5*N_v);
                % 第四行：靠泊位置
                vessel_ready{1, m}(4, k) = pop_repaired_i(j+N_v);
                % 第五行：船身长度
                vessel_ready{1, m}(5, k) = L(j);
                % 第六行：岸桥数量
                vessel_ready{1, m}(6, k) = pop_repaired_i(j+3*N_v);
                % 第七行：岸桥起始编号
                vessel_ready{1, m}(7, k) = pop_repaired_i(j+4*N_v);
                % 第八行：到港时刻松弛
                vessel_ready{1, m}(8, k) = pop_repaired_i(j+6*N_v);
                % 第九行：岸桥效率松弛
                vessel_ready{1, m}(9, k) = pop_repaired_i(j+7*N_v);
                
            end
        end
    end
end