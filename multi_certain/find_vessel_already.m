
% 输出：已完成调度的各码头船舶集合

function vessel_already = find_vessel_already ()
    
    global N_m N_0 d_0 C_0 Q_first_0 x_0 L_0 M_0
    
    % 各码头船舶信息集合
    vessel_already = cell(1,N_m);
    
    for m = 1 : N_m
        k = 0;
        for j = 1 : N_0
            % 码头 i 的已调度船舶集合
            if M_0(j) == m
                
                k = k + 1;
                
                % 第一行：船舶编号
                vessel_already{1, m}(1, k) = -j;
                % 第二行：靠泊时刻
                vessel_already{1, m}(2, k) = 0;
                % 第三行：离港时刻
                vessel_already{1, m}(3, k) = d_0(j);
                % 第四行：靠泊位置
                vessel_already{1, m}(4, k) = x_0(j);
                % 第五行：船身长度
                vessel_already{1, m}(5, k) = L_0(j);
                % 第六行：岸桥数量
                vessel_already{1, m}(6, k) = C_0(j);
                % 第七行：岸桥起始编号
                vessel_already{1, m}(7, k) = Q_first_0(j);
                % 第八行：到港时刻松弛
                vessel_already{1, m}(8, k) = 0;
                % 第九行：岸桥效率松弛
                vessel_already{1, m}(9, k) = 0;
                
            end
        end
    end
    
end
