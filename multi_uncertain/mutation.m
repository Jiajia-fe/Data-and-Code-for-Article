
% 变异操作

function pop_muta = mutation (pop_cross, p_m)
    
    global N_v L_m L C_min C_max N_q W_u W_d v_m g A_sigma v_sigma  A ED  
    
    pop_muta = pop_cross;
    
    % 进行变异操作，均匀变异
    
    for j = 1 : N_v
        % 判断每个基因位是否需要变异
        is_muta = rand();
        
        if is_muta < p_m
            pop_muta(j) = randi(3); % 对靠泊码头进行变异
            pop_muta(j+N_v) = randi([0, L_m(pop_muta(j)) - L(j)], 1, 1); % 对靠泊位置进行变异
            pop_muta(j+6*N_v) = randi([-A_sigma(j) * 20, A_sigma(j) * 20 ], 1, 1)/10;
            pop_muta(j+7*N_v) = randi([-20 * v_sigma(j), 20 * v_sigma(j) ], 1, 1)/10;
            pop_muta(j+2*N_v) = randi([max(A(j) * 10 + pop_muta(j + 6*N_v)*10, 0), round(ED(j)*10, 0)], 1, 1)/10; % 对靠泊时间进行变异    % 10*ED(j) - 10*round((W_u(j) + W_d(j))/(C_max(j) * v_m(pop_muta(j)) * g^(C_max(j)- 1)), 1)
            pop_muta(j+3*N_v) = randi([C_min(j), C_max(j)], 1, 1); % 对岸桥数量进行变异
            pop_muta(j+4*N_v) = randi([1, N_q(pop_muta(j)) - pop_muta(j + 3*N_v) + 1], 1, 1); % 对岸桥起始编号进行变异
            pop_muta(j+5*N_v) = round(pop_muta(j + 2*N_v) + (W_u(j) + W_d(j))/(pop_muta(j+3*N_v) * (v_m(pop_muta(j)) - pop_muta(j+7*N_v)) * g^(pop_muta(j+3*N_v) - 1)), 1); % 对离港时刻进行更新

        end
    end
    
end
