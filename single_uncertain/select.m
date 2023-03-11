
% 选择函数

function pop_select = select(pop_old, proba_choose)
    
    global popsize N_v
    
    pop_select = zeros(popsize, 8*N_v);
    
    for j = 1 : popsize
        r_select = rand(); % 确定一个概率，用于判断选择哪一个个体进入下一阶段
        % disp(r_select);
        % disp(proba_choose(k-1));
        
        z = 1;
        for k = 2 : size(proba_choose, 2)
            z = k;
            if r_select >= proba_choose(z-1) && r_select < proba_choose(z)
                break
            end
        end
        % 确定了选择哪一个个体进入下一个阶段
        pop_select(j, :) = pop_old(z, :);
    end
end
