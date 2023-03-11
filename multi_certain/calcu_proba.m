
% 根据适应度值进行选择概率的计算

function proba_choose = calcu_proba(fitness_value)
    
    global popsize
    
    proba_single = zeros(1, popsize);
    proba_choose = zeros(1, popsize);
    
    pro_single(1) = fitness_value(1) / sum(fitness_value);
    
    proba_choose(1) = pro_single(1);
    
    for i = 2 : popsize
        
        proba_single(i) = fitness_value(i) / sum(fitness_value);
        
        proba_choose(i) = proba_single(i) + proba_choose(i - 1);
        
    end
    
    % 至此，得到各染色体被选择的概率

end
