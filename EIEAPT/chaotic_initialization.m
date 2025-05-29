function [population,PF] = chaotic_initialization(pop_size, f_num,dim, map_type, lower_bound, upper_bound ,fun)
    % pop_size: 种群大小
    % dim: 每个个体的维�?
    % map_type: 选择的混沌映射类�? ('logistic', 'tent', 'chebyshev', etc.)
    % lower_bound: 每个维度的下�? (向量)
    % upper_bound: 每个维度的上�? (向量)

    % 初始化种群矩�?
    population = zeros(pop_size, dim);

    % 初始混沌�?
    chaotic_seq = rand(1, dim);  % 初始化为随机�?

    % 根据map_type选择不同的混沌映�?
    switch map_type
        case 'logistic'
            r = 3.999;  % Logistic映射的参�?
            for i = 1:pop_size
                chaotic_seq = r .* chaotic_seq .* (1 - chaotic_seq);  % Logistic映射公式
                population(i, :) = chaotic_seq;
            end

        case 'tent'
            for i = 1:pop_size
                chaotic_seq = tent_map(chaotic_seq);  % Tent映射公式
                population(i, :) = chaotic_seq;
            end

        case 'chebyshev'
            n = 2;  % Chebyshev映射的参�?
            for i = 1:pop_size
                chaotic_seq = cos(n * acos(chaotic_seq));  % Chebyshev映射公式
                population(i, :) = chaotic_seq;
            end

        % 叻加其他混沌映射类型

        otherwise
            error('Unsupported chaotic map type');
    end

    % 将混沌序列缩放到指定的上下界范围
    for j = 1:dim
        population(:, j) = lower_bound(j) + (upper_bound(j) - lower_bound(j)) .* population(:, j);
    end
    for i=1:pop_size
      [population(i,dim+1:dim+f_num), PF]= object_fun(population(i,:),f_num,dim,fun);
    end
end

function x = tent_map(x)
    % Tent映射公式
    for i = 1:length(x)
        if x(i) < 0.5
            x(i) = 2 * x(i);
        else
            x(i) = 2 * (1 - x(i));
        end
    end
end