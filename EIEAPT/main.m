% Clear workspace and initialize variables
clear all
clc
tic;
ea = 1;
F1daxiao = [];

% Parameter settings
fun = 'LSMOP1';  % Function selection
pop = 50;        % Population size for chaotic mapping

% Additional parameters
maxFE = 100000;  % Maximum function evaluations
f_num = 3;       % Number of objectives
x_num = 500;     % Number of decision variables
fz = x_num * 0.1; % Divide decision variables into 10 groups
fzx_num = x_num / fz;

% Boundary settings for LSMOP
x_min = zeros(1, x_num);
x_max = [ones(1, f_num-1), 10.*ones(1, x_num-f_num+1)];

% Genetic algorithm parameters
pc = 1;         % Crossover probability
pm = 1/x_num;   % Mutation probability
yita1 = 20;     % SBX parameter
yita2 = 20;     % Polynomial mutation parameter

fr = 0.1;       % Frequency to call reference vector
JH = 10;

% Generate reference vectors
[Zs, KK] = UniformPoint(100, f_num);
for i = 1:size(Zs, 1)
    Zs(i,:) = Zs(i,:)./norm(Zs(i,:));
end
Z = Zs;

% Calculate neighboring angle for angle normalization
cosineVV = Z*Z';
[scosineVV, neighbor] = sort(cosineVV, 2, 'descend');
acosVV = acos(scosineVV(:,2));
refV = (acosVV);

%% Initialize population
map_type = 'logistic'; % Choose Logistic mapping
[chromo, PF] = chaotic_initialization(pop, f_num, x_num, map_type, x_min, x_max, fun);

i = 1;
PTchromo = chromo(:,1:x_num);
for i = 1:pop
    PTchromo(pop+i,:) = x_max - PTchromo(i,:);
end

%% Initialize weights
x_min1 = zeros(1, fzx_num);
x_max1 = ones(1, fzx_num);
w_parent = initialize1(2*pop, f_num, fzx_num, x_min1, x_max1);
rew_parent = w_parent;

for j = 1:fzx_num
    chromo_problemt(1:2*pop, (j-1)*fz+1:j*fz) = w_parent(1:2*pop,j) .* PTchromo(1:2*pop, (j-1)*fz+1:j*fz);
end

%% Evaluate transformed values
for i = 1:2*pop
    [chromo_problemt(i, x_num+1:x_num+f_num), ~] = object_fun(chromo_problemt(i,1:x_num), f_num, x_num, fun);
    chromo_problemt(i, x_num+f_num+2) = w_parent(i, fzx_num+1);
end

%% Main loop
% SPS: Set penalty parameters for each reference vector
for i = 1:size(Z, 1)
    theta11(i) = exp(0.5*(max(Z(i,:))-min(Z(i,:))));    
end

Gene = 1; % Counter for population information exchange

while ea < maxFE
    chromo_parent = chromo;
    
    %% Non-dominated sorting to select non-dominated individuals
    [F1, chromo_non] = non_domination_sort(size(chromo_parent), chromo, f_num, x_num);
    F1daxiao = [F1daxiao, size(F1(1).ss, 2)];
    
    %% Generate offspring based on direction vectors
    chromo_offspring = cross_mutation6(chromo_parent, f_num, x_num, x_min, x_max, pc, pm, yita1, yita2, fun, F1);
    ea = ea + size(chromo_offspring, 1);
    
    %% Binary crossover and polynomial mutation for weights
    w_parent1 = F_mating(w_parent);
    w_offspring = cross_mutation1(w_parent1, fzx_num, x_min1, x_max1, pc, pm, yita1, yita2);
    chromo_problemt_offspring = [];
    
    for i = 1:size(w_offspring, 1)
        for j = 1:fzx_num
            chromo_problemt_offspring(i, (j-1)*fz+1:j*fz) = w_offspring(i,j) .* PTchromo(w_offspring(i,fzx_num+1), (j-1)*fz+1:j*fz);
        end
        [chromo_problemt_offspring(i, x_num+1:x_num+f_num), ~] = object_fun(chromo_problemt_offspring(i,1:x_num), f_num, x_num, fun);
        chromo_problemt_offspring(i, x_num+f_num+2) = w_offspring(i, fzx_num+1);
    end
    ea = ea + size(w_offspring, 1);
    
    %% Population information exchange
    rechromo_problemt_offspring = [];
    
    if rem(Gene, JH) == 0
        less_than_50 = w_parent(w_parent(:,fzx_num+1) <= 50,:);
        more_than_50 = w_parent(w_parent(:,fzx_num+1) > 50,:);
        unique_wless = unique(less_than_50(:,fzx_num+1))';
        unique_wmore = unique(more_than_50(:,fzx_num+1))';
        unique_w = unique(w_parent(:,fzx_num+1))';
        numbers = 1:pop;
        complementNumbers = setdiff(numbers, unique_wless);
        numbers1 = pop+1:2*pop;
        complementNumbers1 = setdiff(numbers1, unique_wmore);
        numbers2 = 1:2*pop;
        complementNumbers2 = setdiff(numbers2, unique_w);
        w_reparent = initialize2(complementNumbers2, fzx_num, x_min1, x_max1);
        
        if size(chromo_parent,1) >= size(complementNumbers,2)
            random_indices1 = randperm(size(chromo_parent,1));
            selected_numbers1 = random_indices1(1:size(complementNumbers,2));
            for i = 1:size(complementNumbers,2)
                PTchromo(complementNumbers(i),1:x_num) = chromo_parent(selected_numbers1(i),1:x_num);
            end
        else
            for i = 1:size(chromo_parent,1)
                PTchromo(complementNumbers(i),1:x_num) = chromo_parent(i,1:x_num);
            end
        end
        
        if size(chromo_parent,1) >= size(complementNumbers1,2)
            random_indices1 = randperm(size(chromo_parent,1));
            selected_numbers1 = random_indices1(1:size(complementNumbers1,2));
            for i = 1:size(complementNumbers1,2)
                PTchromo(complementNumbers1(i),1:x_num) = x_max - chromo_parent(selected_numbers1(i),1:x_num);
            end
        else
            for i = 1:size(chromo_parent,1)
                PTchromo(complementNumbers1(i),1:x_num) = x_max - chromo_parent(i,1:x_num);
            end
        end
        
        for i = 1:size(w_reparent,1)
            for j = 1:fzx_num
                rechromo_problemt_offspring(i,(j-1)*fz+1:j*fz) = w_reparent(i,j) .* PTchromo(w_reparent(i,fzx_num+1),(j-1)*fz+1:j*fz);
            end
            [rechromo_problemt_offspring(i,x_num+1:x_num+f_num), ~] = object_fun(rechromo_problemt_offspring(i,1:x_num), f_num, x_num, fun);
            rechromo_problemt_offspring(i,x_num+f_num+2) = w_reparent(i,fzx_num+1);
        end
        
        chromo_problemt_offspring = [chromo_problemt_offspring; rechromo_problemt_offspring];
        w_parent = [w_parent; w_reparent];
    end
    
    %% Combine parent and offspring for dimensionality reduction
    combine_chromo_problemt = [];
    combine_w = [];
    combine_w = [w_parent; w_offspring];
    [pop_chromo_problemt, ~] = size(chromo_problemt);
    [pop_chromo_problemt_offspring, ~] = size(chromo_problemt_offspring);
    combine_chromo_problemt(1:pop_chromo_problemt, 1:(f_num+x_num+2)) = chromo_problemt(:,1:(f_num+x_num+2));
    combine_chromo_problemt((pop_chromo_problemt+1):(pop_chromo_problemt+pop_chromo_problemt_offspring), 1:(f_num+x_num+2)) = chromo_problemt_offspring(:,1:(f_num+x_num+2));
    
    %% Combine parent and offspring
    [pop_parent, ~] = size(chromo);
    [pop_offspring, ~] = size(chromo_offspring);
    combine_chromo(1:pop_parent, 1:(f_num+x_num)) = chromo(:,1:(f_num+x_num));
    combine_chromo((pop_parent+1):(pop_parent+pop_offspring), 1:(f_num+x_num)) = chromo_offspring(:,1:(f_num+x_num));
    
    %% Elite selection guidance
    chromof = [];
    [Selection] = F_select11(combine_chromo(:,(x_num+1):x_num+f_num), Z, theta11);
    chromo = combine_chromo(Selection,:);
    chromof = chromo;
    
    %% Select elite solutions after dimensionality reduction
    chromo_problemt = [];
    [Selection1] = F_select11(combine_chromo_problemt(:,(x_num+1):x_num+f_num), Z, theta11);
    chromo_problemt = combine_chromo_problemt(Selection1,:);
    w_parent = combine_w(Selection1,:);
    
    %% Avoid infinite loops
    A = round(chromo(:,x_num+1:x_num+f_num), 4);
    B = unique(A, 'rows');
    if size(B,1) <= 1
        [chromo, ~] = chaotic_initialization(pop, f_num, x_num, map_type, x_min, x_max, fun);
    end
    
    %% Combine solutions with and without dimensionality reduction as next generation
    if rem(Gene, 10) == 0
        combine_chromo_combine_PT = []; 
        combine_chromo_combine_PT = [chromo(:,1:x_num+f_num); chromo_problemt(:,1:x_num+f_num)];
        [Selection2] = F_select11(combine_chromo_combine_PT(:,(x_num+1):x_num+f_num), Z, theta11);
        chromo111 = combine_chromo_combine_PT(Selection2,:);
        chromo = [];
        chromo = chromo111;
    end
    
    %% Reference vector adaptation
    if rem(Gene, 50) == 0  
        % Update the reference vectors
        Zmin = min(chromo(:,(x_num+1):x_num+f_num),[],1);	
        Zmax = max(chromo(:,(x_num+1):x_num+f_num),[],1);	
        Z = Zs;
        Z = Z.*repmat((Zmax - Zmin)*1.0, size(Z,1), 1);
        
        for i = 1:size(Z,1)
            Z(i,:) = Z(i,:)./norm(Z(i,:));
        end
        
        % Update the neighboring angle value for angle normalization
        cosineVV = Z*Z';
        [scosineVV, neighbor] = sort(cosineVV, 2, 'descend');
        acosVV = acos(scosineVV(:,2));
        refV = (acosVV); 
    end
    
    if mod(i,10) == 0
        fprintf('%d gen has completed!\n',i);
    end
    
    %% Calculate IGD
    Distance = min(pdist2(PF, chromo(:,x_num+1:x_num+f_num)), [], 2);
    IGD = mean(Distance)
    
    cla;
    Draw(chromo(:,x_num+1:x_num+f_num));
    
    %% Calculate HV
    if(f_num == 2)
        plot(PF(1:52,1), PF(1:52,2), 'k-');
        plot(PF(53:end,1), PF(53:end,2), 'k-');
    end
    
    if(f_num == 3)
        TURE_PF(fun, PF);
    end
    
    Gene = Gene + 1;
    pause(0.01);
end

toc;
aaa = toc;