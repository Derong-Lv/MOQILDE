function f = Non_domination_sort(x,dim,n_obj) %(frep,V,M)
xx=x; %copy back before crowd_dis_sort
x=xx;     
[N, ~] = size(x);
clear m
        
% Initialize the front number to 1.
front = 1;

% There is nothing to this assignment, used only to manipulate easily in
% MATLAB.
F(front).f = [];
individual = [];

for i = 1 : N 
    % Number of individuals that dominate this individual
    individual(i).n = 0; 
    % Individuals which this individual dominate
    individual(i).p = []; 
    for j = 1 : N 
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        for k = 1 : n_obj
            if (x(i,dim + k) < x(j,dim + k))
                dom_less = dom_less + 1;
            elseif (x(i,dim + k) == x(j,dim + k))
                dom_equal = dom_equal + 1;
            else 
                dom_more = dom_more + 1;
            end
        end
        if dom_less == 0 && dom_equal ~= n_obj 
            individual(i).n = individual(i).n + 1; 
        elseif dom_more == 0 && dom_equal ~= n_obj 
            individual(i).p = [individual(i).p j]; 
        end
    end   
    if individual(i).n == 0 
        xx(i,n_obj + dim + 1) = 1; %rank-1
        F(front).f = [F(front).f i]; %front-1
    end
end
% Find the subsequent fronts
while ~isempty(F(front).f)
   Q = [];
   for i = 1 : length(F(front).f) 
       if ~isempty(individual(F(front).f(i)).p)
        for j = 1 : length(individual(F(front).f(i)).p) 
            	individual(individual(F(front).f(i)).p(j)).n = individual(individual(F(front).f(i)).p(j)).n - 1; 
         if individual(individual(F(front).f(i)).p(j)).n == 0 
               		xx(individual(F(front).f(i)).p(j),n_obj + dim + 1) = front + 1; %rank-2
                    Q = [Q individual(F(front).f(i)).p(j)]; 
         end
        end
      end
   end
   front =  front + 1;
   F(front).f = Q; %front-2
end

[~,index_of_fronts] = sort(xx(:,n_obj + dim + 1)); 

for i = 1 : length(index_of_fronts)
    sorted_based_on_front(i,:) = xx(index_of_fronts(i),:); 
end
current_index = 0;

% Find the crowding distance for each individual in each front
for front = 1 : (length(F) - 1) 
    y = [];
    previous_index = current_index + 1;
    for i = 1 : length(F(front).f)
        y(i,:) = sorted_based_on_front(current_index + i,:);
    end
    current_index = current_index + i;
    for i = 1 : n_obj
        [~, index_of_objectives] = sort(y(:,dim + i)); 
        sorted_based_on_objective = [];
        for j = 1 : length(index_of_objectives) 
            sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
        end
        f_max = sorted_based_on_objective(length(index_of_objectives), dim + i); 
        f_min = sorted_based_on_objective(1, dim + i); 
        
        y(index_of_objectives(length(index_of_objectives)),n_obj + dim + 1 + i)= Inf; 
        y(index_of_objectives(1),n_obj + dim + 1 + i) = Inf; 
         for j = 2 : length(index_of_objectives) - 1 
            next_obj  = sorted_based_on_objective(j + 1,dim + i);
            previous_obj  = sorted_based_on_objective(j - 1,dim + i);
            if (f_max - f_min == 0) 
                y(index_of_objectives(j),n_obj + dim + 1 + i) = Inf;
            else
                y(index_of_objectives(j),n_obj + dim + 1 + i) = (next_obj - previous_obj)/(f_max - f_min);
            end
         end
    end
    distance = [];
    distance(:,1) = zeros(length(F(front).f),1);
    for i = 1 : n_obj
        distance(:,1) = distance(:,1) + y(:,n_obj + dim + 1 + i); 
    end
    y(:,n_obj + dim + 2) = distance;
    z(previous_index:current_index,:) = y;
end
f = z();