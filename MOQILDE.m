function [TempPop,gen]=MOQILDE(gen,m,d,MaxParValue,MinParValue,parent,M)

Population=parent;
parentval= (parent(:,d+1)/max(parent(:,d+1))+parent(:,d+M)/max(parent(:,d+M)))/2;

gen = gen + 1;

for i = 1 : m
    % Generate the mutant "v"
    while true
        r1 = round(m * rand + 0.5);
        if (r1 ~= i), break, end
    end
    while true
        r2 = round(m * rand + 0.5);
        if (r2 ~= r1) && (r2 ~= i), break, end
    end
    while true
        r3 = round(m * rand + 0.5);
        if (r3 ~= r1) && (r3 ~= r2) && (r3 ~= i) , break, end
    end
    
    
    FF = 0.1+0.9*rand;
    CCR = rand;
    
    v(i,:) = Population(r1,:) + FF * (Population(r2,:) - Population(r3,:));
    
    
    X1 = Population(r1,:);
    
%     while true
%         r4 = round(m * rand + 0.5);
%         if (r4 ~= wz && parent(wz1(r4),d+1)<=parent(wz1(m-1),d+1) ), break, end
%     end
%     while true
%         r5 = round(m * rand + 0.5);
%         if (r5 ~= r4 && (r5 ~= wz) && parent(wz1(r5),d+M)<=parent(wz1(m-1),d+M)), break, end
%     end
    
    X2 = Population(r2,:);
    
    X3 = Population(r3,:);
    
    
    val_X1 = parentval(r1,:);
    
    val_X2 = parentval(r2,:);
    
    val_X3 = parentval(r3,:);
    
    fenzi(i,:) = (X2.^2 - X3.^2) * val_X1 + (X3.^2 - X1.^2) * val_X2 + (X1.^2 - X2.^2) * val_X3;
    fenmu(i,:) = (X2 - X3) * val_X1 + (X3 - X1) * val_X2 + (X1 - X2) * val_X3 + 1e-6;
    
    TempPop0(i,:) = 0.5 * fenzi(i,:) ./ fenmu(i,:);
    
    if any(isnan(Population(i,:)))
        Population(i,:) = MinParValue + rand(1,d) .* (MaxParValue - MinParValue);
    end
    if any(isnan(TempPop0(i,:))) || any(isnan(fenzi(i,:))) || any(isnan(fenmu(i,:)))
        TempPop0(i,:) = Population(i,:);
    end
    
    j_rand = ceil(d * rand );
    
    for j = 1 : d
        if (rand < CCR || j == j_rand)
            TempPop(i,j) = v(i,j);
        else
            TempPop(i,j) = TempPop0(i,j);
        end
    end
    
end
end
