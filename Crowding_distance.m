function f=Crowding_distance(chromosome,dim,n_obj,max_chromosome_size) %(rep,V,M,NP)

[~, ~] = size(chromosome);

f=[];
max_rank=chromosome(max_chromosome_size+1,dim+n_obj+1);
previous_index = 0;
%current_index = 0;
for i=1:max_rank-1
    current_index = find(chromosome(:,dim + n_obj + 1) == i, 1, 'last' );
    f(previous_index + 1 : current_index, :) = chromosome(previous_index + 1 : current_index, :);
    previous_index = current_index;
end

remaining = max_chromosome_size - previous_index;

if remaining==0
   % previous_index = current_index;
else No_temp= chromosome(:,dim+n_obj+1)==max_rank; %rank-1
    temp_chromosome = chromosome(No_temp , :); 
    fxmax=max(temp_chromosome(:,dim+1:dim+n_obj),[],1); 
    fxmin=min(temp_chromosome(:,dim+1:dim+n_obj),[],1); 
    [nr,~]=size(temp_chromosome);
    dd=zeros(nr,1);
    for k=1:n_obj
        [y,yy]=sort(temp_chromosome(:,dim+k));
        dd(yy(1),1)=10000;dd(yy(nr),1)=10000; 
        if length(fxmax)<2 
        range = fxmax;
        else 
        range =(fxmax(k)-fxmin(k)); 
        end
        if range==0,range=.00001;end
        for i=2:nr-1 
                dd(yy(i),1)=dd(yy(i),1)+(y(i+1)-y(i-1))/range;
        end
    end
    [~,ss]=sort(-dd);
    f=[f;temp_chromosome(ss(1:remaining)',:)];
end