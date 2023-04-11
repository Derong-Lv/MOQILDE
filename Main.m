clear all
clc

rep_all=[];

lognplot1(); % for plotting lognormal PDF - SOLAR PLANT 1
lognplot2(); % for plotting lognormal PDF - SOLAR PLANT 2
lognplot3(); % for plotting lognormal PDF - SOLAR PLANT 3

for Run=1:1

    rtemp=[];
    frep=[];

    pmax=[80 50 35 60 45 1.1 1.1 1.1 1.1 1.1 1.1 ];
    pmin=[20 0 10 0 0 0.95 0.95 0.95 0.95 0.95 0.95];
    Pmax=repmat(pmax,1,1);
    Pmin=repmat(pmin,1,1);
    Pimax=Pmax;
    Pimin=Pmin;
    [~,D]=size(Pmax);
    N=100; 
    ncon=4; 
    M=2;           

    Maxgen=150;
    gen=0;
    Max_FEs=N*Maxgen;
    FEs=0;
    
    x = initialize_variables(N,D,Pimin,Pimax);
    parent = lm_NEW(x,M,ncon,Pimax,Pimin); % Power flow requires Matpower
    FEs=FEs+N;

    parent  = Caculate_constraint_value(parent,D,M,ncon);
    
    while FEs < Max_FEs

        [u,gen]=MOQILDE(gen,N,D,Pimax,Pimin,parent,M);
        child = lm_NEW(u,M,ncon,Pimax,Pimin);
        FEs=FEs+N;


        child  = Caculate_constraint_value(child,D,M,ncon);
        
        Rep=[parent;child];
        [~, size_c]=size(Rep);
        fea=find(Rep(:,size_c)<=0);
        if length(fea)<N
            [~,c]=sort(Rep(:,size_c));
            Rep=Rep(c(1:N),:);
        else
            Rep=Rep(fea,:);
            Rep=selection_new(Rep,D,M,N,gen,Maxgen);
        end
        parent=Rep;
        
        if sum(Rep(:,end)) == 0
            frep=[Rep;frep];
            frep_temp = frep(:,1:D);
            [frep1,frep_r,frep_c] = unique(frep_temp,'rows');
            frep = frep(frep_r,:);
        end
    end
    
    Rep=Non_domination_sort(frep,D,M);
    if size(Rep,1)>N 
        Rep=Crowding_distance(Rep,D,M,N);
    end
    
    rep_all=[rep_all;Rep];
    
    objpareto = Rep(:,D+1:D+M);
    objpareto = flip(sortrows(objpareto,1));
    [no_par,no_pf] = size(objpareto);
    miu_pf = zeros(no_par,no_pf);
    for jj = 1:no_pf
        max_pareto = max(objpareto(:,jj));
        min_pareto = min(objpareto(:,jj));
        miu_pf(:,jj) = (repmat(max_pareto,no_par,1)-objpareto(:,jj))./(max_pareto-min_pareto);
    end
    miu_sum = sum(miu_pf,2);
    miu_wt = miu_sum./sum(miu_sum);
    [miu_max,ind] = max(miu_wt);
    comp_obj(Run,1:(no_pf+1+no_pf)) = [objpareto(ind,:),miu_max,min(objpareto)];
    [ind_comp,~] = find(Rep(:,D+1)==objpareto(ind,1) & Rep(:,D+M)==objpareto(ind,M));
    
    comp_rep(Run,:)=Rep(ind_comp(1),:); 
    
    [~,ind1]=min(Rep(:,D+1));
    best_cost(Run,:)=Rep(ind1,:);
    
    [~,ind2]=min(Rep(:,D+M));
    best_emiss(Run,:)=Rep(ind2,:);
    
    ob1(Run,:)=objpareto(:,1)';
    ob2(Run,:)=objpareto(:,2)';
    
    PopObj=[ob1(Run,:)' ob2(Run,:)'];
    score(Run) = HV(PopObj);
   
    disp(['HV value£º',num2str(score(Run))]); 
end

TVs=sum(rep_all(:,end));
disp(['TVs£º',num2str(TVs)]);
