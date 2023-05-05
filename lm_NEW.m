function f = lm_NEW(Q,M,ncon,Pimax,Pimin)
% function [x,f,error] = lm_NEW(Q)
[N,D]=size(Q);
K = D+M+3;% 维度 目标数 排序 
f = zeros(N,K+ncon);% 个数维度 目标数  排序  惩罚数

%%%%%%%%潮流计算%%%%%%%%
global SP1 SP2 SP3 nbins mcarlo 

%%%火力机组%%%
% bus_no. alpha beta gama omega miu d e Pmin POZL1 POZH1 POZL2 POZH2
 thgendata = [
 	1	0.04091 -0.05554 0.06490 0.000200 6.667 18 0.037 50 0 0 0 0;
 	2	0.02543 -0.06047 0.05638 0.000500 3.333 16 0.038 20 30 40 55 65;
 	8	0.05326 -0.03550 0.03380 0.002000 2.000 12 0.045 10 0 0 0 0;];

for i=1:N
    for k=1:D
            if Q(i,k) > Pimax(k)
                Q(i,k) = Pimax(k);
            elseif Q(i,k) < Pimin(k)
                Q(i,k) = Pimin(k);
            else
                Q(i,k)=Q(i,k);
            end
    end

% %%%禁止运行区%%%
% % for i=1:N
%     %%%系统1%%%
%     if Q(i,1)>thgendata(2,10) && Q(i,1)<thgendata(2,11)
%         temp = round((Q(i,1)-thgendata(2,10))/(thgendata(2,11)-thgendata(2,10)));
%         Q(i,1) = thgendata(2,10);
%         if temp == 1
%             Q(i,1) = thgendata(2,11);
%         end
%     end
%     if Q(i,1)>thgendata(2,12) && Q(i,1)<thgendata(2,13)
%         temp = round((Q(i,1)-thgendata(2,12))/(thgendata(2,13)-thgendata(2,12)));
%         Q(i,1) = thgendata(2,12);
%         if temp == 1
%             Q(i,1) = thgendata(2,13);
%         end
%     end
end

%%%潮流%%%
for i=1:N
data = loadcase('case30_S');
data.gen(2:6,2)= Q(i,1:5);%功率 %火 风 火 光 风
data.gen(1:6,6) = Q(i,6:11);%电压偏差

mpopt = mpoption('pf.enforce_q_lims',0,'verbose',0,'out.all',0);
result = runpf(data,mpopt);
QH1(i)=result.gen(1,2)';%平衡节点1功率

S1(:,i)=result.branch(:,14);
S2(:,i)=result.branch(:,15);
S3(:,i)=result.branch(:,16);

VII(:,i) = result.bus(:,8);

QG(:,i) = result.gen(:,3)/data.baseMVA;
end

%%%火电
for i=1:N
thpowgen(i,:) = [QH1(i),Q(i,1),Q(i,3)];%火 光 火 光 光

thgencoeff = vertcat(data.gencost(1:2,5:7),data.gencost(4,5:7));

thgencost(i) = sum(thgencoeff(:,1)+thgencoeff(:,2).*thpowgen(i,:)'+thgencoeff(:,3).*(thpowgen(i,:).^2)');
end


for i=1:N
%%%光伏电站1
sgenpar = [1   11  1.60];
Crsj = 3; % Reserve cost for solar power overestimation ($/MW)
Cpsj = 1.4; % Penalty cost for solar power underestimation ($/MW)
schspow = [Q(i,2)]; % solar generator schedule power %火 风 火 光 风
% Segregate over and underestimated power on the power histogram
[~,JJ]=size(schspow);
[histy1,histx1] = hist(SP1,nbins);
for j=1:JJ
Lowind1(j,:)= histx1<schspow(j);
Highind1(j,:) = histx1>schspow(j);
end
for j=1:JJ
allP1und = schspow(j)-histx1(histx1<schspow(j));
allP1over = histx1(histx1>schspow(j))-schspow(j);
ProbP1und = histy1(Lowind1(j,:))./mcarlo;
ProbP1over = histy1(Highind1(j,:))./mcarlo;
% Finding under and over estimation cost
C1und(j) = sum(Crsj*(ProbP1und.*allP1und));
C1over(j) = sum(Cpsj*(ProbP1over.*allP1over));
end
sovundcost = [C1und;C1over];
sgencosta(i,:)= sgenpar(3)*schspow+sum(sovundcost,1); % solar generator cost
sgencost1(i)=sum(sgencosta(i,:),2)';

%%%光伏电站2
sgenpar = [1   11  1.60];
Crsj = 3; % Reserve cost for solar power overestimation ($/MW)
Cpsj = 1.4; % Penalty cost for solar power underestimation ($/MW)
schspow = [Q(i,4)]; % solar generator schedule power %火 风 火 光 风
% Segregate over and underestimated power on the power histogram
[~,JJ]=size(schspow);
[histy1,histx1] = hist(SP2,nbins);
for j=1:JJ
Lowind1(j,:)= histx1<schspow(j);
Highind1(j,:) = histx1>schspow(j);
end
for j=1:JJ
allP1und = schspow(j)-histx1(histx1<schspow(j));
allP1over = histx1(histx1>schspow(j))-schspow(j);
ProbP1und = histy1(Lowind1(j,:))./mcarlo;
ProbP1over = histy1(Highind1(j,:))./mcarlo;
% Finding under and over estimation cost
C1und(j) = sum(Crsj*(ProbP1und.*allP1und));
C1over(j) = sum(Cpsj*(ProbP1over.*allP1over));
end
sovundcost = [C1und;C1over];
sgencosta(i,:)= sgenpar(3)*schspow+sum(sovundcost,1); % solar generator cost
sgencost2(i)=sum(sgencosta(i,:),2)';

%%%光伏电站3
sgenpar = [1   11  1.60];
Crsj = 3; % Reserve cost for solar power overestimation ($/MW)
Cpsj = 1.4; % Penalty cost for solar power underestimation ($/MW)
schspow = [Q(i,5)]; % solar generator schedule power %火 风 火 光 风
% Segregate over and underestimated power on the power histogram
[~,JJ]=size(schspow);
[histy1,histx1] = hist(SP3,nbins);
for j=1:JJ
Lowind1(j,:)= histx1<schspow(j);
Highind1(j,:) = histx1>schspow(j);
end
for j=1:JJ
allP1und = schspow(j)-histx1(histx1<schspow(j));
allP1over = histx1(histx1>schspow(j))-schspow(j);
ProbP1und = histy1(Lowind1(j,:))./mcarlo;
ProbP1over = histy1(Highind1(j,:))./mcarlo;
% Finding under and over estimation cost
C1und(j) = sum(Crsj*(ProbP1und.*allP1und));
C1over(j) = sum(Cpsj*(ProbP1over.*allP1over));
end
sovundcost = [C1und;C1over];
sgencosta(i,:)= sgenpar(3)*schspow+sum(sovundcost,1); % solar generator cost
sgencost3(i)=sum(sgencosta(i,:),2)';
end

%%%约束%%%
for i=1:N
%功率约束%
PGSmax = data.gen(:,9)';
PGSmin = data.gen(:,10)';
PGS(i,:) = [QH1(i),Q(i,1:5)];
% PJZ(i) = Q(i,1)>thgendata(2,12) && Q(i,1)<thgendata(2,13)+Q(i,1)>thgendata(2,10) && Q(i,1)<thgendata(2,11);
for j=1:6
PG(i,j) = (PGS(i,j)<PGSmin(j))*(abs(PGSmin(j)-PGS(i,j))/(PGSmax(j)-PGSmin(j)))+(PGS(i,j)>PGSmax(j))*(abs(PGSmax(j)-PGS(i,j))/(PGSmax(j)-PGSmin(j)));
end
% PGSerr(i)=sum([PG(i,:),PJZ(i)])';%越界惩罚
PGSerr(i)=sum([PG(i,:)])';%越界惩罚

%线路功率约束%
blimit = data.branch(:,6);%支路长期允许功率
Slimit(:,i) = sqrt(S1(:,i).^2+S2(:,i).^2);
Serr(i) = sum((Slimit(:,i)>blimit).*abs(blimit-Slimit(:,i)))/data.baseMVA;

%电压约束%
genbus = data.gen(:,1);
VI=VII(:,i);
VI(genbus)=[];
Vmax=data.bus(:,12);
Vmin=data.bus(:,13);
Vmax(genbus)=[];
Vmin(genbus)=[];
VIerr(i) = sum((VI<Vmin).*(abs(Vmin-VI)./(Vmax-Vmin))+(VI>Vmax).*(abs(Vmax-VI)./(Vmax-Vmin)));
VD(i)=sum(abs(VI-1));%电压偏差

%无功功率约束%
Qmax = data.gen(:,4)/data.baseMVA;
Qmin = data.gen(:,5)/data.baseMVA;
Qerr(i) = sum((QG(:,i)<Qmin).*(abs(Qmin-QG(:,i))./(Qmax-Qmin))+(QG(:,i)>Qmax).*(abs(Qmax-QG(:,i))./(Qmax-Qmin)));

%阀点效应%
valveff(i) = sum(abs(thgendata(:,7).*sin(thgendata(:,8).*(thgendata(:,9)-thpowgen(i,:)')))); % if all have valve effects


end
% ploss = sum(result.branch(:,14)+result.branch(:,16));%网损
% loss=get_losses(result)全部损耗（实虚）real imag

ploss = sum(S1+S3);%网损

error = [Qerr' VIerr' Serr' PGSerr'];%error = [Qerr;VIerr;Serr;PGSerr];——罚函数

%%%%%成本计算%%%%%
cumcost=thgencost+valveff+sgencost1+sgencost2+sgencost3;

% %%%%污染排放计算%%%%
emission = sum(thgendata(:,2)+thgendata(:,3).*thpowgen'/100+thgendata(:,4).*(thpowgen.^2/100^2)'...
     +thgendata(:,5).*exp(thgendata(:,6).*thpowgen'/100));


     z1 = cumcost';%成本 %%罚函数导致搜寻点太少
     z2 = emission';%污染排放

% z3 = ploss';
 
f(:,1:D)=Q;
f(:,D+1:D+M)=[z1,z2];
f(:,D+M+4:end)=error;
% f=[Q,z1',z2']; %包含两个目标值
% p=cumcost'; % 单目标

%  x=Q;
  end