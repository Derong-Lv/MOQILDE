function lognplot2()
global mcarlo SP2 nbins
mu = 5.2; sigma = 0.6;
Psr = 60;
Gstd = 1200;
Rc = 140;
nbins = 50;
mcarlo = 10000; 
G1 = lognrnd(mu,sigma,mcarlo,1);
G1und = G1(G1<=Rc);
G1over = G1(G1>Rc);
P1und = Psr*(G1und ./Gstd);
P1over_site = Psr*(G1over./Gstd); 
P1over = min(repmat(Psr,size(G1over,1),1),Psr*(G1over./Gstd));
SP1_site = vertcat(P1und,P1over_site);
SP2 = vertcat(P1und,P1over);
