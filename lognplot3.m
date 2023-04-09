function lognplot3()
global mcarlo SP3 nbins
mu = 5.2; sigma = 0.6;
Psr = 45; 
Gstd = 900;
Rc = 110; 
nbins = 50;
mcarlo = 10000;
G1 = lognrnd(mu,sigma,mcarlo,1);
G1und = G1(G1<=Rc);
G1over = G1(G1>Rc);
P1und = Psr*(G1und ./Gstd);
P1over_site = Psr*(G1over./Gstd);
P1over = min(repmat(Psr,size(G1over,1),1),Psr*(G1over./Gstd));
SP1_site = vertcat(P1und,P1over_site);
SP3 = vertcat(P1und,P1over);