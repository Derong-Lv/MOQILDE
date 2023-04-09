function Rep_C = Caculate_constraint_value(Rep,D,M,ncon)
cons = Rep(:,(D+M+4):(D+M+3+ncon));
cons_max = max(cons);
cons_den = 1./cons_max;
cons_den(~isfinite(cons_den)) = 0;
cons_norm = bsxfun(@rdivide,cons,cons_max);
cons_norm(~isfinite(cons_norm)) = 0;
tot_cons = sum(cons_norm,2)./sum(cons_den);
tot_cons(~isfinite(tot_cons)) = 0;
Rep_C = [Rep(:,1:(D+M+3)),tot_cons];
end