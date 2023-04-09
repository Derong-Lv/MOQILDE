function x = initialize_variables(N,D,Pimin,Pimax)
  for i=1:N
     for j=1:D
         x(i,j)=Pimin(j)+rand*(Pimax(j)-Pimin(j));
     end
  end  
end