function xf= selection_new(x,dim,n_obj,NP,gen,Max_Gen)

no_g=Max_Gen;

xx=x;

for i=1:n_obj
    fmax(i)=max(x(:,dim+i));
    fmin(i)=min(x(:,dim+i));
    f(:,i)=(x(:,dim+i)-fmin(i))/(fmax(i)-fmin(i));
end


rank2=sum(f,2);

xx(:,dim+n_obj+1)=rank2;

perc=0.5-0.5*gen/Max_Gen;
[~,k]=min(sqrt(sum((ones(size(f,1),n_obj)*perc-f).^2,2))); 

    findex2=[];

    for i=1:n_obj

        findex1=find(f(:,i)<=f(k,i));
        findex2=[findex2;findex1];  
    end

    findex=unique(findex2);

if length(findex)<NP 
    xx(findex,dim+n_obj+1)=0; 
    [~,index_of_fronts] = sort(xx(:,n_obj + dim + 1));
    index_of_fronts=index_of_fronts(1:NP);
    xf = xx(index_of_fronts,:);

    
else
    xx=xx(findex,:);
    x=xx;

    length(xx);

    tim=1;

    saver=[];

    for i=1:n_obj

         [mini(i),bbb]=min(x(:,dim+i));
         maxi(i)=max(x(:,dim+i));
         previrous=bbb;
         n_grid(i)=(maxi(i)-mini(i))/no_g;
         for j=1:no_g*tim 
             a=find((xx(:,dim+i)>(j-1)*n_grid(i)+mini(i)) & (xx(:,dim+i)<=j*n_grid(i)+mini(i)));
             if (isempty(a))==0 

                     [~,c]=min(xx(a,dim+n_obj+1));

                     saver=[saver,a(c)]; 

             else
                 continue;
             end
         end

           saver=[saver,bbb];
    end  
     saver1=unique(saver);
    length(saver1);
     if length(saver1)>NP

        yy=x(saver1,:);
        index1=randperm(length(saver1));
        xf=yy(index1,:);

     else  
        xx(saver1,dim+n_obj+1)=0;
        [~,index_of_fronts] = sort(xx(:,n_obj + dim + 1));
        index_of_fronts=index_of_fronts(1:NP);
        xf = xx(index_of_fronts,:); 
     end
end
