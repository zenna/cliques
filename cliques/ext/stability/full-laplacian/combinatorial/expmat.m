function list = expmat(edge_list, time)

%Reads an edge list in the variable edge_list
%time = str2double(time);
load(edge_list);
l=size(Graph,1);
Graph(:,1) = Graph(:,1)+1;
Graph(:,2) = Graph(:,2)+1;
m=max(max(Graph(:,1)),max(Graph(:,2)));

A=sparse(m,m);

for i=1:l
   A(Graph(i,1),Graph(i,2))=Graph(i,3);
end

%A=A+A';

%We calculate the matrix X(t) for the normalised laplacian for 100 values
%of time

%degrees=sum(A);
diagdeg=sparse((diag(sum(A)))/sum(sum(A)));  %diag matrix with stat distr
trans=sparse(diag(    (sum(A)).^(-1)     ) * A);  %(stochastic) transition matrix
%stat=(sum(A))/sum(sum(A));  % vector of statio proba
clear A;
Id = eye(m);
Lap=sparse(trans-Id);
%sum(Lap');
%degrees*Lap;
clear trans;
clear Id;


exponential=sparse(expm(time.*Lap));
clear Lap;
solution=sparse(diagdeg*exponential);
clear exponential;
clear diagdeg;


%On reecrit sous forme de liste de liens

limit = max(max(solution))/100000;

count = 0;
for i=1:m
   for j=1:m
      if(solution(i,j)>limit)
          count = count + 1;
      end
   end
end

list=zeros(count,3);


count = 0;
for i=1:m
   for j=1:m
      if(solution(i,j)>limit)
          count = count + 1;
          list(count,1)=i;
          list(count,2)=j;
          list(count,3)=solution(i,j);
      end
   end
end
list(:,1) = list(:,1)-1;
list(:,2) = list(:,2)-1;
% name='expmat_out.txt';
% %dlmwrite(name,list, ' ');
% out_id=fopen(name,'w');
% for i=1:count
%     fprintf(out_id,int2str(list(i,1)));
%     fprintf(out_id,'\t');
%     fprintf(out_id,int2str(list(i,2)));
%     fprintf(out_id,'\t');
%     fprintf(out_id,num2str(list(i,3)));
%     fprintf(out_id,'\t');
%     fprintf(out_id,'\n');
% end
% fclose(out_id);