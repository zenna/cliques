function list= expmat_revised(time)

load('A.mat');
m=size(A,1);

diagdeg=sparse((diag(sum(A)))/sum(sum(A)));
trans=sparse(diag(    (sum(A)).^(-1)     ) * A);
clear A;
Id = eye(m);
Lap=sparse(trans-Id);
clear trans;
clear Id;

exponential=sparse(expm(time.*Lap));
clear Lap;
solution=sparse(diagdeg*exponential);
solution=round(1e16*solution)/1e16;
clear exponential;
clear diagdeg;

solution=full(solution);

list=adj2edj(solution);

clear solution


