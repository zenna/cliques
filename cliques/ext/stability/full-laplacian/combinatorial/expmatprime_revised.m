function list= expmatprime_revised(time)

load('A.mat');
m=size(A,1);

diagdeg=sparse(1/m .* eye(m));
Lap=sparse(A-diag(sum(A)));
total=sum(sum(A));
avdegree=total/m;
clear A;

exponential=sparse(expm((time/avdegree).*Lap));
clear Lap;
solution=sparse(diagdeg*exponential);
solution=round(1e16*solution)/1e16;
clear exponential;
clear diagdeg;

solution=full(solution);

list=adj2edj(solution);

clear solution

