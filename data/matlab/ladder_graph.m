clear all;
close all;

addpath('ext/stability_cluster/');

a = 10;
b = 10;
c = 4;
d = 0.7

Q = [0, c, a, a, 0, 0, 0, 0;
     c, 0, a, a, 0, 0, 0, 0;
     a, a, 0, c, d, 0, 0, 0;
     a, a, c, 0, 0, d, 0, 0];

Q2 = [0, c, a, a, 0, 0, 0, 0;
     c, 0, a, a, 0, 0, 0, 0;
     a, a, 0, c, d, 0, 0, 0;
     a, a, c, 0, 0, d, 0, 0;
     0, 0, d, 0, 0, c, a, b];

Q = vertcat(Q,flipud(fliplr(Q)))

%Q = round(rand(15,15))
%Q = triu(Q)
%Q.*(1 - diag(ones(1,length(Q))))

incident = []
for j = 1:length(Q)
    for i = 1:length(Q)
        if (Q(j,i) > 0)
            incident = vertcat(incident, [j i])
        end
    end
end
 
[s_all,curve,maxcurve,numb_comp,maxcurveFrac,numb_comp_Frac] =  AllTimeScales(Q,8,10,100,'output.dat');    