function [s_all,curve,maxcurve,numb_comp,maxcurveFrac,numb_comp_Frac]=AllTimeScales(Q,longu2,tmax,samplingtime,output_files)
%
%[s_all,curve,maxcurve,numb_comp]=AllTimeScales(Q,longu2,tmax)
%This function constructs a sequence of clusterings 
%which are optimal at different time scales.
%
%INPUTS: 
%Q is the adjacency matrix Q of the graph (square matrix)
%longu2 is the maximum number of clusters it looks for.
%tmax is the maximum time for which we'll compute stability
%
%OUTPUTS:
%s_all is a three-dim (longu x longu2 x longu2) array 
% where longu is the size of Q
% It contains the information about the clusters.  
% s_all(:,:,k) is the membership matrix for the k-way clustering
% s_all(i,j,k)=1 if vertex i belongs to cluster j
% In fact, only s_all(:,k,k) is relevant, 
% since s_all(:,k+1:longu2,k) is identically zero 
%
%curve contains the stability curves. It is a (longu2 X tmax) matrix.
%curve(k,t+1) is the stability at time t for the k-way clustering
%
%
%
%maxcurve is max(curve), hence maxcurve(t+1) contains the optimal stability
%(among the clusterings that are computed) for time t
%
%
%numb_comp(t+1) gives the optimal number of clusters at time t. 
%
%A graph with all the curves and the maximum curve is produced as well
%
% NOTE: By default, this function call cluster_shi_r 
%(=spectral clustering by Shi-Malik), contained in the directory
%SPECTRAL_HOME
%

%samplingtime=100; %change here for the sampling time

global SPECTRAL_HOME
SPECTRAL_HOME = getenv('SPECTRAL_HOME') 
if isempty(SPECTRAL_HOME) 
  SPECTRAL_HOME='C:\Documents and Settings\Stefano\Documenti\Imperial\Matlab\work\TimeScalesFiles\SPECTRAL_HOME';
end

addpath(genpath(SPECTRAL_HOME))


PiGoo=(diag(sum(Q)))/sum(sum(Q));  %diag matrix with stat distr
Trans=diag(    (sum(Q)).^(-1)     ) * Q;  %(stochastic) transition matrix
MM=Q/sum(sum(Q));  %symmetric matrix of proba of edges
piGoo= (sum(Q))/sum(sum(Q));  % vector of statio proba
  
longu=  length(Q);

curve=[];

s_all=zeros(longu,longu2,longu2);

atoms_clusters = zeros (longu,longu2); %ste: Na x Nc matrix. rows number = atom index; column number = clusterization way; elements = clusters index the (row) atom belongs to according to the (column) partition way.
%ste savedata stability01 = zeros (longu2,2); %ste: Nc x 2 matric. contains stabilities S(t=0) and S(t=1) for each clustering way.

for i=2:longu2
    i
        
    %[ind,s,s_orth]=SpecClust2(Q,i);
        
    ind=cluster_shi_r(Q,i,'ncut');  %use this one by default
    atoms_clusters (:,i)= ind'; %ste savedata
    %ind=cluster_single_linkage(Q,i);
        
    %%%%%%%%%%ind=cluster_kvv(Q,i,'scale_row','ncut');
    
    %ind=cluster_kvv(Q,i,'scale_row','conductance'); %used in final tests
    
    %ind=cluster_Newman_Girvan(Q,i);
    
    %ind=MultiGooModuBis(Q,i);    %Newman's spectral algo 
    
    %ind=cluster_ward_linkage(Q,i);
    
    s=zeros(i,longu);
    for p=1:i  %convert ind into s
        s(p,:)=(ind==p);
        
    end;
    s=s*1;  %convert to numerical
    s=s';
        
    
    taille=size(s);
    s_all(1:taille(1),1:taille(2),i)=s;
end;

%Now we compute the stability curves
PiGooMinuspiGoo=PiGoo-piGoo'*piGoo;

disp('power')
Transk=Trans^(samplingtime); 
disp('done')

for i=2:longu2% i=[2,3,4,5,6,7,15,31]        %i=2:longu2 %                  
    i
    s=s_all(:,1:i,i);

    %TransPower=eye(longu);
    
    %ste Transk=Trans^(samplingtime); 
    s_power=s;
     
    for j=1:tmax
        %curve(i,j)=trace(s'*(PiGoo-piGoo'*piGoo)*TransPower*s); %too
        %complex to compute
        curve(i,j)=trace(s'*(PiGooMinuspiGoo)*s_power);
        s_power=Transk*s_power; %sampling k by k
        %TransPower=Transk*TransPower;  %sampling k by k
    end;
    
end;

[maxcurve,numb_comp]=max(curve);


%FRACTIONAL TIMES
%%%%%%%%%%%%%%%%%%%%%%%

%general %plot 

for i=2:longu2
    
%ste    if any(numb_comp==i) 
        ccc=curve(i,:);        
        xaxis=0:(length(ccc)-1);
        %plot(xaxis,ccc,'--');
        %ste hold on
        %figure; %ste
%ste    else %ste
%ste        ste=1 %ste
%ste    end;
end;
    
ccc=maxcurve;
xaxis=0:(length(ccc)-1);
%plot(xaxis,ccc); %ste hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%

curveFrac=zeros(longu2,tmax,2);

for     i=2:longu2
    i
    s=s_all(:,1:i,i);
    modul1=trace(s'*(PiGooMinuspiGoo)*Trans*s)
    modul0=trace(s'*(PiGooMinuspiGoo)*s)
    %pause

    for j=1:tmax
        %we look at times 1/j and we interpolate
        %(1/j)*modul0+(1-1/j)*modul1;
        %pause
        curveFrac(i,j,2)=(1/j)*modul1+(1-1/j)*modul0;
        curveFrac(i,j,1)=1/j; %ste: useless*! it fills a longu2 x tmax matrix of trivial times. Furthermore, each line is the same
    end; %for j
end %for i


%curveFrac(i,j,1) = time = 1/j
%curveFrac(i,j,2) = interpolated stability of clustering i at time 1/j


[maxcurveFrac,numb_comp_Frac]=max(curveFrac(:,:,2));

maxcurveFrac(2,:)=maxcurveFrac;
maxcurveFrac(1,:)=curveFrac(2,:,1);  %times %ste: consequence of useless*


%we extracted which clustering is optimal for every time 1/j. The value of stability for time 1/j is in maxcurveFrac(2,:)
%numb_comp_Frac(j) contains the optimal clustering (=how many clusters) for time 1/j

%general %plot which draws all optimal communities, for fractional and integral times.

for i=2:longu2

    if any(numb_comp==i) | any(numb_comp_Frac==i) %OR for frac times
        ccc=curve(i,:)
        xaxis=0:(length(ccc)-1)
        %plot(xaxis,ccc,'--')
        hold on
        %fract times
        ccc = curveFrac(i,:,2);
        xaxis= curveFrac(i,:,1);
        %plot(xaxis,ccc,'--')
        hold on
        %end of Frac times
    end;
end;

disp('end of computation. timescale generation...')
dataout=zeros(tmax/samplingtime-1,3);
 for j=1:tmax/samplingtime-1
     dataout(j,1)=j*samplingtime;
 end
 disp('done. stabilities...')
 dataout(:,2)=maxcurve(2:tmax/samplingtime)';
 disp('done. cluster number...')
 dataout(:,3)=numb_comp(2:tmax/samplingtime)';
 %maxcurve is max(curve), hence maxcurve(t+1) contains the optimal stability
%(among the clusterings that are computed) for time t
%
%
%numb_comp(t+1) gives the optimal number of clusters at time t.

file_timeVSstabVSclust = [output_files '_JC_timeVSstabVSclust.dat'];
file_clusters = [output_files '_JC_atom_cluster.dat'];
 disp('done. write clusters to file...')

dlmwrite(file_clusters, atoms_clusters, 'delimiter', '\t');

 disp('done. write curves to file...')

dlmwrite(file_timeVSstabVSclust, dataout, 'delimiter', '\t');
