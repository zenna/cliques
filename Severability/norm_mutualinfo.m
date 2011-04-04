function [ I ] = norm_mutualinfo( X,Y )
%NORM_MUTUALINFO Calculates the normalised mutual information of two
%partitions
%   Follows exactly the definition given in "Comparing community structure
%   identification", by Danon, Diaz-Guilera, Duch, and Arenas
%   Takes two vectors, each of which specifies the partition that a node
%   belongs in. Labelling of communities is not important. Labelling of
%   nodes is critical though.
%
%   e.g. X = [1 1 1 2 2] would mean that nodes 1-3 belong to community 1
%   and nodes 4-5 to community 2.
%
%   X is the "real partition" while Y is the "found partition", but the
%   order does not really matter
%
%   Copyright (c) 2010 Yun William Yu

%Confusion matrix
X_comm = unique(X);
Y_comm = unique(Y);
N = zeros(length(X_comm),length(Y_comm));
for i=1:length(X_comm)
    if (mod(i,10)==0)
        if i~=10
            fprintf('\b\b\b\b')
        else
            fprintf('.')
        end
        fprintf('%4d',i)    % Progress bar
    end
    for j=1:length(Y_comm)
        % N(i,j) is the number of nodes in the real community i that are
        % also in the found community j
        N(i,j)=length(intersect(find(Y==Y_comm(j)),find(X==X_comm(i))));
    end
end

N2 = zeros(size(N));
for i=1:length(X_comm)
    for j=1:length(Y_comm)
        if N(i,j)==0
            N2(i,j)=0;
        else
            N2(i,j) = N(i,j)*log((N(i,j)*sum(sum(N)))/(sum(N(i,:))*sum(N(:,j))));
        end
    end
end
I_numerator = -2 * sum(sum(N2));

N3 = zeros(1,length(X_comm));
for i=1:length(X_comm)
    N3(i)=sum(N(i,:))*log(sum(N(i,:))/sum(sum(N)));
end

N4 = zeros(1, length(Y_comm));
for j=1:length(Y_comm)
    N4(j)=sum(N(:,j))*log(sum(N(:,j))/sum(sum(N)));
end
I_denominator = sum(N3) + sum(N4);

I = I_numerator/I_denominator;


end

