
path(path,'image-analysis/');
path(path,'image-analysis/jc/');

I=imread('Isi164c.jpg');
%I2=rgb2gray(I);
I2=I(:,:,1);
[W, D] = im2mat(I2,20,20,1);

push_ahead=-1;
radius = 175;
time =24;
orad=1; % Overlap radius

outputname = 'Isi164ct24.mat';

	adj_matrix=W;
	%realcomm=load(realcomm_file);


graph_size=length(adj_matrix);
commrun=sparse(graph_size,graph_size); % Contains all the communities found

covered=zeros(1,graph_size);  % List of all the nodes that have already been covered
node=zeros(1,graph_size);    % Node associated with "i"th entry
commsize=zeros(1,graph_size);   % Sizes of the communities (redundant with commrun)
sev=zeros(1,graph_size);    % Severabilities of the communities found

curr_entry = 0;
for i=1:graph_size;
    
    % Figure out the coverpatch
    coverpatch = zeros(1+2*orad,1+2*orad);
    for c_i = 1:(1+2*orad)
        shifted=(i-(orad+1-c_i)*size(I2,1));
        coverpatch(c_i,:)=shifted-orad:shifted+orad;
    end
    coverpatch=coverpatch(:);
    coverpatch((coverpatch<1))=[];
    coverpatch((coverpatch>graph_size))=[];


    if sum(covered(coverpatch))==0
        curr_entry=curr_entry+1;
        node(curr_entry)=i;
        [comm commsize(curr_entry) sev(curr_entry)]=sev_node(i,adj_matrix,time,push_ahead,radius);

        % Put entries in covered
        if commsize(curr_entry)~=graph_size
            for j=1:(commsize(curr_entry))
                covered(comm(j))=covered(comm(j))+1;
            end
        end
    
        commrun(curr_entry,:)=comm;
        fprintf('\t: %5d %5d %5d %4f \n', ...
        [sum((covered~=0)) i commsize(curr_entry) sev(curr_entry)]);

    end    
end;
% Resize all of the arrays
commrun=commrun(1:curr_entry,:);
node=node(1:curr_entry);
commsize=commsize(1:curr_entry);
sev=sev(1:curr_entry);

%% Partition generation
%{
% Create a partition from the communities found
%graph_size=1000;
foundcomm = zeros(1,graph_size);
for i=1:graph_size
    if (mod(i,10)==0)
        if i~=10
            fprintf('\b\b\b\b\b')
        else
            fprintf('.')
        end
        fprintf('%5d',i)    % Progress bar
    end
    for j=1:size(commrun,1)
        if ismember(i,commrun(j,:))
            foundcomm(i)=j;
            break
        end
        foundcomm(i)=-i;
    end
end
%}
%%

% First make commrun_small smaller again, taking out all communities that aren't used at all
%commrun_small=commrun(unique(foundcomm(foundcomm>0)),:);

	save(outputname)

