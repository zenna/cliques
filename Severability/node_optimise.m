
max_size = 256;

%tcommrun=zeros(16,480);
%tcommsize=zeros(16,1);
%tsev=zeros(16,1);
for time=[1:1:12]
	disp({'time:',time});
	[commrun commsize sev]=sev_node(1, adj, time, max_size);
	disp({'commsize:',commsize,'sev:',sev})
	tcommrun(time,:)=commrun;
	tcommsize(time)=commsize;
	tsev(time)=sev;
end

save multilevelr256k12_16_64_256.mat
