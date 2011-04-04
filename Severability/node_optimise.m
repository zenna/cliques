
push_ahead=-1;
radius = 480;

%tcommrun=zeros(16,480);
%tcommsize=zeros(16,1);
%tsev=zeros(16,1);
for time=[1:1:24]
	disp({'time:',time});
	[commrun commsize sev]=sev_node(1, adj, time, push_ahead, radius);
	disp({'commsize:',commsize,'sev:',sev})
	tcommrun(time,:)=commrun;
	tcommsize(time)=commsize;
	tsev(time)=sev;
end

save multilevelr480_20_40_160.mat
