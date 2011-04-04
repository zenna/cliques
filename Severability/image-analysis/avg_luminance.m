%% Calculates the average luminance of a community
% I2 = image to be analyzed

sev_scaled = (sev-min(sev))/(max(sev)-min(sev));

luminance=zeros(size(commsize));
for i=1:length(luminance)
    luminance(i)=sum(I2(commrun(i,1:commsize(i))))/(commsize(i)*255)*sev(i);
%    luminance(i)=sum(I2(commrun(i,1:commsize(i))))* sev(i);
end

clear sev_scaled;

%% Sorts commrun by average luminance
[luminance_sorted IX] = sort(luminance);
commrun_sorted=commrun(IX,:);
commsize_sorted=commsize(IX);
sev_sorted=sev(IX);
clear IX;