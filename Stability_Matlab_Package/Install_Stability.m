

disp('Building of combinatorial laplacian mex files...');

cd('Combinatorial');

mex stability_louvain_LCL.cpp community.cpp graph_binary.cpp;

cd ..;

disp('Building of normalised laplacian mex files...');

cd('Normalised');

mex stability_louvain_LNL.cpp community.cpp graph_binary.cpp;

cd ..;

disp('Moving build files to bin directory...');

cd('Combinatorial');
files=dir('stability_louvain_LCL.mex*');
for i=1:length(files)
    copyfile(files(i).name,'../bin');
end

cd('../Normalised');
files=dir('stability_louvain_LNL.mex*');
for i=1:length(files)
    copyfile(files(i).name,'../bin');
end

disp('Adding bin directory to path...');

cd('../bin');
path(path,pwd);
savepath;
cd ..;

disp('Done');