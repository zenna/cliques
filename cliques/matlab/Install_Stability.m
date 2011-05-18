function [] = Install_Stability(clique_path)

if ~isdir(clique_path)
    error('The path provided is not a folder. Please provide the path of the source code directory of Clique.');
end

curr_folder = pwd;

cd(clique_path);

disp('Building of mex files...');

mex -DUSE_BOOST -I./ ./louvain_matlab_interface.cpp

disp('Moving build files to bin directory...');

movefile('louvain_matlab_interface.mex*',curr_folder);

copyfile('./stability.m',curr_folder);

copyfile('./varvarinfo.m',curr_folder);

cd(curr_folder);

mexfile = dir('louvain_matlab_interface.mex*');

movefile(mexfile(1).name,['stability_louvain' mexfile(1).name(25:end)]);

disp('Adding bin directory to path...');

path(pwd,path);

savepath;

disp('Done');