function [] = Install_Stability(clique_path)

if ~isdir(clique_path)
    error('The path provided is not a folder. Please provide the path of the source code directory of Clique.');
end

curr_folder = pwd;

cd(clique_path);

disp('Building of mex files...');

mex -DUSE_BOOST -O -output stability_louvain -I./ ./louvain_matlab_interface.cpp

disp('Moving build files to bin directory...');

movefile('louvain_matlab_interface.mex*',curr_folder);

cd(curr_folder);

disp('Adding bin directory to path...');

path(pwd,path);

savepath;

disp('Done');
