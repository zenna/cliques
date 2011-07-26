function [] = Install_Stability(clique_path)

% path to be provided is the cxx folder..
if ~isdir(clique_path)
    error('The path provided is not a folder. Please provide the path of the source code directory of Clique.');
end

curr_folder = pwd;

cd(clique_path);

disp('Building of mex files...');

mex -v -DUSE_BOOST CXX=g++-4.6 CXXFLAGS="\$CXXFLAGS -std=gnu++0x" -O -output stability_louvain -I./ ./cliques/louvain_matlab_interface.cpp

disp('Moving build files to bin directory...');

movefile('louvain_matlab_interface.mex*',curr_folder);

cd(curr_folder);

disp('Adding bin directory to path...');

path(pwd,path);

savepath;

disp('Done');
