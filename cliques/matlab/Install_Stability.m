function [] = Install_Stability(clique_path)

% path to be provided is the cxx folder..
if ~isdir(clique_path)
    error('The path provided is not a folder. Please provide the path of the source code directory of Clique.');
end

curr_folder = pwd;

cd(clique_path);

disp('Compiling files...');

% Comile with :
!g++-4.6 -c  -I./ -I/usr/local/MATLAB/R2010b/extern/include -I/usr/local/MATLAB/R2010b/simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread -std=gnu++0x  -DUSE_BOOST -DMX_COMPAT_32 -O3 -DNDEBUG  ./cliques/louvain_matlab_interface.cpp -lemon -L/usr/local/include/lemon

!gfortran -c  -I./ -I/usr/local/MATLAB/R2010b/extern/include -I/usr/local/MATLAB/R2010b/simulink/include -fexceptions -fPIC -fno-omit-frame-pointer  -DUSE_BOOST -DMX_COMPAT_32 -O3  /home/mts09/repositories/group_repository/graph-codes/cliques/cxx/cliques/expokit.f

!g++-4.6 -O3 -pthread -shared -Wl,--version-script,/usr/local/MATLAB/R2010b/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -o  stability_louvain.mexa64  louvain_matlab_interface.o expokit.o -Wl,-rpath-link,/usr/local/MATLAB/R2010b/bin/glnxa64 -L/usr/local/MATLAB/R2010b/bin/glnxa64 -lmx -lmex -lmat -lm -lblas -llapack -lemon -lgfortran


% setenv('LD_RUN_PATH','/home/mts09/repositories/group_repository/graph-codes/cliques/cxx/cliques')
% 
% mex -v -DUSE_BOOST CXXFLAGS="\$CXXFLAGS -std=gnu++0x" -O -output stability_louvain -I./ ...
%         ./cliques/louvain_matlab_interface.cpp ...
%         CXX=g++-4.6 -lblas -llapack -lemon -lgfortran...
%     /home/mts09/repositories/group_repository/graph-codes/cliques/cxx/cliques/libexpokit.so
%/home/mts09/repositories/group_repository/graph-codes/cliques/cxx/cliques/expokit.f ...

disp('Moving build files to bin directory...');

movefile('stability_louvain.mex*',curr_folder);

% cd(curr_folder);
% 
% disp('Adding bin directory to path...');
% 
% path(pwd,path);
% 
% savepath;

disp('Done');


