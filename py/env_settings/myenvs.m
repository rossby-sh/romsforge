% MATLAB
% Add path

% restoredefaultpath

setenv('OMP_NUM_THREADS', '1')

% === Define PATH ===

% Custom path (Temporary path)
% PATH for roms_roms c_initial_NPZD.m
addpath('/home/shjo/ROMS/src/remap/d_roms_roms/')


% --- pth01 -> COAST Rutergus tools ---------------------------------------  
pth01 = '/home/shjo/ROMS/src/super_obs/COAWST/Tools/mfiles/rutgers/';

% --- pth02 -> Old netcdf tools (mexcdf) ----------------------------------  
pth02 = '/home/shjo/ROMS/src/super_obs/netcdf_Tools/mexcdf/';


% === Add Path ===
addpath(genpath(pth01))
addpath(genpath(pth02))


disp('--- Myenv set ---')
