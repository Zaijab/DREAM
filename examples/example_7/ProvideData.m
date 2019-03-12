function [Extra] = ProvideData;
% Load observational data and define initial and boundary conditions

global EXAMPLE_dir

% Provide observational data
load meas.mat
Extra.water = meas.water(meas.waterind);
Extra.hoy = meas.hoy(meas.waterind);
Extra.ObsNode = 1;

% Provide data needed to modify initial condition
load ProfileDat.mat
Extra.initial = initial_conditions;

% Provide data needed to modify the lower boundary condition
load AtmosphIn.mat
Extra.boundcon = boundcon;

% Modify level_01.dir (needed to execute HYDRUS-1D)
fid = fopen([EXAMPLE_dir '\level_01.dir'],'w+');
fprintf(fid,'%s',[EXAMPLE_dir '\H1D']);
fclose(fid);