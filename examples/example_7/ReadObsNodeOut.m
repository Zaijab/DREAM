function [ModPred] = ReadObsNodeOut
% Read text file containing the time series of simulated soil water contents

global EXAMPLE_dir

% Open OBS_NODE.OUT
fid_2 = fopen([EXAMPLE_dir '\H1D\OBS_NODE.OUT']);

% Go to data section
flag = [];
while isempty(flag)
	str = fgetl(fid_2);
	flag = findstr(str,'time');
end

% Read simulated soil water contents
cols = 4;
rows = Inf;
data = fscanf(fid_2,'%f',[cols rows]);

% Store simulated soil water contents in structure
data = data';
ModPred.hoy = data(:,1);
ModPred.water = data(:,3);

% close file
fclose(fid_2);