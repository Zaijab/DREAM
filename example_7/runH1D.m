function [ModPred] = runH1D(x,Extra)
% Run the HYDRUS-1D model with new parameters and initial and boundary conditions and get simulated time series of soil water contents

global DREAM_dir EXAMPLE_dir

% Modify SELECTOR.IN
ModifySelectorIn(x)

% Modify PROFILE.DAT
Extra.initial(:,3) = x(7);
ModifyProfileDat(Extra);

% Modify ATMOSPH.IN
Extra.boundcon(:,7) = x(7);
ModifyAtmosphIn(Extra);

% Change directory
cd(EXAMPLE_dir)

% Run HYDRUS-1D
[status,output] = dos('H1D_CALC.EXE');

% Change directory
cd(DREAM_dir)

% Read OBS_NODE.OUT
ModPred = ReadObsNodeOut;