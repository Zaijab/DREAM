function [SimRR] = hymodFORTRAN(x,Extra); 
% Runs the FORTRAN HYMOD model and returns the simulated discharge

global DREAM_dir EXAMPLE_dir

% Go to the subdirectory of the HYMOD model
cd(EXAMPLE_dir)

% Write the parameter values to a file Param.in
dlmwrite('Param.in',x,'delimiter',' ');

% Execute the model -- this model reads the current parameter values from Param.in
dos('HYMODsilent.exe');

% Load the output of the model 
SimRR=load('Q.out'); SimRR = Extra.F * SimRR(65:Extra.MaxT,1);

% Return to the main directory 
cd(DREAM_dir)