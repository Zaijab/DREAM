function ModifyProfileDat(Extra);
% Modify text file containing initial condition

global EXAMPLE_dir

% Open PROFILE.DAT
fid_2 = fopen([EXAMPLE_dir '\H1D\PROFILE.DAT'],'r+');

% Go to profile section
flag = [];
while isempty(flag)
	str = fgetl(fid_2);
	flag = findstr(str,'Mat');
end

% Insert new initial condition
fseek(fid_2,0,'cof');
fprintf(fid_2,'%5.0f %15.6e %15.6e %3.0f %3.0f %15.6e %15.6e %15.6e %15.6e\n',Extra.initial');
fseek(fid_2,0,'eof');

% Close file
fclose(fid_2);