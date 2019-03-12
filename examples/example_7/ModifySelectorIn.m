function ModifySelectorIn(x)
% Modify text file containing soil hydraulic parameters

global EXAMPLE_dir

% Open SELECTOR.IN
fid_2 = fopen([EXAMPLE_dir '\H1D\SELECTOR.IN'],'r+');

% Go to water flow section
flag = [];
while isempty(flag)
	str = fgetl(fid_2);
	flag = findstr(str,'thr');
end

% Insert water flow parameters
fseek(fid_2,0,'cof');
fprintf(fid_2,' %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e\r\n',x(1:6));
fseek(fid_2,0,'eof');

% Close file
fclose(fid_2);