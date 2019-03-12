function ModifyAtmosphIn(Extra)
% Modify text file containing atmospheric boundary condition

global EXAMPLE_dir

% Open ATMOSPH.IN
fid_2 = fopen([EXAMPLE_dir '\H1D\ATMOSPH.IN'],'r+');

% Go to atmospheric information section
flag = [];
while isempty(flag)
	str = fgetl(fid_2);
	flag = findstr(str,'tAtm');
end

% Insert new lower boundary condition
fseek(fid_2,0,'cof');
fprintf(fid_2,'%11.0f %11.2f %11.4f %11.0f %11.0f %11.0f %11.2f %11.0f\n',Extra.boundcon');
fseek(fid_2,0,'eof');

% Close file
fclose(fid_2);

