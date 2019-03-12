function NC_Flag=NonConvergenceTest()
%Go through MODFLOW list file to check for non-convergence
NCfid = fopen('Twri.lst','r');      % ... open list file
eofstat = feof(NCfid);
while eofstat == 0                              % ... read file and search for respective message
    string = fgetl(NCfid);                        % ... read one line
    nonconv_fail = findstr('FAILED',string);         % ... search for "non-convergence message"
    nonconv_dry = findstr('DRY(',string);         % ... search for "non-convergence message"
    if ~isempty(nonconv_fail)                        % ... if message was found
        NC_Flag=1;                                  % Error is inf. if the model did not converge
        return
    end
    if ~isempty(nonconv_dry)                        % ... if message was found
        NC_Flag=1;                                  % Error is inf. if the model did not converge
        return
    end    
eofstat = feof(NCfid);
NC_Flag=0;
end
