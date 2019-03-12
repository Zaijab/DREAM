function J = tomokernel_straight(data,x,z)
%
% This function computes the kernel matrix for a straight ray tomographic 
% inversion given the data matrix and the x and z position vectors for the 
% vertical and horizontal cell *boundaries*.
%
% Syntax:  J = tomokernel_straight(data,x,z)
%
% by James Irving
% December 2005

% check that data are within bounds set by x and z
xmin = x(1);  xmax = x(end);  
zmin = z(1);  zmax = z(end);  
if xmin > min([data(:,1); data(:,3)]) | xmax < max([data(:,1); data(:,3)]) | ...
        zmin > min([data(:,2); data(:,4)]) | zmax < max([data(:,2); data(:,4)]);
    disp('Error:  Data outside of range of min and max values');
    return
end

% determine some initial parameters
dx = x(2)-x(1);                         % horizontal discretization
dz = z(2)-z(1);                         % vertical discretization
xmid = (xmin+dx/2):dx:(xmax-dx/2);		% x-coordinates of cell midpoints
zmid = (zmin+dz/2):dz:(zmax-dz/2);	    % z-coordinates of cell midpoints
nrays = size(data,1);	                % number of rays to consider
nx = length(x)-1;                       % number of cells in x-direction
nz = length(z)-1;                       % number of cells in z-direction

% initialize the sparse storage arrays
maxelem = round(nrays*sqrt(nx^2+nz^2));
irow = zeros(maxelem,1);
icol = zeros(maxelem,1);
jaco = zeros(maxelem,1);

% determine elements of Jacobian matrix
count = 1;
for i=1:nrays
    xs = data(i,1);                      % x-position of source
    xr = data(i,3);                      % x-position of receiver
    zs = data(i,2);						 % z-position of source
    zr = data(i,4);						 % z-position of receiver
    if xs==xr; xr=xr+1e-10; end          % if ray is vertical, add for stability
    slope = (zr-zs)/(xr-xs);		     % slope of raypath
    
    % vector containing x-positions of vertical cell boundaries hit by the ray,
    % and also the ray end points
    xcellb = x(:);
    xcellb = xcellb(xcellb > min([xs,xr]) & xcellb < max([xs,xr]));
    xcellb = [xcellb; xs; xr];
    
    % vector containing z-positions of horizontal cell boundaries
    % and also the ray end points
    zcellb = z(:);
    zcellb = zcellb(zcellb > min([zs,zr]) & zcellb < max([zs,zr]));
    zcellb = [zcellb; zs; zr];
    
    % form matrix containing all intersection points of ray with cell boundaries
    % then sort these points in order of increasing x-coordinate
    ipoint(:,1) = [xcellb; xs + (zcellb-zs)*1/(slope+1e-20)];   % x-coords of all intersection points
    ipoint(:,2) = [zs + (xcellb-xs)*slope; zcellb];             % z-coords of all intersection points
    ipoint = sortrows(ipoint,1);
    
    % calculate length and midpoint of the ray bits between the intersection points
    xlength = abs(ipoint(2:end,1)-ipoint(1:end-1,1));      % x-component of length
    zlength = abs(ipoint(2:end,2)-ipoint(1:end-1,2));      % z component of length
    clength = sqrt(xlength.^2 + zlength.^2);
    cmidpt = 0.5*[(ipoint(1:end-1,1)+ipoint(2:end,1)),(ipoint(1:end-1,2)+ipoint(2:end,2))];
    
    % calculate which slowness cell each ray bit belongs to, and place properly in J matrix
    srow = ceil((cmidpt(:,1)-xmin)./dx);
    scol = ceil((cmidpt(:,2)-zmin)./dz);
    srow(srow<1) = 1;  srow(srow>nx) = nx;
    scol(scol<1) = 1;  scol(scol>nz) = nz;
    njaco = length(srow);
    irow(count:(count+njaco-1)) = i*ones(njaco,1);
    icol(count:(count+njaco-1)) = (scol-1)*nx+srow;
    jaco(count:(count+njaco-1)) = clength;
    count = count + njaco;
    clear ipoint   
end

% convert sparse storage arrays to sparse matrix
index = find(jaco);
irow = irow(index); 
icol = icol(index); 
jaco = jaco(index);
J = sparse(irow,icol,jaco,nrays,nx*nz);
