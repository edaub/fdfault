% example using load_output function in MATLAB

% required arguments are problem name and output unit name
% data directory is optional, if no argument provided assumes it is the current working directory

vybody = load_output('hpctest','vybody');

% loads data structure containing information

vybody

% because of differences in how MATLAB and C++ order arrays, field arrays are indexed by (z,y,x,t)
% any singleton dimensions are removed

% plot velocity

pcolor(vybody.x, vybody.y, vybody.vy(:,:,4));
shading flat;
axis image;
colorbar;
