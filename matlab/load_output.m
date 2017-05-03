function output = load_output(probname, name, datadir)
% load_output is a function to load output data from simulation data
% function returns a data structure holding the following information:
%     field (string) = field that was saved from the simulation
%     endian (string) = byte-ordering of the binary data
%     nt (integer) = number of time steps
%     nx (integer) = number of x grid points
%     ny (integer) = number of y grid points
%     nz (integer) = number of z grid points
%     x (array) = x grid values (shape nz*ny*nx)
%     y (array) = y grid values (shape nz*ny*nx)
%     z (array) = z grid values (shape nz*ny*nx)
%     the field data is stored in an array with the name given by the string
%          in the field attribute (shape nz*ny*nx*nt)
%
% All singleton dimensions are flattened in the resulting arrays for the
% grid points and field data

    [result currentdir] = system('pwd');

    if nargin == 2
        datadir = [deblank(currentdir) '/'];
    end
    
    if datadir(end) ~= '/'
        datadir = [ datadir '/'];
    end

    cd(datadir);
    system(['cp ' probname '_' name '.m tmpfile.m']);
    eval('tmpfile');
    system(['rm -f tmpfile.m']);
    cd(deblank(currentdir));

    output.field = field;
    output.endian = endian;
    output.nt = nt;
    output.nx = nx;
    output.ny = ny;
    output.nz = nz;

    if nt > 1
        f = fopen([ datadir probname '_' name '_t.dat'], 'rb');
        output.t = fread(f,[nt],'float64',endian);
        fclose(f);
    end
    f = fopen([ datadir probname '_' name '_x.dat'], 'rb');
    output.x = squeeze(reshape(fread(f,nx*ny*nz,'float64',endian),[nz ny nx]));
    fclose(f);
    f = fopen([ datadir probname '_' name '_y.dat'], 'rb');
    output.y = squeeze(reshape(fread(f,nx*ny*nz,'float64',endian), [nz ny nx]));
    fclose(f);
    try
        f = fopen([ datadir probname '_' name '_z.dat'], 'rb');
        output.z = squeeze(reshape(fread(f,nx*ny*nz,'float64',endian),[nz ny nx]));
        fclose(f);
    end

    f = fopen([ datadir probname '_' name '_' field '.dat'], 'rb');
    eval(['output.' field ' = squeeze(reshape(fread(f,nt*nx*ny*nz,''float64'',endian),[nz ny nx nt]))']);
    fclose(f);

end