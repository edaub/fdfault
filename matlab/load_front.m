function front = load_front(probname, iface, datadir)

    % load_front is a function to load rupture front time data from simulation data
    % function returns a data structure holding the following information:
    %     endian (string) = byte-ordering of the binary data
    %     nx (integer) = number of x grid points
    %     ny (integer) = number of y grid points
    %     nz (integer) = number of z grid points
    %     x (array) = x grid values (shape nz*ny*nx)
    %     y (array) = y grid values (shape nz*ny*nx)
    %     z (array) = z grid values (shape nz*ny*nx)
    %     t (array) = rupture time array (shape nz*ny*nx*nt)
    %
    % If the point in question never satisfies the condition for rupture, the array
    % value for the rupture time is -1.

    [result currentdir] = system('pwd');

    if nargin == 2
        datadir = [deblank(currentdir) '/'];
    end

    if datadir(end) ~= '/'
        datadir = [ datadir '/'];
    end

    cd(datadir);
    system(['cp ' probname '_front_' num2str(iface) '.m tmpfile.m']);
    eval('tmpfile');
    system(['rm -f tmpfile.m']);
    cd(deblank(currentdir));

    front.endian = endian;
    front.nx = nx;
    front.ny = ny;

    f = fopen([ datadir probname '_front_' num2str(iface) '_x.dat'], 'rb');
    front.x = squeeze(reshape(fread(f,nx*ny,'float64',endian),[ny nx]));
    fclose(f);
    f = fopen([ datadir probname '_front_' num2str(iface) '_y.dat'], 'rb');
    front.y = squeeze(reshape(fread(f,nx*ny,'float64',endian), [ny nx]));
    fclose(f);
    try
        f = fopen([ datadir probname '_front_' num2str(iface) '_z.dat'], 'rb');
        front.z = squeeze(reshape(fread(f,nx*ny,'float64',endian),[ny nx]));
        fclose(f);
    end

    f = fopen([ datadir probname '_front_' num2str(iface) '_t.dat'], 'rb');
    front.t = squeeze(reshape(fread(f,nx*ny,'float64',endian),[ny nx]));
    fclose(f);

end