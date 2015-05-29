function front = load_front(probname, iface, datadir)

% function to load output data from simulation data

    if nargin == 2
        datadir = [pwd '/'];
    end

    olddir = pwd;
    cd(datadir);
    eval([probname '_front_' num2str(iface)]);
    cd(olddir);

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
        output.z = squeeze(reshape(fread(f,nx*ny,'float64',endian),[ny nx]));
        fclose(f);
    end

    f = fopen([ datadir probname '_front_' num2str(iface) '_t.dat'], 'rb');
    front.t = squeeze(reshape(fread(f,nx*ny,'float64',endian),[ny nx]));
    fclose(f);

end