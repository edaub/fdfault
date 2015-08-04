function output = load_output(probname, name, datadir)

% function to load output data from simulation data

    if nargin == 2
        datadir = [pwd '/'];
    end

    olddir = pwd;
    cd(datadir);
    eval([probname '_' name]);
    cd(olddir);

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

    if strcmp(field,'vx')
        f = fopen([ datadir probname '_' name '_vx.dat'], 'rb');
        output.vx = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'vy')
        f = fopen([ datadir probname '_' name '_vy.dat'], 'rb');
        output.vy = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'vz')
        f = fopen([ datadir probname '_' name '_vz.dat'], 'rb');
        output.vz = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'sxx')
        f = fopen([ datadir probname '_' name '_sxx.dat'], 'rb');
        output.sxx = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'sxy')
        f = fopen([ datadir probname '_' name '_sxy.dat'], 'rb');
        output.sxy = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'sxz')
        f = fopen([ datadir probname '_' name '_sxz.dat'], 'rb');
        output.sxz = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'syy')
        f = fopen([ datadir probname '_' name '_syy.dat'], 'rb');
        output.syy = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'syz')
        f = fopen([ datadir probname '_' name '_syz.dat'], 'rb');
        output.syz = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'szz')
        f = fopen([ datadir probname '_' name '_szz.dat'], 'rb');
        output.szz = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'lambda')
        f = fopen([ datadir probname '_' name '_lambda.dat'], 'rb');
        output.lambda = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'gammap')
        f = fopen([ datadir probname '_' name '_gammap.dat'], 'rb');
        output.gammap = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'V')
        f = fopen([ datadir probname '_' name '_V.dat'], 'rb');
        output.V = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'U')
        f = fopen([ datadir probname '_' name '_U.dat'], 'rb');
        output.U = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'Vx')
        f = fopen([ datadir probname '_' name '_Vx.dat'], 'rb');
        output.Vx = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'Vy')
        f = fopen([ datadir probname '_' name '_Vy.dat'], 'rb');
        output.Vy = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'Vz')
        f = fopen([ datadir probname '_' name '_Vz.dat'], 'rb');
        output.Vz = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'Ux')
        f = fopen([ datadir probname '_' name '_Ux.dat'], 'rb');
        output.Ux = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'Uy')
        f = fopen([ datadir probname '_' name '_Uy.dat'], 'rb');
        output.Uy = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    elseif strcmp(field,'Uz')
        f = fopen([ datadir probname '_' name '_Uz.dat'], 'rb');
        output.Uz = squeeze(reshape(fread(f,nt*nx*ny*nz,'float64',endian),[nz ny nx nt]));
        fclose(f);
    end

end