function [solv_method,itmax,tol,restart] = read_solver(filename)

ifile = fopen(filename,'r');
C = textscan(fgetl(ifile),'%s'); solv_method = C{1}{1};
C = textscan(fgetl(ifile),'%f'); itmax       = C{1};
C = textscan(fgetl(ifile),'%f'); tol         = C{1};
C = textscan(fgetl(ifile),'%f'); restart     = C{1};
fclose(ifile);

switch lower(solv_method)
    case 'pcg'

    case 'stat_amg'

    case 'bicgstab'

    case 'gmres'

    otherwise
       err_msg = ['Wrong value for solv_method in ' filename];
       error(err_msg);
end

end
