clc;
clear;
close all;

list = {'smoother/MEX_NSY_RFSAI', 'smoother/sym_AFSAI', 'prolong/MEX_EXTI_Prol', ...
        'prolong/EMIN/MEX_EMIN', 'prolong/MEX_BAMG_Prol', 'prolong/MEX_CLAS_prol', ...
        'prolong/MEX_HYBC_prol'};

home = pwd;
for folder = list
    fprintf('Compiling MEX files in %s\n', folder{1});
    cd(folder{1});
    compile
    cd(home);
end
