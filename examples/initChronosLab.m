% ChronosLab path
ChronosLabPath = '../sources';
addpath(genpath(ChronosLabPath));

global DEBINFO;
% GENERAL PART
DEBINFO.flag = true;
DEBINFO.flag = false;
% PROLONGATION
DEBINFO.prol = [];
% PRINT INFO
DEBINFO.prol.prt_flag = false;
% PRINT UNIT
DEBINFO.prol.ofile = 0;
% PROLONGATION: PRINT ITERATIONS COUNT
DEBINFO.prol.it_print = false;
% PROLONGATION: PRINT LIST OF NEIGHBOURS
DEBINFO.prol.neigh_print = false;

% COARSENING
DEBINFO.coarsen = [];
% PRINT INFO
DEBINFO.coarsen.draw_dist = false;
