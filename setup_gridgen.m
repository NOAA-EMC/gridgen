function setup_gridgen
% -------------------------------------------------------------------------
%|                                                                        |
%|                    +----------------------------+                      |
%|                    | setup_gridgen    NOAA/NCEP |                      |
%|                    |                            |                      |
%|                    | Last Update :  07-Dec-2017 |                      |
%|                    +----------------------------+                      |
%|                     Distributed with WAVEWATCH III                     |
%|                                                                        |
%|                 Copyright 2009 National Weather Service (NWS),         |
%|  National Oceanic and Atmospheric Administration.  All rights reserved.|
%|                                                                        |
%| DESCRIPTION                                                            |
%| The setup_gridgen function supports the gridgen                        |
%| 1. downloads the reference data (etopo1 and etopo2) from NCEP's server,|
%| 2. uncompresses the tarball at the appropriate path and                |
%| 3. adds temporarily the grid_gen's paths to the MATLAB's pathdef.      |
%|                                                                        |
%| INPUT                                                                  |
%| None                                                                   | 
%|                                                                        |
%| OUTPUT                                                                 |
%| None                                                                   |
%|                                                                        |
%| NOTES                                                                  |
%| In case of updates update the section "Define paths, files and ftp     |
%| server and paths".                                                     |
%|                                                                        |
%| BUG FIXES                                                              |
%|                                                                        |
%| PRGRMR   :   Stylianos Flampouris                                      |
%| DATE     :                                                             |
%|              v.1.0 - 07-Dec-2017                                       |
% -------------------------------------------------------------------------
%%
display('grid_gen installation!')
%% Define paths, files and ftp server and paths
home = fileparts(which(mfilename)); % grid_gen directory
path_tar='reference_data';          % reference data path
path_bin='bin';                     % bin path
path_exm='examples';                % examples path
%
ftp_svr='polar.ncep.noaa.gov';      % ftp server of reference data
ftp_pth='/waves/gridgen';           % ftp path for reference data
bathy_file='gridgen_addit.tar.gz';     % reference data tarball
%% Downloading
if exist([home,'/',path_tar,'/etopo1.nc'], 'file') ~= 2
    ftp_ind=ftp(ftp_svr);
    cd(ftp_ind,ftp_pth);
    mget(ftp_ind,bathy_file,home);
    close(ftp_ind);
%% Untar the reference data
    untar([home,'/',bathy_file],path_tar);
%
    delete([home,'/',bathy_file]);
end
%% Add the bin, reference path and examples to the user's matlab path
addpath(fullfile(home, path_bin))
addpath(fullfile(home, path_tar));
addpath(fullfile(home, path_exm));

%% Install NetCDF package for octave
vers=ver;
% interpreter='Matlab';
% netcdf_install='false';
for i1=1:1:length(vers)
    if strcmpi (vers(i1).Name, 'Octave')
        pkg install -forge netcdf
    end
end
