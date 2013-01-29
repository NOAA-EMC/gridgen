function [b,usr_cnt] = optional_bound(ref_dir,fname)

% -------------------------------------------------------------------------
%|                                                                        |
%|                    +----------------------------+                      |
%|                    | GRIDGEN          NOAA/NCEP |                      |
%|                    |                            |                      |
%|                    | Last Update :  23-Oct-2012 |                      |
%|                    +----------------------------+                      | 
%|                     Distributed with WAVEWATCH III                     |
%|                                                                        |
%|                 Copyright 2009 National Weather Service (NWS),         |
%|  National Oceanic and Atmospheric Administration.  All rights reserved.|
%|                                                                        |
%| DESCRIPTION                                                            |
%| This routine reads an optional polygon data set that has been created  |
%| to mask out water bodies that do not play a major role in wave         |
%| propagation. The polygons that  are included are based on switches read|
%| from a file                                                            |
%|                                                                        |
%| [b,usr_cnt] = optional_bound(ref_dir,fname,icoords)                    |
%|                                                                        |
%| INPUT                                                                  |
%|   ref_dir : PATH to reference directory that includes the file         |
%|             "optional_coastal_polygons.mat"                            |
%|   fname   : Filename that has a list of switches to choose which of the|
%|             optional polygons need to be switched off or on. The       |
%|             switches in the file should coincide with the polygons in  | 
%|             the "optional_coastal_polygons.mat" file. Example files for|
%|             the two are provided in reference data directory           |
%|                                                                        |
%| OUTPUT                                                                 |
%|   b       : An array of boundary polygon data structures               | 
%|   usr_cnt : Total number of polygons found                             |
%|                                                                        |
% -------------------------------------------------------------------------

fid = fopen(fname,'r');

load([ref_dir,'/optional_coastal_polygons.mat']);

N = length(user_bound);
usr_cnt = 0;

for i = 1:N
    a = fgetl(fid);    
    a1= sscanf(a,'%d%d');
    if (a1(2) == 1)
        b(usr_cnt+1) = user_bound(i);
        usr_cnt = usr_cnt+1;
        b(usr_cnt).west = min(b(usr_cnt).x);
        b(usr_cnt).east = max(b(usr_cnt).x);
    end;
end;
fclose(fid);

if (usr_cnt == 0)
  b(1) = -1;
end;

return
