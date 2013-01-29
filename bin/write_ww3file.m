function [messg,errno] = write_ww3file(fname,d)

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
%| Write the output array into ascii file                                 |
%|                                                                        |
%| [messg,errno] = write_ww3file(fname,d)                                 |
%|                                                                        |
%| INPUT                                                                  |
%|  fname       : Output file name                                        |
%|  d           : Output 2D array                                         | 
%|                                                                        |
%| OUTPUT                                                                 |
%|  messg       : Error message. Is blank if no error occurs              |
%|  errno       : Error number. Is zero for succesful write               |
% -------------------------------------------------------------------------

[Ny,Nx] = size(d);

fid = fopen(fname,'w');

[messg,errno] = ferror(fid);

if (errno == 0)
   for i = 1:Ny
       a = d(i,:);
       fprintf(fid,' %d ',a);
       fprintf(fid,'\n');
   end;
else
   fprintf(1,'!!ERROR!!: %s \n',messg);
end;

fclose(fid);

return;
