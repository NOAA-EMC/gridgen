% THIS IS AN EXAMPLE SCRIPT FOR GENERATING A GLOBAL GRID AND CAN BE USED 
% AS  A TEMPLATE FOR DESIGNING GRIDS

% 0. Initialization

% 0.a Path to directories 

  bin_dir = '/$HOME/gridgen/bin';             % matlab scripts location
  ref_dir = '/$HOME/gridgen/reference_data';  % reference data location
  out_dir = '/OUTPUT-PATH/';                  % output directory

% 0.b Design grid parameters

  fname_poly = 'user_polygons.flag'; % A file with switches for using user 
                                     % defined polygons. An example file 
                                     % has been provided with the reference 
                                     % data for an existing user defined 
                                     % polygon database. A switch of 0 
                                     % ignores the polygon while 1 
                                     % accounts for the polygon.

  fname = 'Global';                  % file name prefix 
  
  grid_box = [-80 0 89.5 360];       % starting and ending lat,lon for grid 
                                     % domain. Longitude data is always
                                     % defined from 0 to 360
  dx = 30.0/60.0;                    % grid resolution in x (degrees)
  dy = 30.0/60.0;                    % grid resolution in y (degrees)
  ref_grid = 'etopo1';               % reference grid source 
                                     %      etopo1 = Etopo1 grid
                                     %      etopo2 = Etopo2 grid
  boundary = 'full';                 % Option to determine which GSHHS 
                                     %  .mat file to load
                                     %         full = full resolution
                                     %         high = 0.2 km
                                     %         inter = 1 lm
                                     %         low   = 5 km
                                     %         coarse = 25 km 

% 0.c Setting the paths for subroutines

  addpath(bin_dir,'-END');

% 0.d Reading Input data

  read_boundary = 1;              % flag to determine if input boundary 
                                  % information needs to be read boundary 
                                  % data files can be significantly large 
                                  % and needs to be read only the first 
                                  % time. So when making multiple grids 
                                  % the flag can be set to 0 for subsequent 
                                  % grids. (Note : If the workspace is
                                  % cleared the boundary data will have to 
                                  % be read again)

  opt_poly = 1;                   % flag for reading the optional user 
                                  % defined polygons. Set to 0 if you do
                                  % not wish to use this option

  if (read_boundary == 1)

      fprintf(1,'.........Reading Boundaries..................\n');


      load([ref_dir,'/coastal_bound_',boundary,'.mat']); 

      N = length(bound);

      if (opt_poly == 1)
          [bound_user,Nu] = optional_bound(ref_dir,fname_poly);       
      end;                                              
      if (Nu == 0)
          opt_poly = 0;
      end;

  end;

% 1. Generate the grid

% Keep in mind that depending upon the available memory, you may not be  
% able to load the Global Reference grid on your workspace all at once. In 
% which case you will have to build the bathymetry data set in segments. 
% Both options are given below. Uncomment the one you want to use

  fprintf(1,'.........Creating Bathymetry..................\n');

  lon = [grid_box(2):dx:grid_box(4)];
  lat = [grid_box(1):dy:grid_box(3)];

% Uncomment the following lines if you want to build the bathymetry in 
% segments. We have set it up ro read in chunks of 30 x 31 degrees
% --->
%  for j = 1:6                                              
%      for k = 1:12                                        
%          lons = grid_box(2) + (k-1)*30;
%          lone = grid_box(2) + k*30;
%          lats = grid_box(1) + (j-1)*31;
%          late = grid_box(1) + j*31;
%          if late > grid_box(3)
%            late = grid_box(3);
%          end;
%          grid_boxt = [lats lons late lone]
%
%          [lont,latt,deptht] = generate_grid(ref_dir,ref_grid,grid_boxt,...
%                                  dx,dy,0.1,0,999); 
%                                                           
%          [tmp,xpos1] = min(abs(lon-lons));
%          [tmp,xpos2] = min(abs(lon-lone));
%          [tmp,ypos1] = min(abs(lat-lats));
%          [tmp,ypos2] = min(abs(lat-late));
%          depth(ypos1:ypos2,xpos1:xpos2) = deptht;
%      end;
%  end;
%-->

% Use the function below if you can load the Reference grid in your workspace 
% all at once

   [lon,lat,depth] = generate_grid(ref_dir,ref_grid,grid_box,dx,dy,0.1,0,999);
         
% 2. Computing boundaries within the domain

  fprintf(1,'.........Computing Boundaries..................\n');

% 2.a Split the domain into two parts to avoid polygons improperly 
% wrapping around

  grid_box1 = [-80 0 89.5 180.0];
  grid_box2 = [-80 180.0 89.5 360];


% 2.b Extract the boundaries from the GSHHS and the optional databases
%     the subset of polygons within the grid domain are stored in b and b_opt
%     for GSHHS and user defined polygons respectively

  [ba,N1a] = compute_boundary(grid_box1,bound);
  [bb,N1b] = compute_boundary(grid_box2,bound);
  
  N1 = 0;
  if (N1a ~= 0)
      b = ba;
      N1 = N1a;
  end;
  if (N1b > 0)
      if (N1 == 0)
          b = bb;
          N1 = N1b;
      else
          b = [ba bb];
          N1 = N1+N1b;
      end;
  end;

  if (opt_poly == 1)
      [ba_opt,N2a] = compute_boundary(grid_box1,bound_user);
      [bb_opt,N2b] = compute_boundary(grid_box2,bound_user);

      N2 = 0;
      if (N2a > 0)
          b_opt = ba_opt;
          N2 = N2a;
      end;
      if (N2b > 0)
          if (N2 == 0)
              b_opt = bb_opt;
              N2 = N2b;
          else
              b_opt = [ba_opt bb_opt];
              N2 = N2+N2b;
          end;
      end;
      clear ba_opt bb_opt;
  end;

  clear ba bb; 


% 3. Set up Land - Sea Mask

% 3.a Set up initial land sea mask. The cells can either all be set to wet
%      or to make the code more efficient the cells marked as dry in 
%      'generate_grid' can be marked as dry cells

  m = ones(size(depth)); 
  loc = depth == 999; 
  m(loc) = 0;

% 3.b Split the larger GSHHS polygons for efficient computation of the 
%     land sea mask. This is an optional step but recomended as it  
%     significantly speeds up the computational time. Rule of thumb is to
%     set the limit for splitting the polygons at least 4-5 times dx,dy

  fprintf(1,'.........Splitting Boundaries..................\n');

  b_split = split_boundary(b,10);

% 3.c Get a better estimate of the land sea mask using the polygon data sets.
%     (NOTE : This part will have to be commented out if cells above the 
%      MSL are being marked as wet, like in inundation studies)

  fprintf(1,'.........Cleaning Mask..................\n');
  
 % GSHHS Polygons. If 'split_boundary' routine is not used then replace
 % b_split with b  
 
  m2 = clean_mask(lon,lat,m,b_split,0.5);
  
 % Masking out regions defined by optional polygons
 
  if (opt_poly == 1 && N2 ~= 0)
      m3 = clean_mask(lon,lat,m2,b_opt,0.5);       
  else                                              
      m3 = m2;
  end;

% 3.d Remove lakes and other minor water bodies

  fprintf(1,'.........Separating Water Bodies..................\n');

  [m4,mask_map] = remove_lake(m3,-1,1); 
  
% 4. Generate sub - grid obstruction sets in x and y direction, based on 
%    the final land/sea mask and the coastal boundaries
    
  fprintf(1,'.........Creating Obstructions..................\n');

  [sx1,sy1] = create_obstr(lon,lat,b,m4,1,1);      
  
% 5. Output to ascii files for WAVEWATCH III
%    For global grids the first and last points in longitude are same and 
%    so the last x point is omitted

  depth = depth(:,1:end-1);                                 
  lon = lon(1:end-1);                                       
  m4 = m4(:,1:end-1);                                       
  mask_map = mask_map(:,1:end-1);
  sx1 = sx1(:,1:end-1);
  sy1 = sy1(:,1:end-1);

  depth_scale = 1000;
  obstr_scale = 100;

  d = round((depth)*depth_scale);
  write_ww3file([out_dir,'/',fname,'.depth_ascii'],d);                 

  write_ww3file([out_dir,'/',fname,'.maskorig_ascii'],m4);             

  d1 = round((sx1)*obstr_scale);
  d2 = round((sy1)*obstr_scale);
  write_ww3obstr([out_dir,'/',fname,'.obstr_lev1'],d1,d2);            

  write_ww3meta([out_dir,'/',fname],lon,lat,1/depth_scale,...
                                  1/obstr_scale);          

% 6. Vizualization (this part can be commented out if resources are limited)

  figure(1);
  clf;
  loc = find(m4 == 0);
  d2 = depth;
  d2(loc) = NaN;

  pcolor(lon,lat,d2);
  shading interp;
  colorbar;
  title(['Bathymetry for ',fname],'fontsize',14);
  set(gca,'fontsize',14);
  clear d2;

  figure(2);
  clf;
  d2 = mask_map;
  loc2 = find(mask_map == -1);
  d2(loc2) = NaN;

  pcolor(lon,lat,d2);
  shading flat;
  colorbar;
  title(['Different water bodies for ',fname],'fontsize',14);
  set(gca,'fontsize',14);
  clear d2;
   
  figure(3);
  clf;

  pcolor(lon,lat,m4);
  shading flat;
  colorbar;
  title(['Final Land-Sea Mask ',fname],'fontsize',14);
  set(gca,'fontsize',14);

  figure(4);
  clf;
  d2 = sx1;
  d2(loc) = NaN;

  pcolor(lon,lat,d2);
  shading flat;
  colorbar;
  title(['Sx obstruction for ',fname],'fontsize',14);
  set(gca,'fontsize',14);
  clear d2;

  figure(5);
  clf;
  d2 = sy1;
  d2(loc) = NaN;

  pcolor(lon,lat,d2);
  shading flat;
  colorbar;
  title(['Sy obstruction for ',fname],'fontsize',14);
  set(gca,'fontsize',14);
  clear d2;

% END OF SCRIPT
