clearvars; clc; 

path = '/home/msg/DATA/GLORYS/Phy/20year/';
path20 = '/home/msg/DATA/GLORYS/Phy/20year/2020/';

list = dir([path20,'*.nc']);   

dep = ncread([path20,list(1,1).name],'depth');
lon = ncread([path20,list(1,1).name],'longitude');
lat = ncread([path20,list(1,1).name],'latitude');

[xlon,ylat] = meshgrid(lon,lat);
xlon = xlon'; ylat = ylat';

% Define the target longitude and latitude range
longitude_range = [100, 160];
latitude_range = [10, 55];

% Find the index of the longitudes and latitudes that match the range
lon_idx = find(lon >= longitude_range(1) & lon <= longitude_range(2));
lat_idx = find(lat >= latitude_range(1) & lat <= latitude_range(2));

% Subset meshgrid of lon/lat
xlon = xlon(lon_idx, lat_idx); longitude = lon(lon_idx);
ylat = ylat(lon_idx, lat_idx); latitude  = lat(lat_idx);

[M,N] = size(xlon); % M and N are now updated to the new subset size
L = length(dep);    % Depth remains the same

T = single(zeros(M,N,L,length(list)*21));
S = single(zeros(M,N,L,length(list)*21));
U = single(zeros(M,N,L,length(list)*21));
V = single(zeros(M,N,L,length(list)*21));
Z = single(zeros(M,N,length(list)*21));

TIME  = 1;
for iyear = 2000 : 2020
    new_path = [path,num2str(iyear),'/'];
    list = dir([new_path,'*.nc']);   
   
    for ii = 1 : length(list)
        to = single(ncread([new_path, list(ii,1).name], 'thetao', [lon_idx(1), lat_idx(1), 1, 1], [length(lon_idx), length(lat_idx), 50, 1]));
        so = single(ncread([new_path, list(ii,1).name], 'so', [lon_idx(1), lat_idx(1), 1, 1], [length(lon_idx), length(lat_idx), 50, 1]));
        uo = single(ncread([new_path, list(ii,1).name], 'uo', [lon_idx(1), lat_idx(1), 1, 1], [length(lon_idx), length(lat_idx), 50, 1]));
        vo = single(ncread([new_path, list(ii,1).name], 'vo', [lon_idx(1), lat_idx(1), 1, 1], [length(lon_idx), length(lat_idx), 50, 1]));
        zo = single(ncread([new_path, list(ii,1).name], 'zos', [lon_idx(1), lat_idx(1), 1], [length(lon_idx), length(lat_idx), 1]));

        T(:,:,:,TIME) = to;
        S(:,:,:,TIME) = so;
        U(:,:,:,TIME) = uo;
        V(:,:,:,TIME) = vo;
        Z(:,:,TIME) = zo;

        disp(['running time: ',num2str(TIME)]);

        TIME = TIME + 1;
    end
end

zeta_std = single(zeros(M,N,12));

parfor ix = 1 : M
    for iy = 1 : N
          for imon = 1 : 12
              zeta_std(ix,iy,imon) = std(Z(ix,iy,imon:12:end),1);              
          end
        disp(['X-pos: ',num2str(ix), '  Y-pos: ',num2str(iy)]);
    end
end             
clear Z

temp_std = single(zeros(M,N,L,12));
salt_std = single(zeros(M,N,L,12));
u_std = single(zeros(M,N,L,12));
v_std = single(zeros(M,N,L,12));

parfor ix = 1 : M
    for iy = 1 : N
        for iz = 1 : L
          for imon = 1 : 12
              temp_std(ix,iy,iz,imon) = std(T(ix,iy,iz,imon:12:end),1);
              salt_std(ix,iy,iz,imon) = std(S(ix,iy,iz,imon:12:end),1);
              u_std(ix,iy,iz,imon) = std(U(ix,iy,iz,imon:12:end),1);
              v_std(ix,iy,iz,imon) = std(V(ix,iy,iz,imon:12:end),1);
          end
        disp(['X-pos: ',num2str(ix), '  Y-pos: ',num2str(iy) ...
            '  Z-pos: ',num2str(iz), '  T-pos: ',num2str(imon)]);
        end
    end
end             
% % clear T S 


% u_std = single(zeros(M,N,L,12));
% v_std = single(zeros(M,N,L,12));
% 
% parfor ix = 1 : M
%     for iy = 1 : N
%         for iz = 1 : L
%           for imon = 1 : 12              
%               u_std(ix,iy,iz,imon) = std(U(ix,iy,iz,imon:12:end),1);
%               v_std(ix,iy,iz,imon) = std(V(ix,iy,iz,imon:12:end),1);              
%           end
%         disp(['X-pos: ',num2str(ix), '  Y-pos: ',num2str(iy) ...
%             '  Z-pos: ',num2str(iz), '  T-pos: ',num2str(imon)]);
%         end
%     end
% end             
% clear U V

%% Save to NetCDF file
output_file = './../glorys_21y_Phy_std2.nc';

% Create dimensions
nccreate(output_file, 'longitude', 'Dimensions', {'lon', M}, 'Datatype', 'single');
nccreate(output_file, 'latitude', 'Dimensions', {'lat', N}, 'Datatype', 'single');
nccreate(output_file, 'depth', 'Dimensions', {'depth', L}, 'Datatype', 'single');
nccreate(output_file, 'time', 'Dimensions', {'time', 12}, 'Datatype', 'single');

% Write longitude, latitude, and depth
ncwrite(output_file, 'longitude', longitude);
ncwrite(output_file, 'latitude', latitude);
ncwrite(output_file, 'depth', dep);
ncwrite(output_file, 'time', 1:12);

% Create variables for each standard deviation
variables = {'temp_std', 'salt_std', 'u_std', 'v_std'};
unit = {'degrees_C', '1e-3', 'm s-1', 'm s-1'};
var_name = { 'temperature', 'salinity', 'u', 'v'};
std_data = {temp_std, salt_std, u_std, v_std};

for i = 1:length(variables)
    nccreate(output_file, variables{i}, 'Dimensions', {'lon', M, 'lat', N, 'depth', L, 'time', 12}, 'Datatype', 'single');
    ncwrite(output_file, variables{i}, std_data{i});
    ncwriteatt(output_file, variables{i}, 'units', unit{i});
    ncwriteatt(output_file, variables{i}, 'long_name', [var_name{i}, ' standard deviation']);
end

nccreate(output_file, 'zeta_std', 'Dimensions', {'lon', M, 'lat', N, 'time', 12}, 'Datatype', 'single');
ncwrite(output_file, 'zeta_std', zeta_std);
ncwriteatt(output_file, 'zeta_std', 'units', 'm');
ncwriteatt(output_file, 'zeta_std', 'long_name', ['zeta', ' standard deviation']);

disp('Data successfully saved to NetCDF file.');