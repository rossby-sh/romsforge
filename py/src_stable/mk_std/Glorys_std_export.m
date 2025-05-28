clearvars; clc; 

path = '/home/msg/DATA/GLORYS/Bio/';
list = dir([path,'*.nc']);

dep = ncread([path,list(1,1).name],'depth');
lon = ncread([path,list(1,1).name],'longitude');
lat = ncread([path,list(1,1).name],'latitude');

[xlon,ylat] = meshgrid(lon,lat);

xlon = xlon';
ylat = ylat';

[M,N] = size(xlon); L = length(dep);

CHL = single(zeros(M,N,L,length(list)));
NO3 = single(zeros(M,N,L,length(list)));
O2  = single(zeros(M,N,L,length(list)));
PH  = single(zeros(M,N,L,length(list)));
PHYC = single(zeros(M,N,L,length(list)));
PO4 = single(zeros(M,N,L,length(list)));
SI  = single(zeros(M,N,L,length(list)));

date_sta = datetime(1993,01,01,00,00,00);
date_end = datetime(2022,12,01,00,00,00);

date = date_sta : calmonths : date_end;

date_snum = find(date==datetime(2000,01,01));
date_enum = find(date==datetime(2020,12,01));

TIME = 1;
for ii = date_snum : date_enum
    chl = single(ncread([path,list(ii,1).name],'chl',[1 1 1 1],[241 181 75 1]));
    no3 = single(ncread([path,list(ii,1).name],'no3',[1 1 1 1],[241 181 75 1]));
    o2  = single(ncread([path,list(ii,1).name],'o2',[1 1 1 1],[241 181 75 1]));
    ph  = single(ncread([path,list(ii,1).name],'ph',[1 1 1 1],[241 181 75 1]));
    phyc = single(ncread([path,list(ii,1).name],'phyc',[1 1 1 1],[241 181 75 1]));
    po4 = single(ncread([path,list(ii,1).name],'po4',[1 1 1 1],[241 181 75 1]));
    si  = single(ncread([path,list(ii,1).name],'si',[1 1 1 1],[241 181 75 1]));

    CHL(:,:,:,TIME) = chl;
    NO3(:,:,:,TIME) = no3;
    O2(:,:,:,TIME) = o2;
    PH(:,:,:,TIME) = ph;
    PHYC(:,:,:,TIME) = phyc;
    PO4(:,:,:,TIME) = po4;
    SI(:,:,:,TIME) = si;   

    disp(['running time: ',num2str(TIME)]);
    TIME = TIME + 1;
end

chl_std = single(zeros(M,N,L,12));
no3_std = single(zeros(M,N,L,12));
o2_std = single(zeros(M,N,L,12));
ph_std = single(zeros(M,N,L,12));
phyc_std = single(zeros(M,N,L,12));
po4_std = single(zeros(M,N,L,12));
si_std = single(zeros(M,N,L,12));

parfor ix = 1 : M
    for iy = 1 : N
        for iz = 1 : L
          for imon = 1 : 12
              chl_std(ix,iy,iz,imon) = std(CHL(ix,iy,iz,imon:12:end));
              no3_std(ix,iy,iz,imon) = std(NO3(ix,iy,iz,imon:12:end));
              o2_std(ix,iy,iz,imon) = std(O2(ix,iy,iz,imon:12:end));
              ph_std(ix,iy,iz,imon) = std(PH(ix,iy,iz,imon:12:end));
              phyc_std(ix,iy,iz,imon) = std(PHYC(ix,iy,iz,imon:12:end));
              po4_std(ix,iy,iz,imon) = std(PO4(ix,iy,iz,imon:12:end));
              si_std(ix,iy,iz,imon) = std(SI(ix,iy,iz,imon:12:end)); 
          end
        disp(['X-pos: ',num2str(ix), '  Y-pos: ',num2str(iy) ...
            '  Z-pos: ',num2str(iz), '  T-pos: ',num2str(imon)]);
        end
    end
end
              
% clear CHL NO3 O2 PH
% parfor ix = 1 : M
%     for iy = 1 : N
%         for iz = 1 : L
%           for imon = 1 : 12
%               phyc_std(ix,iy,iz,imon) = std(PHYC(ix,iy,iz,imon:12:end));
%               po4_std(ix,iy,iz,imon) = std(PO4(ix,iy,iz,imon:12:end));
%               si_std(ix,iy,iz,imon) = std(SI(ix,iy,iz,imon:12:end));                         
%           end
%         disp(['X-pos: ',num2str(ix), '  Y-pos: ',num2str(iy) ...
%             '  Z-pos: ',num2str(iz), '  T-pos: ',num2str(imon)]);
%         end
%     end
% end
% clear PHYC PO4 SI

%% Save to NetCDF file
output_file = './glorys_21y_Bio_std.nc';

% Create dimensions
nccreate(output_file, 'longitude', 'Dimensions', {'lon', M}, 'Datatype', 'single');
nccreate(output_file, 'latitude', 'Dimensions', {'lat', N}, 'Datatype', 'single');
nccreate(output_file, 'depth', 'Dimensions', {'depth', L}, 'Datatype', 'single');
nccreate(output_file, 'time', 'Dimensions', {'time', 12}, 'Datatype', 'single');

% Write longitude, latitude, and depth
ncwrite(output_file, 'longitude', lon);
ncwrite(output_file, 'latitude', lat);
ncwrite(output_file, 'depth', dep);
ncwrite(output_file, 'time', 1:12);

% Create variables for each standard deviation
variables = {'chl_std', 'no3_std', 'o2_std', 'ph_std', 'phyc_std', 'po4_std', 'si_std'};
unit = {'mg m-3', 'mmol m-3', 'mmol m-3', '1', 'mmol m-3', 'mmol m-3', 'mmol m-3'};
var_name = {'chl', 'no3', 'o2', 'ph', 'phyc', 'po4', 'si'};
std_data = {chl_std, no3_std, o2_std, ph_std, phyc_std, po4_std, si_std};

for i = 2 : length(variables)
    nccreate(output_file, variables{i}, 'Dimensions', {'lon', M, 'lat', N, 'depth', L, 'time', 12}, 'Datatype', 'single');
    ncwrite(output_file, variables{i}, std_data{i});
    ncwriteatt(output_file, variables{i}, 'units', unit{i});
    ncwriteatt(output_file, variables{i}, 'long_name', [var_name{i}, ' standard deviation']);
end

disp('Data successfully saved to NetCDF file.');