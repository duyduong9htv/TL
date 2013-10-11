%read_ggrid.m simply reads g_grid.dat and puts g_grid into workplace
%environment

%open g_grid:
% fid = fopen('g_grid.dat');
% junk = fread(fid, 8, 'int');
% nz = junk(3)+4;
% g_grid = fread(fid, [2 * nz, range_max/range_inc], 'float32');
% fclose(fid);
% g_grid = g_grid(1:nz,:) + i*g_grid(nz+1:end,:);
% g_grid = g_grid(1:nz-5,:);              %Greens function for a single frequency
% clear junk nz fid

fid=fopen('g_grid.dat');
junk=fread(fid,3,'int');
nz=junk(2)+2;
%g_grid=fread(fid, [2*nz, nr], 'float32');
g_grid=fread(fid, [2*nz, inf], 'float32');
g_grid=g_grid(1:nz,:)+i*g_grid(nz+1:end,:);%this addition step is the most time-consuming one 
%g_grid=g_grid(1:nz-5,:);
g_grid=g_grid(6:nz,1:end-1);
fclose(fid);

%% testing if this gets pushed to github


%  figure(1); imagesc(10*log10(abs(g_grid)))