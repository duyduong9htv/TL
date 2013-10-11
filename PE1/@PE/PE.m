classdef PE < TL 
    %class description of the Parabolic Equation (PE) model class 
    
    properties %additional properties other than the publicly inherited from class TL. Only functions of this class can use these properties, other siblings can't have access.         
          
        fStart; %fstart
        fStop; %fstop, for broadband PE runs
        df; %frequency step               
        PSgrid; %total energy grid (Parseval sum), for broadband signals
   
        
        greenArray;         
        fs; %sampling frequency; 
        Td; %duration of signal (in case want to calculate received field in time and freq domain)
        NoHydrophone; 
        incidentAng; %incident angle on the array (for AI applications) 
        range; %real range between source and receiver (for AI broadband calculation)
        spacing; %hydrophone spacing; 
        x1; %source UTM X
        y1; %source UTM Y
        x2; %receiver UTM X
        y2; %receiver UTM Y
    end %properties 
    
    
    %% METHODS 
    methods
        %% class constructor 
        function r = PE() %for each run
            %default source depth,max PE range and max PE depth for a run 
            r.zs = 65; 
            r.NoHydrophone = 64; %default number of hydrophones 
            r.maxRange = 20e3; 
            r.zmax = 400; 
            load ssps_basin
            r.ssps_basin = ssps_basin; 
            load ssps_bank
            r.ssps_bank = ssps_bank; 
            r.ranges = [0:500:r.maxRange]; 
            r.depth = [0.2:0.2:200];
            r.bathy = [0:50:r.maxRange; 
                        200*ones(size(0:50:r.maxRange))]'; %default bathy, flat 200 m 
            r.frequency = 415; %default frequency 
            r.fs = 8000; %default sampling rate 
            r.spacing = 1.5; 
            r.zr = 105; 
            r.ssps = []; 
            for k = 1:length(r.ranges)
                ssp = 1500*ones(size(r.depth))'; 
                r.ssps = [r.ssps ssp]; 
            end
            r.dr = 1; 
            r.dz = 0.2; 
            
        end        
        
       
        %% get a transect for bathymetry 
        
        function r = getTransectUTM(r)
               x1 = r.x1; x2 = r.x2; 
               y1 = r.y1; y2 = r.y2; 
               dr = 50;
              fid=fopen('bathymetry.arr','rb','ieee-le');
              Xsize = fread(fid,1,'int32');
              Ysize = fread(fid,1,'int32');
              Imagedata = zeros(Xsize,Ysize);

              Imagedata=fread(fid,[Xsize,Ysize],'float32');
              fclose(fid);

              Imagedata2=rot90(Imagedata); 
              clear Imagedata;

              bathymetry;

              datax=[grid_xmin:grid_inc:grid_xmax];
              datay=[grid_ymax:-grid_inc:grid_ymin];

              len = ddist([x1 y1], [x2 y2]); 

              rtemp=[0:dr:round(len)];

              phii=atan2((y2-y1),(x2-x1));

              xtemp=x1+rtemp*cos(phii);
              ytemp=y1+rtemp*sin(phii);

              bathytemp=interp2(datax,datay,Imagedata2,xtemp,ytemp);

              data=[rtemp.' bathytemp.'];

              data=[data(:,1) abs(data(:,2))];  
              r.bathy = data; 
 
 
        end 
        
        function r = getTransectUTM2(r)
%                [lat1, lon1] = utm2deg(r.x1, r.y1, '19 T'); 
%                [lat2, lon2] = utm2deg(r.x2, r.y2, '19 T'); 
               
               load bathy_large; % X(lon), Y(lat), Z(depth)
               datax = zeros(size(X)); 
               datay = zeros(size(Y)); 
               for ii = 1:length(X)
                   [xx, yy, zone] = deg2utm(Y(1), X(ii));
                   if ~strcmp(zone, '19 T')
                    disp(zone);
                   end
                   
                    datax(ii) = xx;  
               end
               
               for ii = 1:length(Y)
                   [xx, yy, ~] = deg2utm(Y(ii), X(1));
                    datay(ii) = yy;  
               end
               
%                datay = flipud(datay); 
               
               len = ddist([r.x1 r.y1], [r.x2 r.y2]); 
               dr = 50;
%               fid=fopen('bathymetry.arr','rb','ieee-le');
%               Xsize = fread(fid,1,'int32');
%               Ysize = fread(fid,1,'int32');
%               Imagedata = zeros(Xsize,Ysize);
% 
%               Imagedata=fread(fid,[Xsize,Ysize],'float32');
%               fclose(fid);
% 
%               Imagedata2=rot90(Imagedata); 
%               clear Imagedata;
% 
%               bathymetry;
% 
%               datax=[grid_xmin:grid_inc:grid_xmax];
%               datay=[grid_ymax:-grid_inc:grid_ymin];

              
              rtemp=[0:dr:round(len)];

              phii=atan2((r.y2-r.y1),(r.x2-r.x1));

              xtemp=r.x1+rtemp*cos(phii);
              ytemp=r.y1+rtemp*sin(phii);

              bathytemp=interp2(datax,datay,Z',xtemp,ytemp);

              data=[rtemp.' bathytemp.'];

              data=[data(:,1) abs(data(:,2))];  
              r.bathy = data; 
 
 
        end 
        
        
        
       %% 
        
        function r = randomSSP2013(r)
            load ssp2013.mat; 
            rough_range = r.ranges; 
            for ss = 1:length(rough_range)
                k = floor(size(ssp, 2)*rand + 1); 
                rough_svp(:, ss) = ssp(:, k); 
            end
            r.ssps = rough_svp; 
        end
        
        
        %% calculate Green's function: run PE model 
        function r = calculateGreen(r)
            disp('Writing RAM input...'); 
            write_R_MF3_new(r.ranges,r.ssps,r.depth,...
                             r.zs,100,...
                                r.maxRange,r.dr,r.frequency,r.dz,r.bathy);% write ram.in
            disp(['Done!' char(10) 'Running RAM...'])
            %replace the path below with path to your executable RAM file 
            !/home/dtran/PEmodel/ramDD > temp 
            disp(['Done!' char(10) 'Reading data out...']); 
            read_ggrid;
            r.gGrid = g_grid; 
            clear g_grid;
            disp('Done!'); 
        end
        
        %% calculate Green's function for a source array 
        
        function r = calculateGreenSourceArray(r)
            disp('writing RAM input...'); 
            write_R_MF3_new(r.ranges,r.ssps,r.depth,...
                             r.zs,100,...
                                r.maxRange,r.dr,r.frequency,r.dz,r.bathy);% write ram.in
            disp(['Done!' char(10) 'Running RAM...'])
            
            if r.frequency < 800
                !scp /home/dtran/MatlabCodes/TLbox/source_spacing.in ./source_spacing.in 
                
            else 
				!scp /home/dtran/test/Source/source_spacing_mod30.in ./source_spacing.in             
            end 

            %replace the path below with path to your executable RAM file 
            !/home/dtran/PEmodel/ram1.5.4_src_array.out     
            disp(['Done!' char(10) 'Reading data out...']); 
            read_ggrid;
            r.gGrid = g_grid; 
            clear g_grid;
            disp('Done!'); 
        end
        
        
        %% calculate Green's function, Pekeris waveguide
        
        function r = calculateGreenPekeris(r)
            
        end
        
        %% calculate broadband Transmission loss 
        
        function r = calculatePS(r)
            r.PSgrid = 0; 
            count = 0; 
            for f = r.fStart:r.df:r.fStop
                count = count + 1; 
                r1 = r; 
                r1.frequency = f; 
                r1.calculateGreen; 
                r.PSgrid = r.PSgrid + abs(r1.gGrid).^2;                 
            end 
            r.PSgrid = r.PSgrid/count; 
        end
        
        %% calculate Green's function to an array 
        %this is to investigate the beamforming capability and range
        %localization using the array invariant
        function r = calculateTLarray(r)
            r.greenArray = zeros(r.fs*r.Td, r.NoHydrophone); 
            theta = r.incidentAng; 
            range_step = r.spacing*sind(theta); 
            depth_ind = round(r.zr/r.dz);
            if range_step ~=0 
                range_h = (r.range - (r.NoHydrophone-1)/2*range_step):...
                            range_step:(r.range + (r.NoHydrophone-1)/2*range_step); %ranges to 64 hydrophones 
            else 
                range_h = ones(r.NoHydrophone, 1)*r.range; %case of broadside arrival (0 incident angle) 
            end
             
            %sweeping through all frequencies
            
            for f = r.fStart:r.df:r.fStop
                disp((f-r.fStart)/(r.fStop - r.fStart)); 
                r1 = r; 
                r1.frequency = f; 
                r1.maxRange = max(range_h)+r1.spacing + 1; %just to cover up to the farthest hydrophone 
                r1.calculateGreen; 
                f_ind = round(f*r.Td); 
                for hy = 1:r.NoHydrophone
                    range_ind = round(range_h(hy)/r.dr); 
                    r.greenArray(f_ind, hy) = r1.gGrid(depth_ind, range_ind); 
                end
            end 
                
            %time shift for corresponding incident angle            
            delay_distance = range_h - r.range;  
            delay_matrix = exp(j*2*pi*[linspace(0, r.fs, size(r.greenArray,1))]'*delay_distance/1500);
            r.greenArray = r.greenArray.*delay_matrix;
        
        end
        
        
        %% PLOTTING METHODS 
        
        %plot Transmission Loss        
        function r = plotTL(r)
            range = 0:r.dr:r.dr*size(r.gGrid, 2); 
            depth = 0:r.dz:(r.dz*size(r.gGrid, 1)); 
            imagesc(range, depth, 20*log10(abs(r.gGrid))); 
            caxis([-100 -50]); 
            setFont(18); setFigureAuto; 
            xlabel('Range (m)'); 
            ylabel('Depth (m)'); 
            colorbar; 
            setFont(18); setFigureAuto; 
        end
        
        %plot broadband transmission loss (Parseval sum)
        function r = plotPS(r)
            range = 0:r.dr:size(r.PSgrid, 2); 
            depth = 0:r.dz:(dz*size(r.gGrid, 1)); 
            imagesc(range, depth, 10*log10(abs(r.PSgrid))); 
            caxis([-100 -50]); 
            setFont(18); setFigureAuto; 
            xlabel('Range (m)'); 
            ylabel('Depth (m)'); 
            colorbar; 
            setFont(18); setFigureAuto; 
        end
        
        
        
        
    end %methods 
end %classdef
