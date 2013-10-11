classdef TL < handle 
    %class description of the TL (Transmission Loss) goes here 
    
    properties (Access = public) %public access for inheritance purposes 
        zs; %source depth 
        maxRange; %max range 
        zmax; %maximum depth for RAM output
        dr; %range increment for RAM
        dz; %depth increment for RAM (should be 1/10 wavelength)
        frequency; %frequency for each run
        ssps_bank; %all the random SSPs on the bank
        ssps_basin; %SSPs in basin 
        zr;  %receiver depth
        gGrid; %Green's function grid, for one single run of single frequency
        ssps; %sound speed profiles at each updated range interval   
        depth; %depth vector, eg. 0:0.2:200
        ranges; 
        bathy; %bathymetry (Nx2) matrix, consisting of range and depth           
        
    end
    
    methods 
        
        %% class constructor 
        function r = TL() %default class constructor 
                %default source depth,max PE range and max PE depth for a run 
                r.zs = 65; 
                r.maxRange = 20e3; 
                r.ranges = [0:500:r.maxRange]; 
                load ssps_basin
                r.ssps_basin = ssps_basin; 
                load ssps_bank
                r.ssps_bank = ssps_bank; 

                r.depth = [0.2:0.2:200];
                r.bathy = [0:50:r.maxRange; 
                        200*ones(size(0:50:r.maxRange))]'; 
                r.frequency = 415; %default frequency 
                r.zr = 105; 
                r.ssps = []; 
                for k = 1:length(r.ranges)
                    ssp = 1500*ones(size(r.depth))'; 
                    r.ssps = [r.ssps ssp]; 
                end
                r.dr = 1; 
                r.dz = 0.2; 

        end  
        
        
         %% randomize the sound speed profile
        function r = randomSSP(r)
            shidx = find(r.bathy(:,2)<=80); %shallow parts of waveguide (bank)
            bathy_data = r.bathy; 
            ssps_basin = r.ssps_basin; 
            ssps_bank = r.ssps_bank; 
            rough_range = r.ranges; 
            corr_length = 500;
            if ~isempty(shidx)
                basin = floor(bathy_data(shidx(1),1)/corr_length);
                rand_ssps = ceil((size(ssps_basin,2)-1).*rand(basin,1))+1;
                rand_ssps = [rand_ssps;ceil((size(ssps_bank,2)-1).*rand(length(rough_range)-basin,1))+1];
                for ss = 1:basin
                    rough_svp(:,ss) = ssps_basin(:,rand_ssps(ss));
                end
                for ss = basin+1:length(rough_range)
                    rough_svp(:,ss) = ssps_bank(:,rand_ssps(ss));
                end
            else
                rand_ssps=ceil((size(ssps_basin,2)-1).*rand(length(rough_range),1))+1; 
                for ss=1:length(rough_range)
                    rough_svp(:,ss)=ssps_basin(:,rand_ssps(ss));
                end
            end
            
            r.ssps = rough_svp; 
        end

    end %methods 
    
    
end 