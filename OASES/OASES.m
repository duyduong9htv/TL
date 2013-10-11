classdef OASES < TL %TL class to run OASES is a child class of the TL class using the PE model 
% This class handles the interfacing to the fortran wavenumber integration
% model OASES, written by Prof. Henrik Schmidt at MIT. 
% Documentation goes here: 
% 
% 
    properties
        pathToOASES; %path to OASES binary executable file 
        cb; %bottom sound speed 
        NW; %number of wavenumber to integrate over (should be power of 2) 
    end %properties 
    
    
    methods
        
        function r = OASES()
            r.ssps = r.ssps(:, 1); 
            r.maxRange = 1000; %by default 1 km = max range, because the wave number integration is very 
            %time consuming at larger ranges. 
            r.cb = 1700; 
            r.NW = -1;
            r.zr = [0:2:300]; 
        end
        
%         %set path to OASES binary files 
%         function r = setPathToOASES(r, dir_path)
%             r.pathToOASES = dir_path; 
%         end

        
        %set the waveguide to be isovelocity with sound speed c 
        function r = setIsovelocity(r, c)
            %the whole water column is considered 1 single layer now 
            r.depth = [0 r.depth(end)]; 
            r.ssps = c*ones(2, 1); 
        end

        
        %% Input file writing functions 
        
        function r = writeEnvOASP(r)
            fid = fopen('input.dat', 'wt'); 
            fprintf(fid, 'Wave number integration model OASP module input file\n');
            fprintf(fid, 'N J F \n'); 
            fprintf(fid, '%4.0f  0\n', r.frequency); 
            fprintf(fid, '%4.0f\n', length(r.depth)+1); %add one layer for vaccuum (air) half space   
            %vacuum layer on top, all values are zeroes  
            fprintf(fid, '%4.0f %6.0f %4.0f %3.1f %3.1f %3.1f %1.0f\n', 0, 0, 0, 0, 0, 0, 0); 
            
            %write depth-dependent sound speed profile 
            for k = 1:(length(r.depth)-1)
                fprintf(fid, '%4.1f %6.2f %4.0f %3.1f %3.1f %3.1f %1.0f\n',...
                            r.depth(k), r.ssps(k), 0, 0.0, 0, 1, 0); 
            end 
            
            %write bottom parameters 
            fprintf(fid, '%4.1f %6.2f 400 0.8 0.5 1.9 0\n', r.depth(end), r.cb); 
            
            %write source depth (single point source as of now) 
            fprintf(fid, '%4.0f\n', r.zs); 
            
            %write receiver depths (at all depth sample points) 
            fprintf(fid, '%6.2f %6.2f %5.0f\n', 0, max(r.zr), length(r.zr)); 
            
            %write min/max sound speed allowed 
            fprintf(fid, '%6.0f %6.0f\n', 750, 1e8); 
           
            fprintf(fid, '%7.0f 1 %7.0f 1\n', r.NW, r.NW*0.9); 
            fprintf(fid, '%5.0f %5.1f %5.1f %10.5f %10.0f %10.4f %10.0f\n',...
                          r.NW, r.frequency, r.frequency, 0.0001, 0, r.dr/1000, r.maxRange/r.dr); 
            fclose(fid); 
        end
        
        
        %% calculates the complex field using the wavenumber integration model 
        %OASP module
        function r = calculateField(r)
            export_path_str = ['!export PATH=' r.pathToOASES ':$PATH']; 
            eval(export_path_str); 
            !/home/dtran/OASES/bin/oasp input 
            
            %%reading out results 
            [out,sd,z,range,f,fc,omegim]=trf_reader_oases('input.trf');
            
            
            r.gGrid = out;           
            
            
        end
        
        
        
        
        %% calculate the transmission loss (for fixed receiver depths) 
        %OAST module. This does not output the complex field, but a wavenumber spectrum is
        %obtained 
        function r = calculateTL(r)
            
        end
        
        
        
        %% plotting results 
        
        function r = plotTL(r)
            
        end
        
    end %methods 
    
     
    
end % end classdef 
