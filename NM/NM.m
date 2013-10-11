%class definition for the Kraken normal mode TL model 

classdef NM < TL
    
    properties
	 bottomDepth; %bottom depth
	 cb; %bottom velocity
     nmesh; %number of mesh points 
     depths; %depth points where SSP measured
     SSP; %sound speed profile vector
     clow;
     chigh;     
     kr; %horizontal wavenumber vector for all modes 
     uzAll; %mode shapes at all depths 
     uz; %mode shapes at receiver depth 
     Greens; %Greens fucntion (sum of all modes)
     gr; %"Green's function for each mode
     rho; %ranges where TL will be calculated
     TL; %transmission loss calculated at all ranges 
     range; 
      
     %output from Kraken 
     xi1; 
     uz1; 
     
     srcModeShape; 
     rcvModeShape; 

	 
    end %properties 
    
    methods 
        %class constructor 
        function r = NM()
            r.zs = 65; 
            r.maxRange = 20e3; 
            r.bottomDepth = 200; 
            r.cb = 1700;
            r.nmesh = 1001; 
            r.depths = [0:0.2:200]; 
            r.dz = 0.2; 
            r.SSP = 1500*ones(size(r.depths)); %default: pekeris waveguide
            r.clow = 0; 
            r.chigh = 2000;
            r.zr = 100; 
            r.range = 2000; 
            r.rho = 1:1:r.maxRange; 
        end
        
        %set sound speed profile according to input 
        function r = setSSP(r, z, c)
            %r = getSSP(r, z, c)
            %note: z must include 0 
            r.SSP = c; 
            r.depths = z;             
        end
        
        function r = setPekeris(r)
            
        end
        
        function r = writeInputFile(r)
            fid = fopen('gen_modes1.in', 'wt'); 
            fprintf(fid, 'env_file.txt\n'); 
            fprintf(fid, 'input_freq.dat\n'); 
            fprintf(fid, 'xi1test.dat\n'); 
            fprintf(fid, 'u_z1test.dat\n'); 
            fprintf(fid, 'gvtest.dat\n'); 
            fprintf(fid, '1\n'); 
            fprintf(fid, '100\n'); 
            fclose(fid);     
        end
        
        
        function r = writeEnvFile(r)
            r.writeInputFile; %create a gen_modes1.in as input file 
            NMESH = length(r.depths); 
            fid = fopen('env_file.txt', 'wt'); 
            fprintf(fid, '''Normal mode simulation'' \t !TITLE\n');
            fprintf(fid, [num2str(r.frequency) '\t\t !FREQ (Hz) \n']);
            fprintf(fid, '1 \t\t !NMEDIA \n'); 
            fprintf(fid, '''NVW'' \t\t !OPTIONS \n'); 
            fprintf(fid, [num2str(NMESH) ' 0.0 ' num2str(r.bottomDepth) ' \t \t !NMESH SIGMA(m) Z (bottom depth)\n']);
            fprintf(fid, '%10.3f\t\t%10.3f\t\t0.0 1.0 0.00006 0.0\n', r.depths(1), r.SSP(1));
            
            for k = 2:length(r.SSP)
                fprintf(fid, '%10.3f\t\t%10.3f/ \n'  ,r.depths(k) , r.SSP(k));
            end
            
            fprintf(fid, '''A'' 0.0 !BOTOPT SIGMA (m) \n'); 
            fprintf(fid, '%10.3f %10.3f\t 0.0 1.9 0.8 0.0\n', r.depths(k), r.cb); 
            fprintf(fid, '%10.3f %10.3f !CMIN CMAX \n', r.clow, r.chigh); 
            fprintf(fid, '%10.1f \t \t !RMAX (m)\n', r.maxRange); 
            fprintf(fid, '1 \n'); 
            fprintf(fid, '1e-3 /                ! NSD  SD(1:NSD)\n'); 
            noReceivers = r.bottomDepth/r.dz; 
            fprintf(fid, '%5.0f\n', noReceivers); 
            fprintf(fid, '%7.1f\t %7.1f  /                ! NRD  RD(1:NRD)\n', r.dz, r.bottomDepth); 
            fclose(fid); 
            
            %write frequency input file
            fid = fopen('input_freq.dat', 'wt');
            fprintf(fid, '%10.3f', r.frequency); 
            fclose(fid); 
        end
        
        function r = getModalInfo(r)
            r.writeEnvFile; 
            
            %execute fortran code to get the modal info 
            !/home/dtran/MatlabCodes/NM/modal_info_vs_depth.out < gen_modes1.in
%  >suppressed_log
            constant=4*pi*exp(-j*pi/4)/sqrt(8*pi);
            
            %reading results 
            load xi1test.dat %horizontal wavenumber for each mode, separated in real and imaginary parts
            xi1=xi1test(:,1)+j*xi1test(:,2);
            nmax=length(xi1);           
                    
            r.xi1 = xi1; 
            r.kr = xi1; 
               %maximum number of important modes

            load u_z1test.dat; %mode shapes 
            u=u_z1test(:,1)+j*u_z1test(:,2);
            uz=transpose(reshape(u,length(u)/nmax,nmax));  %separates u into its modes
            uz(:,1)=[];   
            
            rcv_ind = round(r.zr/r.dz); 
            r.uzAll = uz; 
            
            src_ind = round(r.zs/r.dz); 
            r.srcModeShape = uz(:, src_ind); 
            r.rcvModeShape = uz(:, rcv_ind);
            
            
            uz = uz(:, rcv_ind); 
            r.uz1 = uz; 
            
            
            r.gr = constant*r.rcvModeShape.*exp(j*xi1*r.range).*r.srcModeShape./sqrt(xi1*r.range); 
            
            %padding 
            
            if length(xi1) >100
                xi1 = xi1(1:100); 
            end
                                    
            %pad in zeros to have maximum of 100 modes
            while length(xi1) < 100
                xi1 = [xi1; 0]; 
            end    
            
            if size(uz, 1) >100
                uz = uz(1:100, :); 
            end
            
            while size(uz, 1) <100
                uz = [uz; zeros(ones, size(uz, 2))]; 
            end
            
            if size(r.gr, 1) >100
                r.gr = r.gr(1:100, :); 
            end
            
            while size(r.gr, 1) <100
                r.gr = [r.gr; zeros(ones, size(r.gr, 2))]; 
            end
            
            r.uz = uz; 
            r.kr = xi1; 
            r.Greens = sum(r.gr); 
            
%             r.calculateTL; 
               
        end
        
        
        %calculate the complex Green's function at all ranges specified by
        %r.rho 
        function r = calculateTL(r)
            
%             load xi1test.dat %horizontal wavenumber for each mode, separated in real and imaginary parts
%             xi1=xi1test(:,1)+j*xi1test(:,2);
%             nmax=length(xi1);           
%                     
%      
%             load u_z1test.dat; %mode shapes 
%             u=u_z1test(:,1)+j*u_z1test(:,2);
%             uz=transpose(reshape(u,length(u)/nmax,nmax));  %separates u into its modes
%             uz(:,1)=[];   
%             
%             rcv_ind = round(r.zr/r.depth_inc); 
%             uz = uz(:, rcv_ind); 
            
            r.TL = zeros(size(r.rho));
            constant=((4*pi)*exp(-j*(pi/4)))/(sqrt(8*pi));
            
%             for k = 1:length(r.rho)                
%                 gr = constant*uz.*(exp(j*xi1.*r.rho(k)))./(sqrt(xi1.*r.rho(k))); 
%                 r.TL(k) = sum(gr); 
%             end 
            
            for k = 1:length(r.rho)        
                gr= constant*r.rcvModeShape.*exp(j*r.xi1*r.rho(k)).*r.srcModeShape./sqrt(r.xi1*r.rho(k)); 
                r.TL(k) = sum(gr); 
            end 
            
        end
        
        
        %calculate a 2-D Green's function grid at all ranges and depths 
        function r = getTLgrid(r)
%             r.getModalInfo(); 
            constant=4*pi*exp(-j*pi/4)/sqrt(8*pi);
            r.gGrid = zeros(length(r.depths), length(r.rho)); 
            
            %reading results 
            load xi1test.dat %horizontal wavenumber for each mode, separated in real and imaginary parts
            xi1=xi1test(:,1)+j*xi1test(:,2);
            nmax=length(xi1);           
                    
            r.xi1 = xi1; 
            r.kr = xi1; 
               %maximum number of important modes

            load u_z1test.dat; %mode shapes 
            u=u_z1test(:,1)+j*u_z1test(:,2);
            uz=transpose(reshape(u,length(u)/nmax,nmax));  %separates u into its modes
            uz(:,1)=[];   
            
            rcv_ind = round(r.zr/r.dz); 
            r.uzAll = uz; 
                       
            
            src_ind = round(r.zs/r.dz); 
            r.srcModeShape = uz(:, src_ind);
            
            
            r.gGrid = 0; 
            
            for K = 1:length(r.xi1)
                disp('calculating Green''s function grid'); 
                disp(K/length(r.xi1))
                kr = r.xi1(K); 
                denom = repmat(sqrt(kr*r.rho), length(r.depths)-1, 1); 
                gr = constant*r.srcModeShape(K)*transpose(uz(K, :))*exp(j*kr*r.rho)./denom;
                r.gGrid = r.gGrid + gr; 
            end
            
            
            
            
                   
        end
        
        
    end 

end

