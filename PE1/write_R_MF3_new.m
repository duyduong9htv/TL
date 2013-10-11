function write_R_MF3_new(range,svp,DEP,source_depth,receiver_depth,rmax,range_inc,frequency, dz, bathy)
  
divider=[-1 -1];
botspeed=[0.0 1700.0];
botden=[0.0 1.9];
botatt1=[280.0 0.8];
botatt2=[300.0 10.0];
cp_and_other=[1500.0 8 1 0.0];
%bathy1=[0.0 85.0];
%bathy2=[rmax 85.0];

%zmax=400;
zmax = 400;


ndz=1;
% zmplt=100;
zmplt=400;

fid=fopen('ram.in','w');
title=sprintf('Internal Wave Monte Carlo \t title');
fprintf(fid,'%s\n',title);
fprintf(fid,'%6.2f %6.3f %6.3f  freq zs zr\n',[frequency source_depth receiver_depth]);
fprintf(fid,'%6.1f %8.4f %3.0f rmax dr ndr\n',[rmax,range_inc,1]);
fprintf(fid,'%5.1f %8.4f %3.0f %5.1f     zmax dz ndz zmplt\n',[zmax dz ndz zmplt]);
fprintf(fid,'%6.1f %1i %1i %3.1f  c0 np ns rs\n', cp_and_other);

for ii=1:size(bathy,1)
    fprintf(fid,'%7.4f %7.4f \n',bathy(ii,:));
end

fprintf(fid,'%4i  %4i\n',divider) ;

fprintf(fid,'%7.1f %7.1f   z cw\n',[DEP(1) svp(1,1)]);

for ii=2:length(DEP)
	fprintf(fid,'%7.1f %7.1f\n',[DEP(ii) svp(ii,1)]);
end

fprintf(fid,'%4i  %4i\n',divider) ;
fprintf(fid,'%7.1f %7.1f   z cb\n',botspeed);
fprintf(fid,'%4i  %4i\n',divider) ;
fprintf(fid,'%7.1f %7.1f   z rhob\n',botden);
fprintf(fid,'%4i  %4i\n',divider); 
fprintf(fid,'%7.1f %7.1f   z attn\n',botatt1);
fprintf(fid,'%7.1f %7.1f\n',botatt2);
fprintf(fid,'%4i  %4i\n',divider) ;

for jj=2:length(range)
	fprintf(fid,'%7.1f  rp\n',range(jj));  
	fprintf(fid,'%7.1f %7.1f   z cw\n',[DEP(1) svp(1,jj)]);	

	for ii=2:length(DEP)
		fprintf(fid,'%7.1f %7.1f\n',[DEP(ii) svp(ii,jj)]);
	end

fprintf(fid,'%4i  %4i\n',divider) ;
fprintf(fid,'%7.1f %7.1f   z cb\n',botspeed);
fprintf(fid,'%4i  %4i\n',divider) ;
fprintf(fid,'%7.1f %7.1f   z rhob\n',botden);
fprintf(fid,'%4i  %4i\n',divider); 
fprintf(fid,'%7.1f %7.1f   z attn\n',botatt1);
fprintf(fid,'%7.1f %7.1f\n',botatt2);
fprintf(fid,'%4i  %4i\n',divider) ;
end;

fclose(fid);



