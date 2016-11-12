fiopen ='wget_ruD6ONxY.txt';

fid=fopen(fiopen);
ss=fscanf(fid,'%s\n');
fclose(fid);

for i=1:31
   urls{i}=ss((i-1)*438+(1:438));
end

for i=1:31
   urlwrite(urls{i},['test_' num2str(i) '.nc'])
   i
end

% nc_varget from mexcdf.sourceforge.net
u=nc_varget('test_2.nc','u');
