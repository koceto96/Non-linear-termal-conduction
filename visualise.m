function visualise

fid = fopen('solution.txt');

m = fscanf(fid,'%d',1);
u = fscanf(fid,'%g',[m+1,m+1]);

fclose(fid);

x = (1/m)*[0:m]';

warning off MATLAB:hg:AutoSoftwareOpenGL

mesh(x,x,u);
view([150 30]);
colorbar;

print(gcf,'-depsc2', 'solution.ps');
system('gs -q -sPAPERSIZE=a4 -sDEVICE=pdfwrite -o solution.pdf -c "<</PageOffset [40 200]>> setpagedevice" -f solution.ps');
delete('solution.ps');

fprintf('Success! solution.pdf was created\n');
exit;


