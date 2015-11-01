fid=fopen('einx.txt','wt');
fprintf(fid,'%g\n',Etranx);
fclose(fid);
fid=fopen('einy.txt','wt');
fprintf(fid,'%g\n',Etrany);
fclose(fid);
fid=fopen('einz.txt','wt');
fprintf(fid,'%g\n',Etranz);
fclose(fid);

fid=fopen('etx2.txt','wt');
fprintf(fid,'%g\n',Etranx);
fclose(fid);
fid=fopen('ety2.txt','wt');
fprintf(fid,'%g\n',Etrany);     
fclose(fid);
fid=fopen('etz2.txt','wt');
fprintf(fid,'%g\n',Etranz);
fclose(fid);

