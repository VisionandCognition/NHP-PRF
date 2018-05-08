function cktest(in1)
fid=fopen(in1,'w');
fprintf(fid, 'This is a test');
fclose(fid);