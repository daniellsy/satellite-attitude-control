clear;clc;
[filename,filepath]=uigetfile('*.01','��ѡ��rpt�����ļ�'); 
file_path_name=strcat(filepath,filename);
fid=fopen(file_path_name);
line=fgetl(fid);
while (~feof(fid))
   if strcmp(line,'Node Locations') 
      line=fgetl(fid);line=fgetl(fid); 
      while(~strcmp(line,'Node Locations') )
          [node,x,y,z]= strread(line, '%f %f %f %f %*s %*s');
          v(node,1:3)=[x y z];
          line=fgetl(fid);
      end
   end
   line=fgetl(fid);
end
save coord v;