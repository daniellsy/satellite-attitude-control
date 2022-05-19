clc;clear;
[filename,filepath]=uigetfile('*.f06','请选择f06数据文件'); 
file_path_name=strcat(filepath,filename);
fid=fopen(file_path_name);
i=1;j=1;
line=fgetl(fid);
while (~feof(fid))
 if length(line)>=74
    if strcmp(line(54:74),'L O A D   V E C T O R')
      line=fgetl(fid); line=fgetl(fid); line=fgetl(fid);       
      while(line(1)~='1'&&line(2)~='*')
         [pointID,~,T1,T2,T3,R1,R2,R3]=strread(line,'%d %s %f %f %f %f %f %f');
         if abs(T1)>0.00001
         m(pointID,:)=abs(T1);
         else if abs(T2)>0.00001
                 m(pointID,:)=abs(T2);
             else m(pointID,:)=abs(T3);
             end
         end
         j=j+1;
         line=fgetl(fid);
      end 
       i=i+1; 
       j=1;
    end
 end
 line=fgetl(fid);
end
M=sum(m);
save M m;