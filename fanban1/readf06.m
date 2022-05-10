clc;clear;
[filename,filepath]=uigetfile('*.f06','请选择f06数据文件'); 
file_path_name=strcat(filepath,filename);
fid=fopen(file_path_name);
i=1;j=1;X=zeros(10,2);
line=fgetl(fid);
while (~feof(fid))
 if length(line)>80
    if strcmp(line(42:80),'R E A L   E I G E N V E C T O R   N O .')
     [~, CYCLE,ORDER]= strread(line, '%s=%f %*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s %f');       
       X(i,1)=ORDER;
       X(i,2)=CYCLE;
      line=fgetl(fid); line=fgetl(fid); line=fgetl(fid);       
      while(line(1)~='1')
         [pointID,~,T1,T2,T3,R1,R2,R3]=strread(line,'%d %s %f %f %f %f %f %f');
         Modes(i).mode(j,:)=[T1,T2,T3,R1,R2,R3];
         j=j+1;
         line=fgetl(fid);
      end 
       i=i+1; 
       j=1;
       
    end
 end
 line=fgetl(fid);
end

k=1;
while X(k,1)==X(k+1,1)
    k=k+1;
end
K=k;
n=1;m=length(X(:,1))/K;
ORDER=zeros(m,2);
while n<=m
   ORDER(n,:)=X(n*K,:); 
    M(n).m=Modes(n*K).mode;
   while k>1
      k=k-1;
      M(n).m=[ Modes((n-1)*K+k).mode;M(n).m];
   end
   k=K;
   n=n+1;
end
y1=ORDER;
y2=M;
y3=pointID;
%振型矩阵
nksi=size(ORDER,1);%模态数
mode=struct2cell(y2);
node=length(cell2mat(mode(:,:,1)));
fi=zeros(node*6,nksi);
for i=1:1:nksi      
    zhenxing=cell2mat(mode(:,:,i));  
    k=1;
    for j=1:6:node*6-5
fi(j:j+5,i)=zhenxing(k,:);
k=k+1;
    end
end
save needdata3 fi ORDER;