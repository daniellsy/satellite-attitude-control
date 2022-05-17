clear;clc;
load needdata3;
n=length(fi)/6;
for i=1:1:n
  fi((i-1)*3+4:(i-1)*3+6,:)=[];
end  