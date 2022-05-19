clear;clc;
load needdata.mat;
nksi=size(ORDER,1);%模态数
for i=1:nksi
    omega(i,i)=2*pi*ORDER(i,2);               %omega对角阵
    omega2(i,i)=(2*pi*ORDER(i,2))^2;          %omega平方对角阵
    zeta(i,i)=0.01;                      %z阻尼比矩阵系数
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%状态空间参数求值
% omega(1:6,:)=0;
% omega2(1:6,:)=0;    %去除刚体模态频率
%%%%%%%%%%%%%%%%%广义MCK系数阵
Mg=eye(nksi);       %广义质量阵
Cg=2*zeta*omega;    %广义阻尼阵
Kg=omega2;          %广义刚度阵
node=length(fi);
%%%%%%%%%%%%%%%%%连续系统ABCD系数阵
A=[zeros(nksi)  eye(nksi); -Kg -Cg];
% B3=[zeros(nksi,node);fi'];
% C=blkdiag(eye(nksi),eye(nksi));
% D=zeros(2*nksi,nksi);
% innum=5676;%施加外力节点编号
% inn=3;%施加力或力矩的自由度
save para.mat A Kg Cg fi;