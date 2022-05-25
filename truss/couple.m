clear;clc;
load fi.mat;
load M.mat;
load coord.mat;%附件节点坐标、质量、振型
% fix=[(3.245+2.452)/2 -0.13 (0.42-0.372)/2];%连接点坐标
mp=[-0.4958 0.2 0.2];%整体质心坐标
[p,q]=size(fi);num=p/3;
for i=1:1:num
r(i,:)=v(i,:)-mp;
end
H=zeros(3,q);
 for i=1:1:num
     for j=1:1:q
t1=r(i,:);t2=fi([3*i-2:3*i],j)';
t(:,j)=cross(t1,t2)'*m(i);
end
 H=H+t;
 end
 save couple.mat H;
% w1=t1(1);w2=t1(2);w3=t1(3);
% b1=t2(1);b2=t2(2);b3=t2(3);
% t2=[w2*b3-w3*b2 w1*b3-w3*b1 w2*b1-w1*b2]; 