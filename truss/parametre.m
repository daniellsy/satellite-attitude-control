clear;clc;
load needdata.mat;
nksi=size(ORDER,1);%ģ̬��
for i=1:nksi
    omega(i,i)=2*pi*ORDER(i,2);               %omega�Խ���
    omega2(i,i)=(2*pi*ORDER(i,2))^2;          %omegaƽ���Խ���
    zeta(i,i)=0.01;                      %z����Ⱦ���ϵ��
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%״̬�ռ������ֵ
% omega(1:6,:)=0;
% omega2(1:6,:)=0;    %ȥ������ģ̬Ƶ��
%%%%%%%%%%%%%%%%%����MCKϵ����
Mg=eye(nksi);       %����������
Cg=2*zeta*omega;    %����������
Kg=omega2;          %����ն���
node=length(fi);
%%%%%%%%%%%%%%%%%����ϵͳABCDϵ����
A=[zeros(nksi)  eye(nksi); -Kg -Cg];
% B3=[zeros(nksi,node);fi'];
% C=blkdiag(eye(nksi),eye(nksi));
% D=zeros(2*nksi,nksi);
% innum=5676;%ʩ�������ڵ���
% inn=3;%ʩ���������ص����ɶ�
save para.mat A Kg Cg fi;