clear;clc;
load needdata3.mat;
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
Mg3=eye(nksi);       %����������
Cg3=2*zeta*omega;    %����������
Kg3=omega2;          %����ն���
node=length(fi);
%%%%%%%%%%%%%%%%%����ϵͳABCDϵ����
A3=[zeros(nksi)  eye(nksi); -Mg3^-1*Kg3 -Mg3^-1*Cg3];
B3=[zeros(nksi,node);Mg3^-1*fi'];
% C=blkdiag(eye(nksi),eye(nksi));
% D=zeros(2*nksi,nksi);
% innum=5676;%ʩ�������ڵ���
% inn=3;%ʩ���������ص����ɶ�
fi3=fi;
save para3.mat A3 Kg3 Cg3 fi3;