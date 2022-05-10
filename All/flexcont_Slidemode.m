clear;clc;
load para1.mat;load para2.mat;
load para3.mat;load para4.mat;
load couple1.mat;load couple2.mat;
load couple3.mat;load couple4.mat;
J=[3.497673e6 -2.643113e4 -3.377337e1;...
  -2.643113e4  2.629181e4 6.409160e-1;...
  -3.377337e1  6.409160e-1 3.509507e6];
Jn=J-H*H'-H2*H2'-H3*H3'-H4*H4';
A=zeros(83,83);
mc=mat2cell(A,[3 20 20 20 20],[3 20 20 20 20]);
mc{2,2}=[zeros(10) eye(10);-Kg -Cg];
mc{3,3}=[zeros(10) eye(10);-Kg2 -Cg2];
mc{4,4}=[zeros(10) eye(10);-Kg3 -Cg3];
mc{5,5}=[zeros(10) eye(10);-Kg4 -Cg4];
mc{2,1}=[-H';Cg*H'];
mc{3,1}=[-H2';Cg2*H2'];
mc{4,1}=[-H3';Cg3*H3'];
mc{5,1}=[-H4';Cg4*H4'];
mc{1,1}=[Jn\(-H*Cg*H'-H2*Cg2*H2'-H3*Cg3*H3'-H4*Cg4*H4')];
mc{1,2}=[Jn\(H*Kg) Jn\(H*Cg)];
mc{1,3}=[Jn\(H2*Kg2) Jn\(H2*Cg2)];
mc{1,4}=[Jn\(H3*Kg3) Jn\(H3*Cg3)];
mc{1,5}=[Jn\(H4*Kg4) Jn\(H4*Cg4)];
A=cell2mat(mc);
w0=0.0011;
base_qd=zeros(4,1);
q1=[cos(-10*pi/180/2);0;sin(-10*pi/180/2);0];
q2=[cos(-10*pi/180/2);sin(-10*pi/180/2);0;0];
q3=[q1(1) -q1(2) -q1(3) -q1(4);
    q1(2) q1(1) -q1(4) q1(3);
    q1(3) q1(4) q1(1) -q1(2);
    q1(4) -q1(3) q1(2) q1(1)]*q2;
wbd=[w0;0;0];%����������ٶ�
wb(:,1)=[0.001;0;0];%��ʼ������ٶ�
w(:,1)=[0.001;0;0];
q=[1;0;0;0];%��ʼ��Ԫ��
% innum=5676;%ʩ�������ڵ���
% inn=3;%ʩ���������ص����ɶ�
outnum=15557;%����ڵ�
outn=1;%����ڵ����ɶ�
outnum4=10543;%����ڵ�
outn4=2;%����ڵ����ɶ�
nksi=length(A1)/2;
y(:,1)=zeros(83,1);
y(1:3,1)=w(:,1);
h=0.01;i=1;
 for j=0:h:1000
t(i)=j;
base_qd(1)=cos(w0*t(i)/2);
base_qd(2)=sin(w0*t(i)/2);
qd=[q3(1) -q3(2) -q3(3) -q3(4);
         q3(2) q3(1) -q3(4) q3(3);
         q3(3) q3(4) q3(1) -q3(2);
         q3(4) -q3(3) q3(2) q3(1)]*base_qd;%������Ԫ��

we=wb(:,i)-wbd;
%������ٶȲ�˾���
W=[0 -w(3,i) w(2,i);w(3,i) 0 -w(1,i);-w(2,i) w(1,i) 0];
%��Ԫ�����
qe=[qd(1),qd(2),qd(3),qd(4);
    -qd(2),qd(1),qd(4),-qd(3);
    -qd(3),-qd(4),qd(1),qd(2);
    -qd(4),qd(3),-qd(2),qd(1)]*q;
%ȡʸ������
qe1=qe(2:4);
qe1_dao=-0.5*W*qe1+0.5*qe(1)*we;
%����
miu=1000;
k=eye(3);
%��ģ��
s=we+k*qe1;
M(:,i)=W*J*w(:,i)-Jn*k*qe1_dao-miu*s-H*(Cg*y(14:23,i)+Kg*y(4:13,i)-Cg*H'*w(:,i))-H2*(Cg2*y(34:43,i)+Kg2*y(24:33,i)-Cg2*H2'*w(:,i))...
-H3*(Cg3*y(54:63,i)+Kg3*y(44:53,i)-Cg3*H3'*w(:,i))-H4*(Cg4*y(74:83,i)+Kg4*y(64:73,i)-Cg4*H4'*w(:,i));
% -Jn*G*qv(2:4);
f=Jn\M(:,i);%��Ч����
% f=zeros(3,1);
F=[f;zeros(80,1)];
K1=A*y(:,i)+F;
K2=A*(y(:,i)+h*K1/2)+F;
K3=A*(y(:,i)+h*K2/2)+F;
K4=A*(y(:,i)+h*K3)+F;
y(:,i+1)=y(:,i)+h/6*(K1+2*K2+2*K3+K4);
w(:,i+1)=y(1:3,i+1);

qv=0.5*[-q(2),-q(3),-q(4);
         q(1),q(4),-q(3);
         -q(4),q(1),q(2);
         q(3),-q(2),q(1)]*w(:,i+1);
%����������ٶȣ��ڱ���ϵ�£�
wb(:,i+1)=2*[-q(2),q(1),q(4),-q(3);
             -q(3),-q(4),q(1),q(2);
             -q(4),q(3),-q(2),q(1)]*qv;
q=q+qv*h;
%��������̬�任����
A_OBd=[1-2*qd(3)^2-2*qd(4)^2,2*(qd(2)*qd(3)-qd(1)*qd(4)),2*(qd(2)*qd(4)+qd(1)*qd(3));
    2*(qd(2)*qd(3)+qd(1)*qd(4)),1-2*qd(2)^2-2*qd(4)^2,2*(qd(4)*qd(3)-qd(1)*qd(2));
    2*(qd(2)*qd(4)-qd(1)*qd(3)),2*(qd(4)*qd(3)+qd(1)*qd(2)),1-2*qd(2)^2-2*qd(3)^2];
%ʵ�ʵĻ�����̬�任����(����ڹ�������ϵ)
A_OB=[1-2*q(3)^2-2*q(4)^2,2*(q(2)*q(3)-q(1)*q(4)),2*(q(2)*q(4)+q(1)*q(3));
    2*((2)*q(3)+q(1)*q(4)),1-2*q(2)^2-2*q(4)^2,2*(q(4)*q(3)-q(1)*q(2));
    2*((2)*q(4)-q(1)*q(3)),2*(q(4)*q(3)+q(1)*q(2)),1-2*q(2)^2-2*q(3)^2];
A_dB=A_OBd\A_OB;%�������Ԫ�����ɵ���̬�任����
%���ݱ��������֮�����̬�任�������ŷ����
sita_ze=atan(-A_dB(1,2)/A_dB(1,1));
sita_ye=asin(A_dB(1,3));
sita_xe=-asin(A_dB(2,3)/cos(sita_ye));
%����Ƕȵ�λ
sita_ctrl_eff(:,i)=[sita_xe;sita_ye;sita_ze]*180/pi;
qdsave(:,i)=qd;qsave(:,i)=q;
qplus(:,i)=q(1)^2+q(2)^2+q(3)^2+q(4)^2;
% eff(:,i)=H*r(11:20,:);
i=i+1;
 end
wb(:,i)=[];
y(:,i)=[];
plot(t,wb);
xlabel('t/s');ylabel('rad/s');
z=fi(6*(outnum-1)+3,:)*y(4:13,:);
% x(i)=coord(6*(outnum-1)+outn,:);
% y1(i)=coord(6*(outnum-1)+2,:);
y4=fi4(6*(outnum4-1)+2,:)*y(64:73,:);
figure(2)
plot(t,M);
xlabel('t/s');ylabel('Nm');
figure(3)
plot(t,sita_ctrl_eff);
xlabel('t/s');ylabel('��');
figure(4)
plot(t,z);
xlabel('t/s');ylabel('z/m');
figure(5)
plot(t,y4);
xlabel('t/s');ylabel('y/m');
% plot(t,eff);
% plot(t,qdsave);
% hold on;plot(t,qsave);
