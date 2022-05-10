clear;clc;
load para1.mat;
load para2.mat;
load para3.mat;
load para4.mat;
load couple1.mat;
load couple2.mat;
load couple3.mat;
load couple4.mat;
Jn=[3.497673e6 -2.643113e4 -3.377337e1;...
  -2.643113e4  2.629181e4 6.409160e-1;...
  -3.377337e1  6.409160e-1 3.509507e6];
w0=0.0011;kesi=1;tn=300;
base_qd=zeros(4,1);
q1=[cos(-10*pi/180/2);0;sin(-10*pi/180/2);0];
q2=[cos(-10*pi/180/2);sin(-10*pi/180/2);0;0];
q3=[q1(1) -q1(2) -q1(3) -q1(4);
    q1(2) q1(1) -q1(4) q1(3);
    q1(3) q1(4) q1(1) -q1(2);
    q1(4) -q1(3) q1(2) q1(1)]*q2;
wbd=[w0;0;0];%期望本体角速度
w(:,1)=[0.001;0;0];%初始本体角速度
q=[1;0;0;0];%初始四元数
% innum=5676;%施加外力节点编号
% inn=3;%施加力或力矩的自由度
outnum=15557;%输出节点
outn=1;%输出节点自由度
outnum4=10543;%输出节点
outn4=2;%输出节点自由度
nksi=length(A1)/2;
b1r=zeros(2*nksi,1);b1r2=zeros(10,1);
b2r=zeros(2*nksi,1);b2r2=zeros(10,1);
b3r=zeros(2*nksi,1);b3r2=zeros(10,1);
b4r=zeros(2*nksi,1);b4r2=zeros(10,1);
h=0.01;i=1;
 for j=0:h:400
t(i)=j;
base_qd(1)=cos(w0*t(i)/2);
base_qd(2)=sin(w0*t(i)/2);
qd=[q3(1) -q3(2) -q3(3) -q3(4);
         q3(2) q3(1) -q3(4) q3(3);
         q3(3) q3(4) q3(1) -q3(2);
         q3(4) -q3(3) q3(2) q3(1)]*base_qd;%期望四元数
we=w(:,i)-wbd;
%构造角速度叉乘矩阵
W=[0 -w(3,i) w(2,i);w(3,i) 0 -w(1,i);-w(2,i) w(1,i) 0];
%四元数误差
qe=[qd(1),qd(2),qd(3),qd(4);
    -qd(2),qd(1),qd(4),-qd(3);
    -qd(3),-qd(4),qd(1),qd(2);
    -qd(4),qd(3),-qd(2),qd(1)]*q;
%取矢量部分
qe1=qe(2:4);
%%计算固有频率
wn=10/(kesi*tn);
%%%%%%%%%%%%%%控制率设计
%参数选取与计算
k=wn^2*2;
d=2*kesi*wn;
D=d*Jn;
K=k*Jn;
M(:,i)=(H*b1r2+H2*b2r2+H3*b3r2+H4*b4r2)/4+W*Jn*w(:,i)-K*qe1-D*we;
% M(:,i)=W*Jn*w(:,i)-K*qe1-D*we;
% K=50*eye(3,3);D=100*eye(3,3);
% M(:,i)=-K*qe1+D*we;
k1=Jn\(-cross(w(:,i),Jn*w(:,i))-H*b1r2-H2*b2r2-H3*b3r2-H4*b4r2+M(:,i));
k2=Jn\(-cross((w(:,i)+h*k1/2),Jn*(w(:,i)+h*k1/2))-H*b1r2-H2*b2r2-H3*b3r2-H4*b4r2+M(:,i));
k3=Jn\(-cross((w(:,i)+h*k2/2),Jn*(w(:,i)+h*k2/2))-H*b1r2-H2*b2r2-H3*b3r2-H4*b4r2+M(:,i));
k4=Jn\(-cross((w(:,i)+h*k3),Jn*(w(:,i)+h*k3))-H*b1r2-H2*b2r2-H3*b3r2-H4*b4r2+M(:,i));
% k1=Jn\(-cross(w(:,i),(Jn*w(:,i)+H*r(11:20,:)))-H*r2+M(:,i));
% k2=Jn\(-cross((w(:,i)+h*k1/2),(Jn*(w(:,i)+h*k1/2)+H*r(11:20,:)))-H*r2+M(:,i));
% k3=Jn\(-cross((w(:,i)+h*k2/2),(Jn*(w(:,i)+h*k2/2)+H*r(11:20,:)))-H*r2+M(:,i));
% k4=Jn\(-cross((w(:,i)+h*k3),(Jn*(w(:,i)+h*k3)+H*r(11:20,:)))-H*r2+M(:,i));
w(:,i+1)=w(:,i)+h/6*(k1+2*k2+2*k3+k4);
% wv(:,i)=h/6*(k1+2*k2+2*k3+k4);
wv(:,i)=Jn\(-cross(w(:,i+1),Jn*w(:,i+1))-H*b1r2-H2*b2r2-H3*b3r2-H4*b4r2+M(:,i));
f=-H'*wv(:,i);%附件等效外力
f2=-H2'*wv(:,i);
f3=-H3'*wv(:,i);
f4=-H4'*wv(:,i);
F=[zeros(nksi,1);f];
F2=[zeros(nksi,1);f2];
F3=[zeros(nksi,1);f3];
F4=[zeros(nksi,1);f4];

K1=A1*b1r+F;%第一附件
K2=A1*(b1r+h*K1/2)+F;
K3=A1*(b1r+h*K2/2)+F;
K4=A1*(b1r+h*K3)+F;
b1r=b1r+h/6*(K1+2*K2+2*K3+K4);
b1r1=b1r(1:nksi,:);b1r3=b1r(11:20,:);
b1r2=-Kg*b1r1-H'*wv(:,i)-Cg*b1r3;

coord=fi*b1r1;
x(i)=coord(6*(outnum-1)+outn,:);
y(i)=coord(6*(outnum-1)+2,:);
z(i)=coord(6*(outnum-1)+3,:);

b2K1=A2*b2r+F2;%第二附件
b2K2=A2*(b2r+h*b2K1/2)+F2;
b2K3=A2*(b2r+h*b2K2/2)+F2;
b2K4=A2*(b2r+h*b2K3)+F2;
b2r=b2r+h/6*(b2K1+2*b2K2+2*b2K3+b2K4);
b2r1=b2r(1:nksi,:);b2r3=b2r(11:20,:);
b2r2=-Kg2*b2r1-H2'*wv(:,i)-Cg2*b2r3;

b3K1=A3*b3r+F3;%第三附件
b3K2=A3*(b3r+h*b3K1/2)+F3;
b3K3=A3*(b3r+h*b3K2/2)+F3;
b3K4=A3*(b3r+h*b3K3)+F3;
b3r=b3r+h/6*(b3K1+2*b3K2+2*b3K3+b3K4);
b3r1=b3r(1:nksi,:);b3r3=b3r(11:20,:);
b3r2=-Kg3*b3r1-H3'*wv(:,i)-Cg3*b3r3;

b4K1=A4*b4r+F4;%第四附件
b4K2=A4*(b4r+h*b4K1/2)+F4;
b4K3=A4*(b4r+h*b4K2/2)+F4;
b4K4=A4*(b4r+h*b4K3)+F4;
b4r=b4r+h/6*(b4K1+2*b4K2+2*b4K3+b4K4);
b4r1=b4r(1:nksi,:);b4r3=b4r(11:20,:);
b4r2=-Kg4*b4r1-H4'*wv(:,i)-Cg4*b4r3;

coord=fi4*b4r1;
z4(i)=coord(6*(outnum4-1)+outn4,:);

qv=0.5*[q(1),-q(2),-q(3),-q(4);
         q(2),q(1),-q(4),q(3);
         q(3),q(4),q(1),-q(2);
         q(4),-q(3),q(2),q(1)]*[0;w(:,i+1)];
q=q+qv*h;
%期望的姿态变换矩阵
A_OBd=[1-2*qd(3)^2-2*qd(4)^2,2*(qd(2)*qd(3)-qd(1)*qd(4)),2*(qd(2)*qd(4)+qd(1)*qd(3));
    2*(qd(2)*qd(3)+qd(1)*qd(4)),1-2*qd(2)^2-2*qd(4)^2,2*(qd(4)*qd(3)-qd(1)*qd(2));
    2*(qd(2)*qd(4)-qd(1)*qd(3)),2*(qd(4)*qd(3)+qd(1)*qd(2)),1-2*qd(2)^2-2*qd(3)^2];
%实际的基体姿态变换矩阵(相对于惯性坐标系)
A_OB=[1-2*q(3)^2-2*q(4)^2,2*(q(2)*q(3)-q(1)*q(4)),2*(q(2)*q(4)+q(1)*q(3));
    2*((2)*q(3)+q(1)*q(4)),1-2*q(2)^2-2*q(4)^2,2*(q(4)*q(3)-q(1)*q(2));
    2*((2)*q(4)-q(1)*q(3)),2*(q(4)*q(3)+q(1)*q(2)),1-2*q(2)^2-2*q(3)^2];
A_dB=A_OBd\A_OB;%由误差四元数构成的姿态变换矩阵
%根据本体和期望之间的姿态变换矩阵求解欧拉角
sita_ze=atan(-A_dB(1,2)/A_dB(1,1));
sita_ye=asin(A_dB(1,3));
sita_xe=-asin(A_dB(2,3)/cos(sita_ye));
%输出角度单位
sita_ctrl_eff(:,i)=[sita_xe;sita_ye;sita_ze]*180/pi;
qdsave(:,i)=qd;qsave(:,i)=q;
% eff(:,i)=H*r(11:20,:);
i=i+1;
 end
w(:,i)=[];
plot(t,w);
figure(2)
plot(t,M);
figure(3)
plot(t,z);
figure(4)
plot(t,sita_ctrl_eff);
figure(5)
plot(t,z4);
% plot(t,eff);
% plot(t,qdsave);
% hold on;plot(t,qsave);
qplus=q(1)^2+q(2)^2+q(3)^2+q(4)^2;