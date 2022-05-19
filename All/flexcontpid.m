clear;clc;
load para.mat;load couple.mat;load delta.mat;
Af=A;
% J=[210 0 0;...
%   0  4200 0;...
%   0  0 4800];
% J=[3.497673e6 -2.643113e4 -3.377337e1;...
%   -2.643113e4  2.629181e4 6.409160e-1;...
%   -3.377337e1  6.409160e-1 3.509507e6];
J=inv(eye(3)-delta)*(H*H');
Jn=J-H*H';
Jnn=inv(Jn);
A=zeros(23,23);
mc=mat2cell(A,[3 20],[3 20]);
mc{2,2}=[zeros(10) eye(10);-Kg -Cg];
mc{2,1}=[-H';Cg*H'];
mc{1,1}=Jn\(-H*Cg*H');
mc{1,2}=[Jn\(H*Kg) Jn\(H*Cg)];
A=cell2mat(mc);
w0=0.0011;kesi=1;tn=100;
%计算固有频率
wn=10/(kesi*tn);
%%%%%%%%%%%%%%控制率设计
%参数选取与计算
k=wn^2*2;
d=2*kesi*wn;
D=d*J;
K=k*J;
% D=500*eye(3);K=30*eye(3);
base_qd=zeros(4,1);
q1=[cos(-10*pi/180/2);0;sin(-10*pi/180/2);0];
q2=[cos(-10*pi/180/2);sin(-10*pi/180/2);0;0];
q3=[q1(1) -q1(2) -q1(3) -q1(4);
    q1(2) q1(1) -q1(4) q1(3);
    q1(3) q1(4) q1(1) -q1(2);
    q1(4) -q1(3) q1(2) q1(1)]*q2;
wbd=[w0;0;0];%期望本体角速度
wb(:,1)=[0.001;0;0];%初始本体角速度
w(:,1)=[0.001;0;0];%本体相对于惯性系速度
q=[1;0;0;0];%初始四元数
% innum=5676;%施加外力节点编号
% inn=3;%施加力或力矩的自由度
outnum=100;%输出节点
outn=1;%输出节点自由度
outnum4=163;%输出节点
outn4=2;%输出节点自由度
nksi=length(Af)/2;
y(:,1)=zeros(23,1);
y(1:3,1)=w(:,1);
h=0.01;i=1;
 for j=0:h:100
t(i)=j;
base_qd(1)=cos(w0*t(i)/2);
base_qd(2)=sin(w0*t(i)/2);
qd=[q3(1) -q3(2) -q3(3) -q3(4);
         q3(2) q3(1) -q3(4) q3(3);
         q3(3) q3(4) q3(1) -q3(2);
         q3(4) -q3(3) q3(2) q3(1)]*base_qd;%期望四元数

we=wb(:,i)-wbd;
%构造角速度叉乘矩阵
W=[0 -w(3,i) w(2,i);w(3,i) 0 -w(1,i);-w(2,i) w(1,i) 0];
%四元数误差
qe=[qd(1),qd(2),qd(3),qd(4);
    -qd(2),qd(1),qd(4),-qd(3);
    -qd(3),-qd(4),qd(1),qd(2);
    -qd(4),qd(3),-qd(2),qd(1)]*q;
%取矢量部分
qe1=qe(2:4);
% D=20000*eye(3);
% K=10000*eye(3);
% M(:,i)=(H*b1r2+H2*b2r2+H3*b3r2+H4*b4r2)/4+W*Jn*w(:,i)-K*qe1-D*we;
M(:,i)=W*J*w(:,i)-K*qe1-D*we;
f=Jn\M(:,i);%等效外力
% f=zeros(3,1);
F=[f;zeros(20,1)];
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
%计算星体角速度（在本体系下）
wb(:,i+1)=2*[-q(2),q(1),q(4),-q(3);
             -q(3),-q(4),q(1),q(2);
             -q(4),q(3),-q(2),q(1)]*qv;
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
qplus(:,i)=q(1)^2+q(2)^2+q(3)^2+q(4)^2;
% eff(:,i)=H*r(11:20,:);
i=i+1;
 end
wb(:,i)=[];
y(:,i)=[];
plot(t,wb);
xlabel('t/s');ylabel('rad/s');%title('随时间变化曲线');
z=fi(6*(outnum-1)+3,:)*y(4:13,:);
% x(i)=coord(6*(outnum-1)+outn,:);
% y1(i)=coord(6*(outnum-1)+2,:);
figure(2)
plot(t,M);
xlabel('t/s');ylabel('Nm');
figure(3)
plot(t,sita_ctrl_eff);
xlabel('t/s');ylabel('°');
figure(4)
plot(t,z);
xlabel('t/s');ylabel('z/m');
% plot(t,eff);
% plot(t,qdsave);
% hold on;plot(t,qsave);
