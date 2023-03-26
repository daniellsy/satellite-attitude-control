clear;clc;
load para.mat;
Af=A;
J=[5600 0 0;...
  0 5400 0;...
  0 0 4300];
H=[-19.911 0 5.311 -0.094 2.217 0 0 0 0 0;
    0.745 0.039 -0.375 1.441 -0.197 0 0 0 0 0;
    0 -20.738 0 0 0 0 0 0 0 0];
Jn=J-H*H';
Jnn=inv(Jn);
A=zeros(23,23);
mc=mat2cell(A,[3 20],[3 20]);
mc{2,2}=[zeros(10) eye(10);-Kg -Cg];
mc{2,1}=[-H';Cg*H'];
mc{1,1}=Jn\(-H*Cg*H');
mc{1,2}=[Jn\(H*Kg) Jn\(H*Cg)];
A=cell2mat(mc);
theta=8;fii=10;psi=-6;
a1=sin(theta/2/180*pi);b1=cos(theta/2/180*pi);%theta
a2=sin(fii/2/180*pi);b2=cos(fii/2/180*pi);%fi
a3=sin(psi/2/180*pi);b3=cos(psi/2/180*pi);%psi
q0=b1*b2*b3-a1*a2*a3;
q1=b1*a2*b3-a1*b2*a3;q2=a1*b2*b3+b1*a2*a3;
q3=a1*a2*b3+b1*b2*a3;
wbd=[0;0;0];%期望本体角速度
wb(:,1)=[0;0;0];%初始本体角速度
w(:,1)=[0;0;0];%本体相对于惯性系速度
q=[q0;q1;q2;q3];%初始四元数
qd=[1;0;0;0];%期望四元数
%参数
miu=450;
k=eye(3);
% innum=5676;%施加外力节点编号
% inn=3;%施加力或力矩的自由度
outnum=163;%输出节点
outn=1;%输出节点自由度
nksi=length(Af)/2;
y(:,1)=zeros(23,1);
y(1:3,1)=w(:,1);
h=0.01;i=1;
 for j=0:h:100
t(i)=j;
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
fii=asin(2*(q(3)*q(4)+q(1)*q(2)));
theta=atan(2*(q(1)*q(3)-q(4)*q(2))/(q(4)^2+q(1)^2-q(2)^2-q(3)^2));
psi=atan(2*(q(1)*q(4)-q(3)*q(2))/(q(3)^2+q(1)^2-q(2)^2-q(4)^2));
sita_ctrl_eff(:,i)=[fii;theta;psi]*180/pi;
%取矢量部分
qe1=qe(2:4);
qe1_dao=-0.5*W*qe1+0.5*qe(1)*we;

%滑模面
s=we+k*qe1;
M(:,i)=W*J*w(:,i)-Jn*k*qe1_dao-miu*s-H*(Cg*y(14:23,i)+Kg*y(4:13,i)-Cg*H'*w(:,i));
% -Jn*G*qv(2:4);
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
i=i+1;
 end
wb(:,i)=[];
y(:,i)=[];
plot(t,wb);
xlabel('t/s');ylabel('rad/s');
z=fi(6*(outnum-1)+outn,:)*y(4:13,:);
y1=fft(z);
f1=(0:length(y1)-1)*100/length(y1);
plot(t,wb(1,:),'k',t,wb(2,:),'--',t,wb(3,:),'-.','LineWidth',1.2);
xlabel('t/s');ylabel('rad/s');title('本体角速度');
legend('wx','wy','wz');
% x(i)=coord(6*(outnum-1)+outn,:);
% y1(i)=coord(6*(outnum-1)+2,:);
figure(2)
plot(t,M(1,:),'k',t,M(2,:),'--',t,M(3,:),'-.','LineWidth',1.2);
xlabel('t/s');ylabel('Nm');title('姿态控制力矩');
legend('Mx','My','Mz');
figure(3)
plot(t,sita_ctrl_eff(1,:),'k',t,sita_ctrl_eff(2,:),'--',t,sita_ctrl_eff(3,:),'-.','LineWidth',1.2);
xlabel('t/s');ylabel('°');title('欧拉角变化曲线');
legend('偏航角','滚转角','俯仰角');
figure(4)
plot(t,z);
xlabel('t/s');ylabel('z/m');title('桁架末端z向位移');
figure(5)
plot(f1,abs(y1));
save dataslide.mat wb M sita_ctrl_eff z t