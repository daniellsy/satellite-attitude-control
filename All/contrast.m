lear;clc;
load datacollab.mat
wbc=wb; Mc=M; sita_ctrl_effc=sita_ctrl_eff; zc=z;
load dataslide.mat;
plot(t,wb(1,:),'k',t,wb(2,:),'--',t,wb(3,:),'-.','LineWidth',1);
hold on
plot(t,wbc(1,:),'k',t,wbc(2,:),'--',t,wbc(3,:),'-.','LineWidth',1.2);
xlabel('t/s');ylabel('rad/s');title('������ٶ�');
legend('wx','wy','wz');
% x(i)=coord(6*(outnum-1)+outn,:);
% y1(i)=coord(6*(outnum-1)+2,:);
figure(2)
plot(t,M(1,:),'k',t,M(2,:),'--',t,M(3,:),'-.','LineWidth',0.2);
% hold on
% plot(t,Mc(1,:),'k',t,Mc(2,:),'--',t,Mc(3,:),'-.','LineWidth',1.2);
xlabel('t/s');ylabel('Nm');title('��̬��������');
legend('Mx','My','Mz');
figure(3)
plot(t,sita_ctrl_eff(1,:),'k',t,sita_ctrl_eff(2,:),'--',t,sita_ctrl_eff(3,:),'-.','LineWidth',1.2);
xlabel('t/s');ylabel('��');title('ŷ���Ǳ仯����');
legend('ƫ����','��ת��','������');
figure(4)
plot(t,z,'--k','LineWidth',0.6);
hold on
plot(t,zc,'r','LineWidth',1);
xlabel('t/s');ylabel('z/m');title('���ĩ��z��λ��');