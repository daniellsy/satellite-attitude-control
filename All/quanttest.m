clear;clc;
a1=sin(4/180*pi);b1=cos(4/180*pi);%theta
a2=sin(5/180*pi);b2=cos(5/180*pi);%fi
a3=sin(3/180*pi);b3=cos(3/180*pi);%psi
q0=b1*b2*b3-a1*a2*a3;
q1=b1*a2*b3-a1*b2*a3;
q2=a1*b2*b3+b1*a2*a3;
q3=a1*a2*b3+b1*b2*a3;
q=[q0;q1;q2;q3];%初始四元数
fii=asin(2*(q(3)*q(4)+q(1)*q(2)));
theta=atan(2*(q(1)*q(3)-q(4)*q(2))/(q(4)^2+q(1)^2-q(2)^2-q(3)^2));
psi=atan(2*(q(1)*q(4)-q(3)*q(2))/(q(3)^2+q(1)^2-q(2)^2-q(4)^2));
sita_ctrl_eff=[fii;theta;psi]*180/pi;
% fi=asin(2*(q2*q3+q0*q1))*180/pi;
% theta=atan(2*(q0*q2-q3*q1)/(q3^2+q0^2-q1^2-q2^2))*180/pi;
% psi=atan(2*(q0*q3-q2*q1)/(q2^2+q0^2-q1^2-q3^2))*180/pi;