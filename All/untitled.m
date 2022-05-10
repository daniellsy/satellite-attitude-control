function f = ballw( K,ki )
%ballw.m��ʾ��ɫС������һ������������˶���ʵʱ����
%����ʾʵʱ�����ĵ��Ը�ʽΪballw(K)
%����ʾʵʱ������������Ƭ�ĵ��Ը�ʽΪf = ballw(K,ki)
%K�����˶���ѭ����������С��1��
%kiָ��������Ƭ��˲�䣬ȡ1��1034֮�����������
%f�洢�������Ƭ���ݣ�����image(f.cdata)�۲���Ƭ
%������յ��˶��켣

t1 = (0:1000)/1000 * 10 * pi;
x1 = cos(t1);
y1 = sin(t1);
z1 = -t1;

t2 = (0:10)/10;
x2 = x1(end) * (1-t2);
y2 = y1(end) * (1-t2);
z2 = z1(end) * ones(size(x2));

t3 = t2;
z3 = (1-t3)* z1(end);
x3 = zeros(size(z3));
y3 = x3;

t4 = t2;
x4 = t4;
y4 = zeros(size(x4));
z4 = y4;
load movedata;
% x = [x1 x2 x3 x4];
% y = [y1 y2 y3 y4];
% z = [z1 z2 z3 z4];
%data = [x',y',z']              %�鿴������ߵ�������ֵ
plot3(x,y,z,'r','Linewidth',0.1)  %��������
hold on;plot3(0,0,0,'b','Linewidth',1)
axis off;                       %����������
%���塰�ߡ�ɫ�����㡱�ͣ��㣩����Ĵ�С��40����������ʽ��xor)
h = line('Color',[0.1 0.1 0.1],'Marker','.','MarkerSize',10,'EraseMode','xor');
%ʹС���˶�
n = length(x);
i = 1;
j = 1;
while 1
    set(h,'xdata',x(i),'ydata',y(i),'zdata',z(i));
    %bw = [x(i),y(i),z(i)]     %�鿴С��λ��
    drawnow;                    %ˢ����Ļ
    pause(0.00001)               %��������
    i = i+1;
    if nargin == 2 && nargout == 1  %�������������Ϊ2�����������1��ʱ��������Ƭ
        if (i == ki && j == 1)
            f = getframe(gcf);      %����i = kiʱ����Ƭ
        end
    end
    if i > n
        i = 1;
        j = j+1;
        if j > K
            break;
        end
    end
end