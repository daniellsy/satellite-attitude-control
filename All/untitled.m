function f = ballw( K,ki )
%ballw.m演示红色小球沿着一条封闭螺旋线运动的实时动画
%仅演示实时动画的调试格式为ballw(K)
%既演示实时动画又拍摄照片的调试格式为f = ballw(K,ki)
%K红球运动的循环次数（不小于1）
%ki指定拍摄照片的瞬间，取1到1034之间的任意整数
%f存储拍摄的照片数据，可用image(f.cdata)观察照片
%产生封闭的运动轨迹

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
%data = [x',y',z']              %查看封闭曲线的坐标数值
plot3(x,y,z,'r','Linewidth',0.1)  %绘制曲线
hold on;plot3(0,0,0,'b','Linewidth',1)
axis off;                       %不画坐标轴
%定义“线”色、“点”型（点）、点的大小（40）、擦除方式（xor)
h = line('Color',[0.1 0.1 0.1],'Marker','.','MarkerSize',10,'EraseMode','xor');
%使小球运动
n = length(x);
i = 1;
j = 1;
while 1
    set(h,'xdata',x(i),'ydata',y(i),'zdata',z(i));
    %bw = [x(i),y(i),z(i)]     %查看小球位置
    drawnow;                    %刷新屏幕
    pause(0.00001)               %控制球速
    i = i+1;
    if nargin == 2 && nargout == 1  %当输入变量个数为2并且输出变量1个时才拍摄照片
        if (i == ki && j == 1)
            f = getframe(gcf);      %拍摄i = ki时的照片
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