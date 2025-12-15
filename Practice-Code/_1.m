% x1n=[1 1 1 1 1 1 1 1 zeros(1,50)];
% x2n=ones(1,128);
% 
% A=[1,-0.9];
% B=[0.05,0.05];
% hn=impz(B,A,58);%impz计算回执数字滤波器的冲激响应
% n=0:length(hn)-1;
% subplot(1,3,1);
% stem(n,hn,'.');
% 
% xlabel('n');
% ylabel('h(n)');
% title('(a)系统单位脉冲响应');
% 
% 
% %qiu x1n=r8(n)的系统响应
% y1n=filter(B,A,x1n);
% 
% n=0:length(y1n)-1;
% 
% subplot(1,3,2);
% stem(n,y1n,'.');
% 
% xlabel('n');
% ylabel('y1(n)');
% title('(b)系统对于R8(n)的响应y1(n)')
% 
% 
% %qiu x2(n)=u(n)的系统响应
% y2n=filter(B,A,x2n);
% 
% n=0:length(y2n)-1;
% subplot(1,3,3);
% stem(n,y2n,'.');
% 
% xlabel('n');
% ylabel('y2(n)');
% title('(c)系统对于u(n)的响应y2(n)')

% 线性卷积法
x1n=[1 1 1 1 1 1 1 1];%R8(n)
h1n=[ones(1,10) zeros(1,10)];
h2n=[1 2.5 2.5 1 zeros(1,10)];

%计算线性卷积
y21n=conv(h1n,x1n);
y22n=conv(h2n,x1n);

n=0:length(y21n)-1;
subplot(1,2,1);
stem(n,y21n,'.');

xlabel('n');
ylabel('y21(n)');
title('(a)系统对于x1(n)的响应y21(n)')

n=0:length(y22n)-1;
subplot(1,2,2);
stem(n,y22n,'.');

xlabel('n');
ylabel('y22(n)');
title('(b)系统对于x1(n)的响应y22(n)')




