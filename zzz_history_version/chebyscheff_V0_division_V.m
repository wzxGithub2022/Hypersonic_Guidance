%% 使用三阶切比雪夫多项式在区间[1000,2000]上逼近1/V
% f(x) = c0/2 + sum(ck*Tk(x)); k = 1,2,...
% ck = 2/pi * integral(f(cos(theta)) * cos(k*theta)) d_theta; [0,pi]
% x = 500t + 1500; t belongs to [-1,+1]
% t = x/500 - 3;

clear;clc;
syms theta
c0_half = (1/pi) * int( 1 / (500*cos(theta)+1500) ,0,pi);
c0h_double = double(c0_half);

for k = 1:3
%     fprintf(' k = %d ',k);
    c(k) = (2/pi) * int(cos(k*theta) / (500*cos(theta)+1500) ,0,pi);
    c_double(k) = double(c(k));
end

% 0.001414213562373 7.071067811865476e-04 
% -2.426406871192852e-04 4.163056034261583e-05 -7.142674936409832e-06

% T0_x = 1;
% T1_x = @(x) x;
% T2_x = @(x) 2*x^2 - 1;
% T3_x = @(x) 4*x^3 - 3*x;
% T4_x = @(x) 8*x^4 - 8*x^2 + 1;

% f_V = @(x) 7.0711e-04 - 2.4264e-04 * x + 4.1631e-05 * (2*x.^2 - 1) + -7.1427e-06 * (4*x.^3 - 3*x);

% f_V = @(x) 7.0711e-04 +x-x - 2.4264e-04 * x;

% x = linspace(1000,2000,100);
% V_division = f_V(x);
% y = 1./x;
% figure(1)
% plot(x,V_division,x,y);
% legend('V_division','1/x')

f_V = @(t) 7.0711e-04 - 2.4264e-04 * t + 4.1631e-05 * (2*t.^2 - 1) + -7.1427e-06 * (4*t.^3 - 3*t);

t = linspace(-1,1,100);
x = 500*t + 1500;
V_division = f_V(t);
y = 1./x;
figure(1)
subplot(121),plot(x,V_division,x,y);
legend('1/V chebysheff','1/x original');
title('切比雪夫三阶多项式逼近效果');

%%
% close all;clear;clc;

syms x t f_V

t = x/500 - 3;
f_V = 7.0711e-04 - 2.4264e-04 * t + 4.1631e-05 * (2*t.^2 - 1) + -7.1427e-06 * (4*t.^3 - 3*t);

f_V_simple = expand(f_V);

% expand(f_V, 'ArithmeticOnly', true)

% disp(f_V);
 
% f_V_simple =
%  
% - (4216305884649127*x^3)/18446744073709551616000000000
% + (50234055402363781*x^2)/36893488147419103232000000
% - (880836751886113983*x)/295147905179352825856000
% + 1682274762297052657/590295810358705651712

a3 = - 4216305884649127/18446744073709551616000000000;  % -2.2857e-13 -2.285664000000000e-13
a2 = + 50234055402363781/36893488147419103232000000;    % 1.3616e-09 1.361596800000000e-09
a1 = - 880836751886113983/295147905179352825856000;     % -2.9844e-06 -2.984391000000000e-06
a0 = + 1682274762297052657/590295810358705651712;       % 0.0028 0.002849884300000

x = linspace(1000,2000,100);
g_V = a0 + a1*x + a2*x.^2 + a3*x.^3;

g_V_simple = 0.0028498843 - 2.984391e-06*x + 1.3615968e-09*x.^2 - 2.285664e-13*x.^3;

y = 1./x;
figure(1)
subplot(122),plot(x,g_V_simple,x,g_V,x,y);
legend('g(V) simple','g(V) 1/V chebysheff','1/x original');
title('切比雪夫三阶多项式逼近效果');

figure(2)
subplot(121),plot(x,g_V_simple - g_V);
title('二次计算误差');
subplot(122),plot(x,g_V_simple - y);
title('一次逼近误差');