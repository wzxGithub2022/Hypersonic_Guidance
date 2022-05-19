%%
clear; clc;
%%
sin_x = @(x) 0.9101*(x./(pi/3)) - 0.04466*(4*(x./(pi/3)).^3 - 3*(x./(pi/3)));
cos_x = @(x) 0.7441 -0.2499*(2*(x./(pi/3)).^2 -1);

x = linspace(-2,2,100);
y_sin = sin_x(x);
y_cos = cos_x(x);
figure(1),
subplot(121),plot(x,y_sin,':x',x,sin(x),'-.s');
title('切比雪夫三阶多项式逼近效果');
xlabel('x'),ylabel('y');
legend('chebysheff','sin');
subplot(122),plot(x,y_cos,'-o',x,cos(x),'--+','LineWidth',2);
title('切比雪夫三阶多项式逼近效果');
xlabel('x'),ylabel('y');
legend('chebysheff','cos');
