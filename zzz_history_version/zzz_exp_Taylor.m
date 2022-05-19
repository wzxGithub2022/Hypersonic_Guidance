%%
clear; clc;
%%
x = linspace(-3,0,100);     %-3到0就反映了21km到地面，这逼近得可以了
y_exp = exp(-x);
hf = -x;
% y_poly = 1 + hf + hf.^2/2 + hf.^3/6 + hf.^4/24 + hf.^5/120 + hf.^6/720 + hf.^7/5040 + hf.^8/40320 + hf.^9/362880 + hf.^10/3628800;
y_poly = 1 + hf + hf.^2/2 + hf.^3/6 + hf.^4/24 + hf.^5/120 + hf.^6/720;
figure(1)
plot(x,y_exp,x,y_poly);
legend('exp','Taylor');