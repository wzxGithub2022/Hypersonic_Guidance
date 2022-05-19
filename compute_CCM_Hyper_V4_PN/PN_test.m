alpha_lim = deg2rad(5);
alpha = linspace(-alpha_lim,alpha_lim,100);
figure(1)
for M=2:0.1:4
    CL  = 0.4172 + 19.41*alpha + 10.17*alpha.^2 - M*(0.1004 + 0.7536*alpha);
    hold on
    plot(alpha,CL);
    legend
end