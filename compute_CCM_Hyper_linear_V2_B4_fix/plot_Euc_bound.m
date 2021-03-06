lambda = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
PVTOL = [71.15,35.57,25.68,20.87,18,16.15,14.9497,14.07,13.43,13.06];
lambda_B4 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
linear_B4 = [32.19,18.91,15.13,13.16,12.19,11.67,11.32,11.10,10.91];

figure(1)
plot(lambda,PVTOL,'-o','Color',[0 0.447 0.741],'Linewidth',2);
hold on,plot(lambda_B4,linear_B4,'--*','Color',[0.85 0.325 0.098],'Linewidth',2);
legend('参考规律','实验数据');
xlabel('收缩率\lambda');
ylabel('性能指标J');
title('性能指标J随收缩率\lambda的变化图');