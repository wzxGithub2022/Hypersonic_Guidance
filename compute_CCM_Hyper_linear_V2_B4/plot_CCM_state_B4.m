%%
load('state_CCM_B4.mat');
x = linspace(1,46,46);
%%
% alpha_w = max(sigma_ThBw(:));       %(:)��ʾ��������Ԫ�ز����б�ʾ
% d_bar = alpha_w/lambda;             %��d_barû�г�w�Ͻ磬������J_CCMҲ�����ٳ�w�Ͻ��ˣ�һ����λ
% disp('d_bar'); 
% disp(d_bar);
figure(1)
plot(x,state_CCM(1,:));
title('d bar');
%%
% disp('Control:'); 
% disp(max(d_bar*delta_u(:)));        %���ǿ����Ͻ�
figure(2)
plot(x,state_CCM(2,:));
title('Control');
%%
% disp('W:'); 
% disp(min(min(min(eig_W(:,:,:,1)))));      %С��С
% disp(max(max(max(eig_W(:,:,:,2)))));      %���д�
figure(3)
subplot(121),plot(x,state_CCM(3,:));
title('W min');
subplot(122),plot(x,state_CCM(4,:));
title('W max');
%%
% disp('min eig CCM:'); 
% disp(min(eig_CCM_min(:)));
% disp('max eig CCM:'); 
% disp(max(eig_CCM_max(:)));
figure(4)
subplot(121),plot(x,state_CCM(5,:));
title('R-CCM eig min');
subplot(122),plot(x,state_CCM(6,:));
title('R-CCM eig max');
%%
% disp('euc_bounds');
% disp(d_bar*sqrt(diag(W_upper)));
figure(5)
subplot(231),plot(x,state_CCM(7,:));
title('euc bounds 1');
subplot(232),plot(x,state_CCM(8,:));
title('euc bounds 2');
subplot(233),plot(x,state_CCM(9,:));
title('euc bounds 3');
subplot(234),plot(x,state_CCM(10,:));
title('euc bounds 4');

