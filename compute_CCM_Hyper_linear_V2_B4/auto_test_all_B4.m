global state_num state_CCM CCW_upper

state_CCM = zeros(10,46);
CCW_upper = zeros(4,4,46);

for state_num = 1:1:46
    fprintf(' state_num = %d , start \n',state_num);
    compute_CCM_Hyper_linear_V2_B4;
    fprintf(' state_num = %d , end \n',state_num);
end

save('state_CCM_B4.mat','state_CCM');
save('CCW_upper_B4.mat','CCW_upper');

CCM_upper = zeros(4,4,46);

for i=1:1:46
    CCM_upper(:,:,i) = inv(CCW_upper(:,:,i));
end

save('CCM_upper_B4.mat','CCM_upper');