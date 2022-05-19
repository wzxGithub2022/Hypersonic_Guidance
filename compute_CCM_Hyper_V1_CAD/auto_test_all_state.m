global state_num state_CCM

state_CCM = zeros(11,46);

for state_num = 1:1:46
    fprintf(' state_num = %d , start \n',state_num);
    compute_CCM_Hyper_V1_CAD;
    fprintf(' state_num = %d , end \n',state_num);
end

% save('state_CCM.mat','state_CCM');