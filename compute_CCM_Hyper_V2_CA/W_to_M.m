load('CCW_upper_CA.mat')

CCM_upper = zeros(4,4,46);

for i=1:1:46
    CCM_upper(:,:,i) = inv(CCW_upper(:,:,i));
end

save('CCM_upper_CA.mat','CCM_upper');