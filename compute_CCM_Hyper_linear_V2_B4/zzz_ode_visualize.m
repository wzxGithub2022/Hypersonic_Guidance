clear;clc;
for t0 = 0:0.1:3

%���н��ȿ��ӻ�
if (t0 < 1)
    t0_disp = 1;
else
    if (t0 > 1 && t0 > t0_disp)
        fprintf('t0 calculating %f second \n',t0_disp);
        t0_disp = t0_disp + 1;
    end
end

end