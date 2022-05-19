function [condn] = CAyyds_find_condn(n,X_d,condn_prev,lambda,ccm_eps)

eps = 1;
return_metric = 0;
solved = 0;

%Determine upper bound，决定上界
cond_l = condn_prev;
cond_u = 1.2*condn_prev;

while (~solved)
    fprintf(' cond_u: %.2f: ', cond_u);
    [sos_prob,~,~,~] = CAyyds_CCM_Opt(n,X_d,cond_u,lambda,ccm_eps,return_metric);
    if (sos_prob == 0)
        solved = 1;
        fprintf('feasible \n');
    else
        %shift up condition number range,解决不了的话就一点点调大condition number，只调大不一定解决问题
        fprintf('你看这解不出来了吧\n');
        cond_l = cond_u;
        cond_u = 1.2*cond_u;
    end
end

if (solved)           
    fprintf(' cond_l: %.2f, cond_u: %.2f\n', cond_l, cond_u);
end

%Now do bisection search，做二分搜索
while(cond_u - cond_l >= eps)     %误差较大，细分优化
    condn = (cond_l+cond_u)/2;
    fprintf(' cond: %.2f', condn);
    [sos_prob,~,~,~] = CAyyds_CCM_Opt(n,X_d,condn,lambda,ccm_eps,return_metric);
    if (sos_prob == 0)
        fprintf(' feasible\n');       %这就可行了，向左再取小一点，向左动的话，euc和d就还能再优化一些
        cond_u = condn;
    else
        fprintf(' infeasible\n');     %不可行的话，向右放大一点
        cond_l = condn;
    end
end

condn = cond_u;      %误差较小，全都相等

end