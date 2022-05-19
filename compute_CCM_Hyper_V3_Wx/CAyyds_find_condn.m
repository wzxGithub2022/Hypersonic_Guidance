function [condn] = CAyyds_find_condn(n,X_d,condn_prev,lambda,ccm_eps)

eps = 1;
return_metric = 0;
solved = 0;

%Determine upper bound�������Ͻ�
cond_l = condn_prev;
cond_u = 1.2*condn_prev;

while (~solved)
    fprintf(' cond_u: %.2f: ', cond_u);
    [sos_prob,~,~,~] = CAyyds_CCM_Opt(n,X_d,cond_u,lambda,ccm_eps,return_metric);
    if (sos_prob == 0)
        solved = 1;
        fprintf('feasible \n');
    else
        %shift up condition number range,������˵Ļ���һ������condition number��ֻ����һ���������
        fprintf('�㿴��ⲻ�����˰�\n');
        cond_l = cond_u;
        cond_u = 1.2*cond_u;
    end
end

if (solved)           
    fprintf(' cond_l: %.2f, cond_u: %.2f\n', cond_l, cond_u);
end

%Now do bisection search������������
while(cond_u - cond_l >= eps)     %���ϴ�ϸ���Ż�
    condn = (cond_l+cond_u)/2;
    fprintf(' cond: %.2f', condn);
    [sos_prob,~,~,~] = CAyyds_CCM_Opt(n,X_d,condn,lambda,ccm_eps,return_metric);
    if (sos_prob == 0)
        fprintf(' feasible\n');       %��Ϳ����ˣ�������ȡСһ�㣬���󶯵Ļ���euc��d�ͻ������Ż�һЩ
        cond_u = condn;
    else
        fprintf(' infeasible\n');     %�����еĻ������ҷŴ�һ��
        cond_l = condn;
    end
end

condn = cond_u;      %����С��ȫ�����

end