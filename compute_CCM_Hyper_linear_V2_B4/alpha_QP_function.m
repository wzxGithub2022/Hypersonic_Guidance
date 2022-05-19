function [k] = alpha_QP_function(t,state1,alpha_std)
%% B4��CAֻ��ֻ��QP�ǲ�һ���ģ���ӿ������Ĺ켣�����Ƕ�һ��
% ˼�룺
% Ҫ��ƫ���׼_std�͵�ǰ
% Ҫ�ж���ѧ���÷��ź�������ʽд���������Է�������
% ��CCW��CCM�����棬�ٲ�ֵ
% ���룺
% ��ʵ����t ��һ��ʱ�� state1 ��һ��״̬ 
% ϵͳ��������ʱ���״̬(��ֵ)�ͱ�ƵĿ���(alpha_std)
% lambda�ǹ̶�0.5�ܵģ�CCM��ֵ�õ���ǰʱ���CCM
%% ��õ�ǰ״̬�ͱ��״̬
global R0 g t_nor y_nor h_nor v_nor gama_nor alpha_nor CCM_nor B_nor df_nor
%��ʵ״̬
y1 = state1(1);
h1 = state1(2);
v1 = state1(3);
gama1 = state1(4);
alpha1 = state1(5);
%���״̬
t_org = t_nor * (sqrt(R0/g));  %��չ���ʱ�����У�������һ��
t0 = t * (sqrt(R0/g));         %��չ��ǰ��һ��ʱ�䣬�õ���ǰ��ʵʱ��
y_std= interp1(t_org,y_nor,t0,'spline');
h_std= interp1(t_org,h_nor,t0,'spline');
v_std= interp1(t_org,v_nor,t0,'spline');
gama_std= interp1(t_org,gama_nor,t0,'spline');
% alpha_stdֱ����
CCM_std = zeros(4,4);
for ccm_line = 1:4
    for ccm_row = 1:4
        ccm_ij_temp = reshape( CCM_nor(ccm_line,ccm_row,:),46,1 );  %reshape�Ǹ��ö���
        CCM_std(ccm_line,ccm_row) = interp1(t_org,ccm_ij_temp,t0,'spline');
    end
end
%% �Ȳ�����ʡ��һ�������ٶ���
B_std = zeros(4,1);
for B_num = 1:4
        B_ij_temp = reshape( B_nor(B_num,1,:),46,1 );  %reshape�Ǹ��ö���
        B_std(B_num,1) = interp1(t_org,B_ij_temp,t0,'spline');
end
df_std = zeros(4,4);
for df_line = 1:4
    for df_row = 1:4
        df_ij_temp = reshape( df_nor(df_line,df_row,:),46,1 );  %reshape�Ǹ��ö���
        df_std(ccm_line,ccm_row) = interp1(t_org,df_ij_temp,t0,'spline');
    end
end
%% ���ڿ��Կ���д��QP����ʽ�ˣ���ǰ�Ĺ�����ô��ã���д״̬

x1 = [y1;h1;v1;gama1];
x_std = [y_std;h_std;v_std;gama_std];
x_e = x1 - x_std;

if norm(x_e) == 0
    k = 0;              %û�����ֱ���������찡
else
    delta_gama = x_e / norm(x_e);   %���û����norm��0����Ϊ�����������޷�����
    
    % CCM_std�Ѿ�д����
    lambda = 0.5;
    
    %��С��Ŀ�꺯��
    H = 2;
    f = 0;
    
    %����ʽԼ��
    A = delta_gama' * CCM_std * B_std;
    b = -lambda * x_e' * CCM_std * x_e - delta_gama' * CCM_std * df_std * x_e;
    
    %û�е�ʽԼ��
    Aeq = [];
    beq = [];
    %û�г�ֵ�²�
    lb = [];
    ub = [];
    
    k = quadprog(H,f,A,b,Aeq,beq,lb,ub);
end

end