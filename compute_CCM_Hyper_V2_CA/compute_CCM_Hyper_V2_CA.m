%% compute_CCM_Hyper
% �汾������˵����V2 CA control_alpha 
% ���׵�alpha bug���޸�
% ��ǰ�Թ��Ǳ仯����Ϊ��������˼·�ǿ�ѧ�ģ��������ϰ��Ҹ�(0.7,10)ʱG(x)���������������������Ͽ�������
% ����ֹͣ״̬��չ���������ù��ǣ�B��B_perp����Ҫ��д������V2 CA
% Ԥ��V2 CA���ٶȻ��V1 CAD�죬��Ч������CAD��

% clear all; close all; clc;
% yalmip('clear'); %��ͼ����Ĺ��ܣ����Ҳ���
warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Current state load������״̬�����Ϣ

load('Trajectory_normalization.mat');

global state_num state_CCM CCW_upper
% state_num = 1;  %������һ������Ҫ����һ�зų���

i = state_num;
h1 = Trajectory_normalization(i,3);
v1 = Trajectory_normalization(i,4);
gama1 = Trajectory_normalization(i,5);
alpha1 = Trajectory_normalization(i,6);     %�����벻��alpha1
state_base = [h1,v1,gama1,alpha1];          %alpha1����ȥֱ��������

global R0 g
R0 = 10*10^3;                    %R0��λ��m
g = 9.81;

%% Constants������

n = 4;                      %����ѧ��4����
ccm_eps = 0.05;             %���ֵ0.05�������ǲ���̫�Ͽ���

%ע�⣬����ÿһ��������
h_lim = 100/R0;                  
v_lim = 30/sqrt(R0*g);
gama_lim = 1*pi/180;           %�������Ƕȵĵ�λ���ǣ�����

%% Uncomment for global optimization of bound�������OPT-CCM
% 
%Ѱ����С�������Ĺ��̣���1.2�����������ٶ��ַ�����
%�Ż����յ���ʲô���ڲ�ͬ���������£�CCM���������ʱ������ָ���ǲ�ͬ�ģ���һϵ��lambda��������Ӧ������ָ�꣬�˹���ѡ���ŵ�lambda
%eucǰ����sqrt(double(w_upper/w_lower))/lambda���ȼ��ڣ�1/lambda^2��*(w_upper/w_lower)��������d_bar * diag(W_upper)^1/2
%d_barsǰ����sqrt(double(1/w_lower))/lambda��������alpha_w/lambda

% lambda_range = linspace(0.1,2,10);              %���Բ�ֵ��lambda�ķ�Χ
lambda_range = 0.7;
lambda_range = (1/100)*round(lambda_range*100);   %��������Ϊ����
euc_bounds = NaN(length(lambda_range),1);         %ע�⣬ǰ���Ż�Ѱ��CCM�ͺ���׼ȷ����ʱ�ĺ��岻�ϸ���ͬ��
d_bars = NaN(length(lambda_range),1);             %bar������һ��
cond_bound = NaN(length(lambda_range),1);         %������

eps = 1;
condn_prev = 5;                                  %Ԥ��������
return_metric = 0;

for ll = 1:length(lambda_range)       %��ÿһ��lambda��������Сд�ġ�LL��,���ۣ�������
    
    lambda = lambda_range(ll);
    fprintf('**********\n');
    fprintf('lambda: %f\n', lambda);
    solved = 0;
    
    %Determine upper bound�������Ͻ�
    cond_l = condn_prev;              
    cond_u = 1.2*condn_prev;
%------------------------������һֱ��ѭ�����ǿ���������䣬sos����һֱ�ò������-----------------------%
    while (~solved) 
        fprintf(' cond_u: %.2f: ', cond_u);
        [sos_prob,~,~] = CCM_Hyper_Opt(n,h_lim,v_lim,gama_lim,state_base,...
                                cond_u,lambda,ccm_eps,return_metric);
        if (sos_prob == 0)
            solved = 1;
            fprintf('feasible \n');
        else
            %shift up condition number range
            %������˵Ļ���һ������condition number��ֻ����һ���������
            fprintf('�㿴��ⲻ�����˰�\n');
            cond_l = cond_u;
            cond_u = 1.2*cond_u;
        end
    end
%---------------------------------------------------------------------------------------------------%
    if (solved)           %ֱ������ˣ���euc_bounds(ll)������
        euc_bounds(ll) = sqrt(cond_u)/lambda;       %���������һ��ֵ����ʱ���޷�����d_bars
        fprintf(' cond_l: %.2f, cond_u: %.2f\n', cond_l, cond_u);
    else
        continue;
    end
    
    %Now do bisection search������������
    while(cond_u - cond_l >= eps)     %���ϴ�ϸ���Ż�
        condn = (cond_l+cond_u)/2;
        fprintf(' cond: %.2f', condn);
        [sos_prob, w_lower, w_upper] = CCM_Hyper_Opt(n,h_lim,v_lim,gama_lim,state_base,...
                                condn,lambda,ccm_eps,return_metric);
        if (sos_prob == 0)
            fprintf(' feasible\n');       %��Ϳ����ˣ�������ȡСһ�㣬���󶯵Ļ���euc��d�ͻ������Ż�һЩ      
            euc_bounds(ll) = sqrt(double(w_upper/w_lower))/lambda;   %euc�ȼ��ڣ�1/lambda^2��*(w_upper/w_lower)
            d_bars(ll) = sqrt(double(1/w_lower))/lambda;             %������״μ���d_bars
            cond_u = condn;
        else
            fprintf(' infeasible\n');     %�����еĻ������ҷŴ�һ��
            cond_l = condn;
        end
    end

%     condn_prev = cond_u;            %���д�������ķ��˼·��lambda����Ԥ��������Ҳ����
    cond_bound(ll) = cond_u;      %����С��ȫ�����

    disp('Euc_bound:'); disp(euc_bounds(ll));
    disp('d_bar:'); disp(d_bars(ll));
    fprintf('**********\n');
    
end
% figure(state_num)
% subplot(211),plot(lambda_range,euc_bounds);title('Euc_bound');
% subplot(212),plot(lambda_range,d_bars);title('d_bars');
% pause; 

%% Pick a (lambda,condn)���������ȭ����һ�׵�

lambda = 0.7;          %������������
condn = cond_u;
return_metric = 1;

save_name = 'metric_Hyper_vectorized.mat';          %֮��Ĳ���ֻ������һ��CCM_OTP������save_name�Ѷ���ȡ������֮ǰ����Ĳ����Ѿ�������Ż�CCM�Ĺ���
%%
[sos_prob, w_lower, w_upper] = CCM_Hyper_Opt(n,h_lim,v_lim,gama_lim,state_base,...
                                condn,lambda,ccm_eps,return_metric,save_name);      %����һ����ã�sos_prob:sum of squares problem

%% Compute aux control bound���������Ʊ߽�
%�����뿴�����������(lambda��condn)����£�CCM��ʲô���ı��֣���Ҫ����u���ܿص�ס
%���ƺ��Ŷ�����ϵģ�����һ���Ŷ��Ͻ磬�㻹�����������ʲô��tube�ڣ����������u�����

load(save_name);
disp('Checking CCM conditions and Computing control bound...');
%% �ƽ��� ����ѧ��΢�ֶ���ѧ
%----------------------BEGIN������ѧ��΢�ֶ���ѧ---------------------------%
% ����PVTOL������ǿ����þ�ȷ�ģ������䣬msspoly�Ǳ��ȶ���ʽ��
syms h v gama alpha           %�ݶ�Ϊ��һ��������д��һ���Ķ���ѧ

%dynamics f������ѧ
S   = 0.5026;                       %�ο����
% rou = 1.225 * exp(-x(2)/7110);    %�ܶ�rou����Ҫ����ʽ�ƽ�
hf  = -h*R0/7110;                    %h fake���м������exp_cheby_bug�Ѹı�ע
rou = 1.225 * (0.996838143712596 + 0.957192272404239*hf + 0.403293676867969*hf^2 + 0.083714145730035*hf^3 + 0.006865365386321*hf^4);
q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %��ѹ
qf  = 0.5 * rou * v*(R0*g);               %q fake��α��ѹ����Լ���ٶ�v��bug��ע
M   = v*sqrt(R0*g) / 340;                     %�����
m   = 1600;                                   %����

CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
L_nor = q*CL*S / (m*g);                        %����
Lf_nor = qf*CL*S / (m*g);                      %L fake��α��������Լ���ٶ�v

Cd  = 0.3042 + 0.02988*CL^2;
D_nor = q*Cd*S / (m*g);                        %����

%handle function����������������б�ѩ�����ʽ���бƽ�
sin_x = @(xx) 0.985400513847416*xx - 0.142496853019355*xx^3;
cos_x = @(xx) 0.999396553656190 - 0.495559134405119*xx^2 + 0.036782872656059*xx^4;

sin_g = sin_x(gama);        
cos_g = cos_x(gama);         %gama

%�б�ѩ��ƽ�1/V���Թ�һ��֮���V������[3,7]
V_division = 1.101028089323823 - 0.473922941895471*v + 0.099720451270527*v^2 - 0.010265506659430*v^3 + 4.140771839924636e-04*v^4;

f1 = -v * cos_g;
f2 = v * sin_g;
f3 = -D_nor - sin_g;
f4 = Lf_nor - cos_g*V_division;

%���㿴͸�������f_mat�ͺ��ڵ�f(1...4)�Ĺ�ϵ����ȫ���ÿ�������compute��opt��д������
f_mat(h, v, gama, alpha) = [f1;f2;f3;f4];               %����ѧ��f���ڵ��������� ƫW/ƫf ��һ��

%΢�ֶ���ѧ��df_mat
df1_y = 0;
df1_h = 0;
df1_v = -cos_g;
df1_gama = v * sin_g;

df1 = [df1_y, df1_h, df1_v, df1_gama];

df2_y = 0;
df2_h = 0;
df2_v = sin_g;
df2_gama = v * cos_g;

df2 = [df2_y, df2_h, df2_v, df2_gama];

df3_y = 0;
df3_h = diff(f3,h);
df3_v = diff(f3,v);
df3_gama = -cos_g;

df3 = [df3_y, df3_h, df3_v, df3_gama];

df4_y = 0;
df4_h = diff(f4,h);
df4_v = diff(f4,v);
df4_gama = sin_g*V_division;

df4 = [df4_y, df4_h, df4_v, df4_gama];

df_mat(h, v, gama, alpha) = [df1;df2;df3;df4];      %4*4����΢��ϵ������
%------------------------END������ѧ��΢�ֶ���ѧ---------------------------%
% d_bar
%    15.9294
% Control:
%    4.6769e+08
% W:
%     1.0000
%     3.7252
% min eig CCM:
%     0.0558
% max eig CCM:
%     1.3816
% euc_bounds
%    19.9095
%    19.6997
%    28.4844
%    15.9297
%% ��ȷ�� ����ѧ��΢�ֶ���ѧ
%----------------------BEGIN������ѧ��΢�ֶ���ѧ---------------------------%
% ����PVTOL������ǿ����þ�ȷ�ģ������䣬msspoly�Ǳ��ȶ���ʽ��
% syms h v gama alpha           %�ݶ�Ϊ��һ��������д��һ���Ķ���ѧ
% 
% %dynamics f������ѧ
% S   = 0.5026;                       %�ο����
% rou = 1.225 * exp(-h*R0/7110);      %�ܶ�rou�����Բ��ö���ʽ�ƽ�
% q   = 0.5 * rou * (v*sqrt(R0*g))^2;           %��ѹ
% qf  = 0.5 * rou * v*(R0*g);               %q fake��α��ѹ����Լ���ٶ�v��bug��ע
% M   = v*sqrt(R0*g) / 340;                     %�����
% m   = 1600;                                   %����
% 
% CL  = 0.4172 + 19.41*alpha + 10.17*alpha^2 - M*(0.1004 + 0.7536*alpha);
% L_nor = q*CL*S / (m*g);                        %����
% Lf_nor = qf*CL*S / (m*g);                      %L fake��α��������Լ���ٶ�v
% 
% Cd  = 0.3042 + 0.02988*CL^2;
% D_nor = q*Cd*S / (m*g);                        %����
% 
% f1 = -v * cos(gama);
% f2 = v * sin(gama);
% f3 = -D_nor - sin(gama);
% f4 = Lf_nor - cos(gama)/v;
% 
% %���㿴͸�������f_mat�ͺ��ڵ�f(1...4)�Ĺ�ϵ����ȫ���ÿ�������compute��opt��д������
% f_mat(h, v, gama, alpha) = [f1;f2;f3;f4];               %����ѧ��f���ڵ��������� ƫW/ƫf ��һ��
% 
% %΢�ֶ���ѧ��df_mat
% df1_y = 0;
% df1_h = 0;
% df1_v = diff(f1,v);
% df1_gama = diff(f1,gama);
% 
% df1 = [df1_y, df1_h, df1_v, df1_gama];
% 
% df2_y = 0;
% df2_h = 0;
% df2_v = diff(f2,v);
% df2_gama = diff(f2,gama);
% 
% df2 = [df2_y, df2_h, df2_v, df2_gama];
% 
% df3_y = 0;
% df3_h = diff(f3,h);
% df3_v = diff(f3,v);
% df3_gama = diff(f3,gama);
% 
% df3 = [df3_y, df3_h, df3_v, df3_gama];
% 
% df4_y = 0;
% df4_h = diff(f4,h);
% df4_v = diff(f4,v);
% df4_gama = diff(f4,gama);
% 
% df4 = [df4_y, df4_h, df4_v, df4_gama];
% 
% df_mat(h, v, gama, alpha) = [df1;df2;df3;df4];      %4*4����΢��ϵ������
%------------------------END������ѧ��΢�ֶ���ѧ---------------------------%
%�����ף���д����d_bar�Ķ���ѧ��������2-3��������
%����exp_cheby�Ĵ�����ô���ûʲô�仯���ⲻ��ѧ
% d_bar
%     0.0116
% Control:
%    1.8822e+06
% W:
%     1.0000
%     3.7252
% min eig CCM:
%     0.0582
% max eig CCM:
%     1.4243
% euc_bounds
%     0.0145
%     0.0143
%     0.0207
%     0.0116
%%
%Bֻ��������Ͻ������ù������Կ����ӵ�alpha_dot�����Ǿ��ǿ���������ֱ�Ӹ�ǿ�����Ŀ��ƿ϶��Ǹ��ÿص�
%����ʹ��alpha1�������Ƿ�����߼��ϵ����⣬Ǩ�ͣ�������ȫ���Ի�Ҫǿ
% B��BwӦ�ö�����״̬�ı���ı��

B_symfnc(h, v, gama, alpha) = [0;0;diff(f3,alpha);diff(f4,alpha)];
% B = double( B_symfnc(h1, v1, gama1, alpha1) );      % 4*1 ��double

Bw_symfnc(h, v, gama, alpha) = [zeros(2,2);         %����ϵ���ϵ��Ŷ������Tube�ʹ������
                                0.01*D_nor, 0;
                                0, 0.01*Lf_nor];     %�Ŷ���0.1��������ϵ�����
% Bw = double( Bw_symfnc(h1, v1, gama1, alpha1) );

B_perp = [eye(2);zeros(2,2)];            %�����Ǹ�B��ֱ
%%
ctrl_N = 6;               %���Բ�ֵ�ĸ��������ǰ���ô�����״̬ȫģ����һ��
h_range = linspace(h1-h_lim, h1+h_lim, ctrl_N);      %ע�⣬����״̬��ֵӦ��ÿ������״̬��������[h0-h_lim,h0+hlim]
v_range = linspace(v1-v_lim, v1+v_lim, ctrl_N);
gama_range = linspace(gama1-gama_lim, gama1+gama_lim, ctrl_N);

delta_u = zeros(ctrl_N,ctrl_N,ctrl_N);       %��3ά���飬��ΪҪ����3��״̬
eig_CCM_min = zeros(ctrl_N, ctrl_N, ctrl_N);
eig_CCM_max = zeros(ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,2);                     %Wֻ��i������Ϊ��Ĺ��ǻ�Ӱ��v��gama����������Wֻ�Ǹ߶�h�ĺ���
sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N);
%%
for i = 1:length(h_range)
    fprintf(' i = %d , calculating \n',i);
    for j = 1:length(v_range)
        fprintf(' j = %d , calculating \n',j);
        for k = 1:length(gama_range)
                
                x = [randn(1,1);h_range(i);v_range(j);gama_range(k)];    %������������״̬����
                
                B = double( B_symfnc(x(2),x(3),x(4), alpha1) );      
                Bw = double( Bw_symfnc(x(2),x(3),x(4), alpha1) );
                
                W = W_eval(w_poly_fnc(x));  %��״̬xִ��W�����ɹ���õ�W
                M = W\eye(n);               %����������
                Theta = chol(M);            %R = chol(A) ���ھ��� A �ĶԽ��ߺ������������������Ǿ��� R�����㷽�� R'*R=A��
                Theta_Bw = Theta*Bw;
                sigma_ThBw(i,j,k) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));   %Ϊ����ƫ���Ͻ��׼��
                
                L = chol(W);        %ע�⣬L_nor�͹�һ������������ͻ��������shift+enter�Ľ����ԭ������L
                
                %alpha1���벻����
                %����ʹ��alpha1�������Ƿ�����߼��ϵ����⣬Ǩ�ͣ�������ȫ���Ի�Ҫǿ
                f = double( f_mat(x(2),x(3),x(4), alpha1) );           %״̬΢��
                df = double( df_mat(x(2),x(3),x(4), alpha1) );         %״̬����΢��
                
                % W�ǶԳƵģ�FҲ�ǶԳƵģ�������G��x����F��Ҫ����������Ķ�Ӧ��ϵ�й�bug���ȸ�Ϊ��
                F = -W_eval(dw_poly_h_fnc(x)) * f(2) + df*W + W*df' + 2*lambda*W;  %B4��ȫ΢�־�ֻ��һ��ƫ΢����
                
                delta_u_den = eig((inv(L))'*(B*B')*inv(L));     %L��W�ķֽ�
                delta_u(i,j,k) = 0.5*max(eig((inv(L))'*F*inv(L)))/...
                    sqrt(min( delta_u_den(delta_u_den>0) ));    %���ǿ�����
                
                R_CCM = -B_perp'*F*B_perp;  %ǰ�󶼳�B��ֱ��R_CCM��Ҫ����
                
                eig_CCM_min(i,j,k) = min(eig(R_CCM));
                eig_CCM_max(i,j,k) = max(eig(R_CCM));
                eig_W(i,1) = min(eig(W));
                eig_W(i,2) = max(eig(W));
                
        end
    end
end
%%
alpha_w = max(sigma_ThBw(:));       %(:)��ʾ��������Ԫ�ز����б�ʾ
d_bar = alpha_w/lambda;             %��d_barû�г�w�Ͻ磬������J_CCMҲ�����ٳ�w�Ͻ��ˣ�һ����λ
disp('d_bar'); 
disp(d_bar);

disp('Control:'); 
disp(max(d_bar*delta_u(:)));        %���ǿ����Ͻ�

disp('W:'); 
disp(min( eig_W(:,1) ));      %С��С
disp(max( eig_W(:,2) ));      %���д�

disp('min eig CCM:'); 
disp(min(eig_CCM_min(:)));

disp('max eig CCM:'); 
disp(max(eig_CCM_max(:)));

disp('euc_bounds');
disp(d_bar*sqrt(diag(W_upper)));

%% auto test all ����
state_CCM(1,state_num) = alpha_w/lambda;
state_CCM(2,state_num) = max(d_bar*delta_u(:));
state_CCM(3,state_num) = min(eig_W(:,1));
state_CCM(4,state_num) = max(eig_W(:,2));
state_CCM(5,state_num) = min(eig_CCM_min(:));
state_CCM(6,state_num) = max(eig_CCM_max(:));
state_CCM(7:10,state_num) = d_bar*sqrt(diag(W_upper));
CCW_upper(:,:,state_num) = W_upper;