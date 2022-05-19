%% ��ģ��e2
clear all; close all; clc;

%% Current state load������״̬�����Ϣ
global R0 g
R0 = 10*10^3;                    %R0��λ��m
g = 9.81;

%% Constants������

n = 4;                      %����ѧ��5����
ccm_eps = 0.05;

%ȫ�������    
v_lim_l = 3;
v_lim_u = 7;
gama_lim_l = deg2rad(-75);           %�������Ƕȵĵ�λ���ǣ�����
gama_lim_u = deg2rad(-5);

state_lim = [v_lim_l,v_lim_u,gama_lim_l,gama_lim_u];

%% Uncomment for global optimization of bound�������OPT-CCM

%Ѱ����С�������Ĺ��̣���1.2�����������ٶ��ַ�����
%�Ż����յ���ʲô���ڲ�ͬ���������£�CCM���������ʱ������ָ���ǲ�ͬ�ģ���һϵ��lambda��������Ӧ������ָ�꣬�˹���ѡ���ŵ�lambda
%eucǰ����sqrt(double(w_upper/w_lower))/lambda���ȼ��ڣ�1/lambda^2��*(w_upper/w_lower)��������d_bar * diag(W_upper)^1/2
%d_barsǰ����sqrt(double(1/w_lower))/lambda��������alpha_w/lambda

lambda_range = linspace(0.1,2,10);              %���Բ�ֵ��lambda�ķ�Χ
lambda_range = (1/100)*round(lambda_range*100);   %��������Ϊ����
euc_bounds = NaN(length(lambda_range),1);         %ע�⣬ǰ���Ż�Ѱ��CCM�ͺ���׼ȷ����ʱ�ĺ��岻�ϸ���ͬ��
d_bars = NaN(length(lambda_range),1);             %bar������һ��
cond_bound = NaN(length(lambda_range),1);         %������

eps = 1;
condn_prev = 10;                                  %Ԥ��������
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
        [sos_prob,~,~] = Opt_Model_e2(n,state_lim,cond_u,lambda,ccm_eps,return_metric);
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
        if cond_u > 100
            break
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
        [sos_prob, w_lower, w_upper] = Opt_Model_e2(n,state_lim,condn,lambda,ccm_eps,return_metric);
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

    condn_prev = cond_u;              %����С��ȫ�����
    cond_bound(ll) = cond_u;

    disp('Euc_bound:'); disp(euc_bounds(ll));
    disp('d_bar:'); disp(d_bars(ll));
    fprintf('**********\n');
    
end


%% Pick a (lambda,condn)���������ȭ����һ�׵�

lambda = 0.7;          %������������
condn = 10;
return_metric = 1;

save_name = 'metric_Model_e2_vectorized.mat';          %֮��Ĳ���ֻ������һ��CCM_OTP������save_name�Ѷ���ȡ������֮ǰ����Ĳ����Ѿ�������Ż�CCM�Ĺ���
%%
[sos_prob, w_lower, w_upper] = Opt_Model_e2(n,state_lim,condn,lambda,ccm_eps,return_metric,save_name);      %����һ����ã�sos_prob:sum of squares problem

%% Compute aux control bound���������Ʊ߽�
%�����뿴�����������(lambda��condn)����£�CCM��ʲô���ı��֣���Ҫ����u���ܿص�ס
%���ƺ��Ŷ�����ϵģ�����һ���Ŷ��Ͻ磬�㻹�����������ʲô��tube�ڣ����������u�����
load(save_name);
disp('Checking CCM conditions and Computing control bound...');
%%
%----------------------BEGIN������ѧ��΢�ֶ���ѧ---------------------------%
% ����PVTOL������ǿ����þ�ȷ�ģ������䣬msspoly�Ǳ��ȶ���ʽ��
syms h v gama alpha          %�ݶ�Ϊ��һ��������д��һ���Ķ���ѧ

%dynamics f������ѧ
m  = 1600;        %����

CL = 0.4 * 180 / pi;
Cd = 0.02988 * CL^2;

f1 = -v * cos(gama);
f2 = v * sin(gama);
f3 = - Cd/(m*g)  - sin(gama);
f4 = - cos(gama)/v;

f_mat(v, gama) = [f1;f2;f3;f4];               %����ѧ��f���ڵ��������� ƫW/ƫf ��һ��

%����΢�ֶ���ѧ��df_mat��ǰ���У����ƻ�Ӱ�쵽��һ�У������õ���
df1_y = 0;
df1_h = 0;
df1_v = -cos(gama);
df1_gama = v * sin(gama);
df1 = [df1_y, df1_h, df1_v, df1_gama];

df2_y = 0;
df2_h = 0;
df2_v = sin(gama);
df2_gama = v * cos(gama);
df2 = [df2_y, df2_h, df2_v, df2_gama];

df3_y = 0;
df3_h = 0;
df3_v = 0;
df3_gama = -cos(gama);
df3 = [df3_y, df3_h, df3_v, df3_gama];

df4_y = 0;
df4_h = 0;
df4_v = cos(gama) / v^2;
df4_gama = sin(gama) / v;
df4 = [df4_y, df4_h, df4_v, df4_gama];

df_mat(v, gama) = [df1;df2;df3;df4];      %df���õ��ϣ�5*5ûë��

%------------------------END������ѧ��΢�ֶ���ѧ---------------------------%
%%
B_symfnc = [zeros(3,1);CL/(m*g*v)];       % 4*1 ��double

% B = double( B_symfnc(h1, v1, gama1) );
            
Bw_symfnc(v, gama) = [zeros(2,2);         %����ϵ���ϵ��Ŷ������Tube�ʹ������
                             0.01*Cd/(m*g), 0;
                             0, 0];     %�Ŷ���0.1��������ϵ���������㣬�����
% Bw = double( Bw_symfnc(h1, v1, gama1, alpha1) );
  
B_perp = [eye(4);zeros(1,4)];            %�����Ǹ�B��ֱ

%%
ctrl_N = 6;               %���Բ�ֵ�ĸ��������ǰ���ô�����״̬ȫģ����һ��
v_range = linspace(v_lim_l, v_lim_u, ctrl_N);
gama_range = linspace(gama_lim_l, gama_lim_u, ctrl_N);

delta_u = zeros(ctrl_N,ctrl_N);       %��2ά���飬��Ϊ��2��״̬
eig_CCM_min = zeros(ctrl_N, ctrl_N);
eig_CCM_max = zeros(ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,2);                     %Wֻ��i,j������Ϊ���u��Ӱ�칥�ǣ���������Wֻ�Ǹ߶ȡ��ٶȡ�������ǵĺ���
sigma_ThBw = zeros(ctrl_N,ctrl_N);
%%
for i = 1:length(v_range)
    fprintf(' i = %d , calculating \n',i);
    for j = 1:length(gama_range)
        fprintf(' j = %d , calculating \n',j);
                
                x = [0;0;v_range(i);gama_range(j)];    %������������״̬����
                
                Bw = double( Bw_symfnc(x(3),x(4)) );
                
                W = W_eval(w_poly_fnc(x));
                M = W\eye(n);               %����������
                Theta = chol(M);            %R = chol(A) ���ھ��� A �ĶԽ��ߺ������������������Ǿ��� R�����㷽�� R'*R=A��
                Theta_Bw = Theta*Bw;
                sigma_ThBw(i,j) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));   %Ϊ����ƫ���Ͻ��׼��
                
                L = chol(W);        %ע�⣬L_nor�͹�һ������������ͻ��������shift+enter�Ľ����ԭ������L
                
                f = double( f_mat(x(3),x(4)) );           %״̬΢��
                df = double( df_mat(x(3),x(4)) );         %״̬����΢��
                
                % W�ǶԳƵģ�FҲ�ǶԳƵģ�������G��x����F��Ҫ����������Ķ�Ӧ��ϵ���������⣿
                F = - W_eval(dw_poly_v_fnc(x)) * f(3) ...      
                            + df*W + W*df' + 2*lambda*W;  %����ƫ΢���Ѿ�������ȫ΢�֣�alpha��ȫ��0��ֱ�Ӳ���д
                
                delta_u_den = eig((inv(L))'*(B*B')*inv(L));     %L��W�ķֽ�
                delta_u(i,j) = 0.5*max(eig((inv(L))'*F*inv(L)))/...
                    sqrt(min( delta_u_den(delta_u_den>0) ));    %���ǿ�����
                
                R_CCM = -B_perp'*F*B_perp;  %ǰ�󶼳�B��ֱ��R_CCM��Ҫ����
                
                eig_CCM_min(i,j) = min(eig(R_CCM));     %����Ϊ��ֵ����������������
                eig_CCM_max(i,j) = max(eig(R_CCM));
                eig_W(i,1) = min(eig(W));
                eig_W(i,2) = max(eig(W));
                
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
disp(min(eig_W(:,1)));      %С��С
disp(max(eig_W(:,2)));      %���д�

disp('min eig CCM:'); 
disp(min(eig_CCM_min(:)));

disp('max eig CCM:'); 
disp(max(eig_CCM_max(:)));

disp('euc_bounds');
disp(d_bar*sqrt(diag(W_upper)));
