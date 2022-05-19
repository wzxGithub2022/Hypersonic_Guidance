%% ����һ���������ϰ������ã�����
% 2022/3/30 19:21

%% compute_CCM_Hyper
% ����������ȫû���õ�spotless�Ķ���������ֻ���Ż��õĹ���
clear all; close all; clc;
% yalmip('clear'); %��ͼ����Ĺ��ܣ����Ҳ���
warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Current state load������״̬�����Ϣ

load('Trajectory_all_information.mat');

i = 1;
h0 = Trajectory_all_information(i,3);
v0 = Trajectory_all_information(i,4);
gama0 = Trajectory_all_information(i,5);
alpha0 = Trajectory_all_information(i,6);
state_base = [h0,v0,gama0,alpha0];

%% Constants������

n = 5;                      %����ѧ��5����
g = 9.81;
ccm_eps = 0.05;

%ע�⣬����ÿһ��������
h_lim = 1;                  %�߶�h��λ��km
v_lim = 100;
gama_lim = pi/18;           %�������Ƕȵĵ�λ���ǣ�����
alpha_lim = pi/36;

%% Uncomment for global optimization of bound�������OPT-CCM

%Ѱ����С�������Ĺ��̣���1.2�����������ٶ��ַ�����
%�Ż����յ���ʲô���ڲ�ͬ���������£�CCM���������ʱ������ָ���ǲ�ͬ�ģ���һϵ��lambda��������Ӧ������ָ�꣬�˹���ѡ���ŵ�lambda
%eucǰ����sqrt(double(w_upper/w_lower))/lambda���ȼ��ڣ�1/lambda^2��*(w_upper/w_lower)��������d_bar * diag(W_upper)^1/2
%d_barsǰ����sqrt(double(1/w_lower))/lambda��������alpha_w/lambda

lambda_range = linspace(0.7,0.95,5);              %���Բ�ֵ��lambda�ķ�Χ
lambda_range = (1/100)*round(lambda_range*100);   %��������Ϊ����
euc_bounds = NaN(length(lambda_range),1);         %ע�⣬ǰ���Ż�Ѱ��CCM�ͺ���׼ȷ����ʱ�ĺ��岻�ϸ���ͬ��
d_bars = NaN(length(lambda_range),1);             %bar������һ��
cond_bound = NaN(length(lambda_range),1);         %������

eps = 1;
condn_prev = 50;                                  %Ԥ��������
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
        [sos_prob,~,~] = CCM_Hyper_Opt(n,g,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
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
        [sos_prob, w_lower, w_upper] = CCM_Hyper_Opt(n,g,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
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

    condn_prev = cond_u;              %����С��ȫ�����
    cond_bound(ll) = cond_u;

    disp('Euc_bound:'); disp(euc_bounds(ll));
    disp('d_bar:'); disp(d_bars(ll));
    fprintf('**********\n');
    
end

pause; 

%% Pick a (lambda,condn)���������ȭ����һ�׵�

lambda = 0.83;          %������������
condn = 132.8;
return_metric = 1;

save_name = 'metric_Hyper_vectorized.mat';          %֮��Ĳ���ֻ������һ��CCM_OTP������save_name�Ѷ���ȡ������֮ǰ����Ĳ����Ѿ�������Ż�CCM�Ĺ���
[sos_prob, w_lower, w_upper] = CCM_Hyper_Opt(n,g,h_lim,v_lim,gama_lim,alpha_lim,state_base,...
                                condn,lambda,ccm_eps,return_metric,save_name);      %����һ����ã�sos_prob:sum of squares problem

%% Compute aux control bound���������Ʊ߽�
%�����뿴�����������(lambda��condn)����£�CCM��ʲô���ı��֣���Ҫ����u���ܿص�ס
%���ƺ��Ŷ�����ϵģ�����һ���Ŷ��Ͻ磬�㻹�����������ʲô��tube�ڣ����������u�����

load(save_name);
disp('Checking CCM conditions and Computing control bound...');
lambda = 0.998*lambda;      %��Сһ����,�ⲻ�Ǻ�������

B = [zeros(4,1);1];       % 5*1 ��double
 
Bw = @(x)[zeros(2,2);
           -1/1600,0;
     0,1/(1600*2000);
                0,0];      %������д��ֱ�Ӿ���Сfunction�����������Ŷ���ϵͳ�����룬��СΪ5*2
  
B_perp = [eye(4);zeros(1,4)];            %�����Ǹ�B��ֱ

ctrl_N = 12;               %���Բ�ֵ�ĸ��������ǰ���ô�����״̬ȫģ����һ��
h_range = linspace(h0-h_lim, h0+h_lim, ctrl_N);      %ע�⣬����״̬��ֵӦ��ÿ������״̬��������[h0-h_lim,h0+hlim]
v_range = linspace(v0-v_lim, v0+v_lim, ctrl_N);
gama_range = linspace(gama0-gama_lim, gama0+gama_lim, ctrl_N);
alpha_range = linspace(alpha0-alpha_lim, alpha0+alpha_lim, ctrl_N);

%----------------------BEGIN������ѧ��΢�ֶ���ѧ---------------------------%

%handle function��������������������б�ѩ�����ʽ
sin_x_cheby = @(xx) 0.9101*(xx/(pi/3)) - 0.04466*(4*(xx/(pi/3))^3 - 3*(xx/(pi/3)));
cos_x_cheby = @(xx) 0.7441 -0.2499*(2*(xx/(pi/3))^2 -1);

sin_g = @(x) sin_x_cheby(x(4));        
cos_g = @(x) cos_x_cheby(x(4));        %gama
sin_a = @(x) sin_x_cheby(x(5));
cos_a = @(x) cos_x_cheby(x(5));        %alpha

%dynamics f������ѧ
S   = @(x) 0.5026;                       %�ο����
hf  = @(x) -x(2)/7110;                   %h fake���м����
rou = @(x) 1.225 * (1 + hf + hf^2/2 + hf^3/6 + hf^4/24 + hf^5/120 + hf^6/720 + hf^7/5040 + hf^8/40320 + hf^9/362880 + hf^10/3628800);
q   = @(x) 0.5 * rou * x(3)^2;           %��ѹ
qf  = @(x) 0.5 * rou * x(3);             %q fake��α��ѹ����Լ���ٶ�v
M   = @(x) x(3) / 340;                   %�����
m   = @(x) 1600;                         %����

CL  = @(x) 0.4172 + 19.41*x(5) + 10.17*x(5)^2 - M*(0.1004 + 0.7536*x(5));
L   = @(x) q*CL*S;                       %����
Lf  = @(x) qf*CL*S;                      %L fake��α��������Լ���ٶ�v

Cd  = @(x) 0.3042 + 0.02988*CL^2;
D   = @(x) q*Cd*S;                       %����

%�����б�ѩ��ƽ�1/V������[1000,2000]
V_division = @(x) 0.0028498843 - 2.984391e-06*x(3) + 1.3615968e-09*x(3)^2 - 2.285664e-13*x(3)^3;

f1 = @(x) -x(3) * cos_a;                 %���Ǻܱ�Ҫ��������������
f2 = @(x) x(3) * sin_a;
f3 = @(x) -D/m - g*sin_g;
% f4 = L/(m*x(3)) - g*cos_g/x(3);   %(msspoly�ĳ��������ǲ��������)
% f4 = Lf/m - g*cos_g/x(3);        %Lf��Լ���ٶ�v��g*cos_g/x(3)������
% ��δ������ȸ���һ�£�2022/03/27/10:10
f4 = @(x) Lf/m - g*cos_g*V_division;
f5 = @(x) 0;

f_mat = @(x) [f2;f3;f4;f5];               %����ѧ

%����΢�ֶ���ѧ��df_mat��ǰ���У����ƻ�Ӱ�쵽��һ�У������õ���
df11 = @(x) 0;
df12 = @(x) 0;
df13 = @(x) -cos_a;
df14 = @(x) 0;
df15 = @(x) x(3) * sin_a;
 df1 = @(x) [df11,df12,df13,df14,df15];

df21 = @(x) 0;
df22 = @(x) 0;
df23 = @(x) sin_a;
df24 = @(x) 0;
df25 = @(x) x(3) * cos_a;
 df2 = @(x) [df21,df22,df23,df24,df25];

df31 = @(x) 0;
df32 = @(x) diff(f3,x(2));
df33 = @(x) diff(f3,x(3));
df34 = @(x) diff(f3,x(4));
df35 = @(x) diff(f3,x(5));
 df3 = @(x) [df31,df32,df33,df34,df35];
 
df41 = @(x) 0;
df42 = @(x) diff(f4,x(2));
df43 = @(x) diff(f4,x(3));
df44 = @(x) diff(f4,x(4));
df45 = @(x) diff(f4,x(5));
 df4 = @(x) [df41,df42,df43,df44,df45];
 
df_mat = @(x) [df1;df2;df3;df4;zeros(1,5)];        %΢�ֶ���ѧ

%------------------------END������ѧ��΢�ֶ���ѧ---------------------------%         

delta_u = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);       %����ά���飬��Ϊ��4��״̬
eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,ctrl_N,ctrl_N,2);                     %Wֻ��i,j,k������Ϊ���u��Ӱ�칥�ǣ���������Wֻ�Ǹ߶ȡ��ٶȡ�������ǵĺ���
sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);

for i = 1:length(h_range)
    for j = 1:length(v_range)
        for k = 1:length(gama_range)
            for l = 1:length(alpha_range)      %����Сд��L��
                
                x = [randn(1,1);h_range(i);v_range(j);gama_range(k);alpha_range(l)];    %������������״̬����
                
                W = W_eval(w_poly_fnc(x));
                M = W\eye(n);               %����������
                Theta = chol(M);            %R = chol(A) ���ھ��� A �ĶԽ��ߺ������������������Ǿ��� R�����㷽�� R'*R=A��
                Theta_Bw = Theta*Bw(x);
                sigma_ThBw(i,j,k,l) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));   %Ϊ����ƫ���Ͻ��׼��
                
                L = chol(W);
                
                f = f_mat(x);           %״̬΢��
                df = df_mat(x);         %״̬����΢��
                
                % W�ǶԳƵģ�FҲ�ǶԳƵģ�������G��x��
                F = -W_eval(dw_poly_h_fnc(x)) * f(2) - W_eval(dw_poly_v_fnc(x)) * f(3) - W_eval(dW_gama_fnc(x)) * f(4)...      
                            + df*W + W*df' + 2*lambda*W;  %����ƫ΢���Ѿ�������ȫ΢�֣�alpha��ȫ��0��ֱ�Ӳ���д
                
                delta_u_den = eig((inv(L))'*(B*B')*inv(L));     %L��W�ķֽ�
                delta_u(i,j,k,l) = 0.5*max(eig((inv(L))'*F*inv(L)))/...
                    sqrt(min( delta_u_den(delta_u_den>0) ));    %���ǿ�����
                
                R_CCM = -B_perp'*F*B_perp;  %ǰ�󶼳�B��ֱ
                
                eig_CCM(i,j,k,l) = min(eig(R_CCM));
                eig_W(i,j,1) = min(eig(W));
                eig_W(i,j,2) = max(eig(W));
                
            end
        end
    end
end

alpha_w = max(sigma_ThBw(:));       %(:)��ʾ��������Ԫ�ز����б�ʾ
d_bar = alpha_w/lambda;             %��d_barû�г�w�Ͻ磬������J_CCMҲ�����ٳ�w�Ͻ��ˣ�һ����λ
disp('d_bar'); 
disp(d_bar);

disp('Control:'); 
disp(max(d_bar*delta_u(:)));        %���ǿ����Ͻ�

disp('W:'); 
disp(min(min(eig_W(:,:,1))));      %С��С
disp(max(max(eig_W(:,:,2))));      %���д�

disp('CCM:'); 
disp(min(eig_CCM(:)));

disp('euc_bounds');
disp(d_bar*sqrt(diag(W_upper)));
