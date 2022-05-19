%为了提高程序的可比性和运行效率，基于在线算法框架，对B4_fix进行重构
function [dx ,alpha] = B4yyds_Dynamics(t,X)
global t_r u_r x_r y_r v_r gamma_r feedback w disturb
global R0 g
t
fprintf('time now is %f second \n',t*sqrt(R0/g));

n=4;
lambda = 0.9;
condn = 10;
ccm_eps = 0.05;
return_metric = 1;

xd = interp1(t_r , x_r ,t);     %按默认方法是线性插值来，线性化的大旗
yd = interp1(t_r , y_r ,t);
vd = interp1(t_r , v_r ,t);
gammad = interp1(t_r , gamma_r ,t);
alphad = interp1(t_r ,u_r, t);

%对每个插值点都保持了相当精度的线性化
df =[ 0,                                                                                                                                                                       0,                                                                                                                                                                                                                                                                                                                                                         -cos(gammad), vd*sin(gammad),                                                                                                                                                                                                          0;
      0,                                                                                                                                                                       0,                                                                                                                                                                                                                                                                                                                                                          sin(gammad), vd*cos(gammad),                                                                                                                                                                                                          0;
      0, (123137*vd^2*exp(-(1000*yd)/711)*((747*((1941*alphad)/100 + (1017*alphad^2)/100 - (3*109^(1/2)*vd*((471*alphad)/625 + 251/2500))/34 + 1043/2500)^2)/25000 + 1521/5000))/45504, (275950017*109^(1/2)*vd^2*exp(-(1000*yd)/711)*((471*alphad)/625 + 251/2500)*((1941*alphad)/100 + (1017*alphad^2)/100 - (3*109^(1/2)*vd*((471*alphad)/625 + 251/2500))/34 + 1043/2500))/27200000000 - (123137*vd*exp(-(1000*yd)/711)*((747*((1941*alphad)/100 + (1017*alphad^2)/100 - (3*109^(1/2)*vd*((471*alphad)/625 + 251/2500))/34 + 1043/2500)^2)/25000 + 1521/5000))/32000,  -cos(gammad), -(91983339*vd^2*exp(-(1000*yd)/711)*((1017*alphad)/50 - (1413*109^(1/2)*vd)/21250 + 1941/100)*((1941*alphad)/100 + (1017*alphad^2)/100 - (3*109^(1/2)*vd*((471*alphad)/625 + 251/2500))/34 + 1043/2500))/800000000;
      0,                              -(123137*vd*exp(-(1000*yd)/711)*((1941*alphad)/100 + (1017*alphad^2)/100 - (3*109^(1/2)*vd*((471*alphad)/625 + 251/2500))/34 + 1043/2500))/45504,                                                                                                                             (123137*exp(-(1000*yd)/711)*((1941*alphad)/100 + (1017*alphad^2)/100 - (3*109^(1/2)*vd*((471*alphad)/625 + 251/2500))/34 + 1043/2500))/64000 + cos(gammad)/vd^2 - (369411*109^(1/2)*vd*exp(-(1000*yd)/711)*((471*alphad)/625 + 251/2500))/2176000, sin(gammad)/vd,                                                                                                                (123137*vd*exp(-(1000*yd)/711)*((1017*alphad)/50 - (1413*109^(1/2)*vd)/21250 + 1941/100))/64000];
A = df(:,1:4);
b = df(:,5);

%对每个插值点算CCM，而不是先算出来46个CCM，之后再对CCM插值，我确实固化在离线计算+在线应用的框架下
X_d = [xd;yd;vd;gammad];
if feedback
    u_d = 0;     %标称的反馈，没有的，直接就是让攻角程序化地改变
    save_name = 'Hyper_B4_Wx_rule.mat';
    [solved,w_lower,w_upper,W_upper] = B4yyds_CCM_Opt(n,A,b,condn,lambda,ccm_eps,return_metric,save_name);
    load(save_name);
    W = W_eval(w_poly_fnc(X-X_d));  %由当前状态X执行W的生成规则得到W
    [k] = B4yyds_QP_Controllor(X_d, u_d, X, lambda, W, A, b);
    alpha = alphad+k;
else
    alpha = alphad;
end

x = X(1);
y = X(2);
v = X(3);
gamma = X(4);

R0 = 10000;
g0 = 9.81;
S = 0.5026;

rou=1.225.*exp(-y.*R0./7110);  %密度rou
q=0.5.*rou.*(v.*sqrt(R0.*g0)).^2;           %动压    
M=v.*sqrt(R0.*g0)./340;    %马赫数
m=1600; %质量
m0=m;
m1=m/m0;
CL=0.4172+19.41*alpha+10.17.*alpha.^2-M.*(0.1004+0.7536.*alpha);   
L1=q.*CL.*S/(m0*g0); %升力
Cd0=0.3042;
Cd=Cd0+0.02988.*CL.^2;           %阻力系数
D1=q.*Cd.*S/(m0*g0); %阻力

if disturb
    D1_dtb = D1./m1 * 0.01;
    Lf1_dtb = L1./v./m1 * 0.01;
else
    D1_dtb = 0;
    Lf1_dtb = 0;
end

dx = [-v.*cos(gamma);
      v.*sin(gamma);
      -D1./m1 - sin(gamma) + D1_dtb;
      L1./v./m1 - cos(gamma)./v + Lf1_dtb];

end