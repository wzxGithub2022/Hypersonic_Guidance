%-------------------------------------------------------------------------%
%--------- BEGIN: function Hyper_dive_trajectory_2D_Continuous.m ---------%
%-------------------------------------------------------------------------%
function phaseout = Hyper_dive_trajectory_2D_Continuous(input)

g0 = input.auxdata.g0;
S  = input.auxdata.S;
R0 = input.auxdata.R0;

t1 = input.phase.time;
y1 = input.phase.state(:,1);
z1 = input.phase.state(:,2);
v1 = input.phase.state(:,3);
gama = input.phase.state(:,4);
u = input.phase.control;            %控制量就是攻角alpha

rou = 1.225.*exp(-z1.*R0./7110);                 %密度rou
q   = 0.5.*rou.*(v1.*sqrt(R0.*g0)).^2;           %动压    
M   = v1.*sqrt(R0.*g0)./340;                     %马赫数
m   = 1600;                                      %质量
m0  = m;                                         %质量的归一化可以说很寂寞
m1  = m/m0;                                      %这m1恒为1，给m也归一化了
CL  = 0.4172+19.41*u+10.17.*u.^2-M.*(0.1004+0.7536.*u);   
L1  = q.*CL.*S/(m0*g0);                          %升力
Cd0 = 0.3042;
Cd  = Cd0+0.02988.*CL.^2;                        %阻力系数
D1  = q.*Cd.*S/(m0*g0);                          %阻力

y1dot =  -v1.*cos(gama) ;
z1dot =  v1.*sin(gama);
v1dot =  -D1./m1-sin(gama);
gamadot = L1./(m1.*v1)-cos(gama)./v1;
phaseout.dynamics = [y1dot, z1dot, v1dot, gamadot];     %给出微分方程就行，系统自己会把微分约束转化为代数约束
phaseout.integrand = u.^2;                              %注意这个积分值小于实际，原因是：时间归一化
%注：u就是攻角，控制代价为phaseout.integrand = alpha.^2
end
%-------------------------------------------------------------------------%
%---------- END: function Hyper_dive_trajectory_2D_Continuous.m ----------%
%-------------------------------------------------------------------------%