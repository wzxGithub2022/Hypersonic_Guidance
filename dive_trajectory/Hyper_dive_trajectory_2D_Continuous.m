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
u = input.phase.control;            %���������ǹ���alpha

rou = 1.225.*exp(-z1.*R0./7110);                 %�ܶ�rou
q   = 0.5.*rou.*(v1.*sqrt(R0.*g0)).^2;           %��ѹ    
M   = v1.*sqrt(R0.*g0)./340;                     %�����
m   = 1600;                                      %����
m0  = m;                                         %�����Ĺ�һ������˵�ܼ�į
m1  = m/m0;                                      %��m1��Ϊ1����mҲ��һ����
CL  = 0.4172+19.41*u+10.17.*u.^2-M.*(0.1004+0.7536.*u);   
L1  = q.*CL.*S/(m0*g0);                          %����
Cd0 = 0.3042;
Cd  = Cd0+0.02988.*CL.^2;                        %����ϵ��
D1  = q.*Cd.*S/(m0*g0);                          %����

y1dot =  -v1.*cos(gama) ;
z1dot =  v1.*sin(gama);
v1dot =  -D1./m1-sin(gama);
gamadot = L1./(m1.*v1)-cos(gama)./v1;
phaseout.dynamics = [y1dot, z1dot, v1dot, gamadot];     %����΢�ַ��̾��У�ϵͳ�Լ����΢��Լ��ת��Ϊ����Լ��
phaseout.integrand = u.^2;                              %ע���������ֵС��ʵ�ʣ�ԭ���ǣ�ʱ���һ��
%ע��u���ǹ��ǣ����ƴ���Ϊphaseout.integrand = alpha.^2
end
%-------------------------------------------------------------------------%
%---------- END: function Hyper_dive_trajectory_2D_Continuous.m ----------%
%-------------------------------------------------------------------------%