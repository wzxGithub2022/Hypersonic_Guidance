%-------------------------------------------------------------------------%
%--------- BEGIN: function Hyper_dive_trajectory_2D_Endpoint.m -----------%
%-------------------------------------------------------------------------%
function output = Hyper_dive_trajectory_2D_Endpoint(input)

global impact_v R0 g0
K  = 1*10^(-1);                          %���KҪ�仯�� ��δ��һ��֮ǰ��ͬ��
vf = input.phase.finalstate(3);
output.objective = K*(impact_v/sqrt(R0*g0)-vf)^2 + input.phase.integral;  %cost��Ҳ������һ��������ƫ��

end
%-------------------------------------------------------------------------%
%---------- END: function Hyper_dive_trajectory_2D_Endpoint.m ------------%
%-------------------------------------------------------------------------%