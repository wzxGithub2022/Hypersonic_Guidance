%-------------------------------------------------------------------------%
%--------- BEGIN: function Hyper_dive_trajectory_2D_Endpoint.m -----------%
%-------------------------------------------------------------------------%
function output = Hyper_dive_trajectory_2D_Endpoint(input)

global impact_v R0 g0
K  = 1*10^(-1);                          %这个K要变化了 与未归一化之前不同了
vf = input.phase.finalstate(3);
output.objective = K*(impact_v/sqrt(R0*g0)-vf)^2 + input.phase.integral;  %cost，也引入了一部分落速偏差

end
%-------------------------------------------------------------------------%
%---------- END: function Hyper_dive_trajectory_2D_Endpoint.m ------------%
%-------------------------------------------------------------------------%