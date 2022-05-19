%% CA一键运行，离线设计+在线跟踪框架

%% 离线设计，计算每个点的CCM
auto_test_all_CA;
% 每个点调用CCM_Hyper_Opt，生成CCW，CCM
% plot_CCM_state_CA;  %可以调用绘图看看计算效果

%% 在线跟踪，仿真
Hyper_dive_Trajectory_Simulation;
% ode时，即Hyper_dive_Dynamitics_2D_function，对CCM进行spline型插值，alpha_QP_function得到控制