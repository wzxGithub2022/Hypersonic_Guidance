%% B4_fix一键运行，离线设计+在线跟踪框架

%% 前期准备
% CCM_Dynamitics_no_chebysheff_V3_safe_B4；
% 生成df，B，D_nor，Lf_nor

%% 离线设计，计算每个点的CCM
auto_test_all_B4_fix;
% 每个点调用CCM_Hyper_Opt_linear，生成CCW，CCM
% plot_CCM_state_B4_fix;  %可以调用绘图看看计算效果

%% 在线跟踪，仿真
Hyper_dive_Trajectory_Simulation;
% ode时，即Hyper_dive_Dynamitics_2D_function，对df，B，D_nor，Lf_nor，CCM进行spline型插值，alpha_QP_function得到控制

%% 在仿真点计算CCM，比利用离线数据插值得到CCM，精度更高
%0426懂了，看论文，推一推，x_e才是正解
% 去控制0.01气动，lambda给0.9，condn给10，delta_gama = x_e
% 射程偏差km
%    -0.0111
% 高度偏差km
%    -0.0501