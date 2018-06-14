% Run a simple unit test on the matlab implementation for the various
% velocity IK solvers that are related to the original SNS-IK papers.
%

clc; clear;

% shared test parameters;
optTol = 1e-4;
cstTol = 1e-6;
fid = 1;
nTest = 100;

%% SNS-IK velocity solver, QP implementation
alpha = 1e-4;  % trade-off between minimum joint speed and max task scale
solver = @(dqLow, dqUpp, dx, J)( snsIk_vel_QP(dqLow, dqUpp, dx, J, alpha) );
result.QP = runTest_snsIk_vel(solver, nTest, optTol, cstTol, fid);

%% SNS-IK solver, Rethink Robotics revised algorithm (Andy Park)
result.rr = runTest_snsIk_vel(@snsIk_vel_rr, nTest, optTol, cstTol, fid);

%% SNS-IK solver, original SNS-IK algorithm
% % This solver currently does not pass the test
% result.basic = runTest_snsIk_vel(@snsIk_vel_basic, nTest, optTol, cstTol, fid);
