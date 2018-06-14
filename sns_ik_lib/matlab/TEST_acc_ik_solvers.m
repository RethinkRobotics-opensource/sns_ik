% This script sets up an runs a set of unit tests on the acceleration IK solvers.
%

clc; clear;

% shared test parameters;
optTol = 1e-4;
cstTol = 1e-6;
fid = 1;
nTest = 100;

%% SNS-IK acceleration solver, QP implementation
alpha = 1e-4;  % trade-off between minimum joint speed and max task scale
solver = @(dqLow, dqUpp, dx, J, dJdq)( snsIk_acc_QP(dqLow, dqUpp, dx, J, dJdq, alpha) );
result.QP = runTest_snsIk_acc(solver, nTest, optTol, cstTol, fid);
