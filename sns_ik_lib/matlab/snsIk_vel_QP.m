function [dq, s, exitCode] = snsIk_vel_QP(dqLow, dqUpp, dxGoal, J, alpha)
% [dq, s, exitCode] = snsIk_vel_QP(dqLow, dqUpp, dxGoal, J, alpha)
%
% This function implements the basic QP version of the SNS-IK velocity
% solver.
%
% INPUTS:
%   dqLow: lower limit for joint velocities
%   dqUpp: upper limit for joint velocities
%   dxGoal: task-space velocity
%   J: jacoabian mapping joint space to task space
%   alpha: trade-off between minimum-joint velocity and maximum task scale
%      default: alpha = 1e-3
%      small alpha: favor solutions that maximize task scale
%      large alpha: favor minimum-joint velocity solutions
%
% OUTPUTS:
%   dq = joint velocity solution with maximum task scale factor
%   s = task scale factor [0, 1]
%
% NOTES:
%   s * dxGoal = J * dq;
%

% TODO: input validation
% TODO: return optimization status

if nargin < 5
   alpha = 1e-3;
end

% problem dimensions
[nTask, nJnt] = size(J);
nDecVar = nJnt + 1;  % [dq; -s]

% set the objective function
problem.H = blkdiag(alpha * eye(nJnt), 1.0 / alpha);
problem.f = zeros(nDecVar, 1);

% set the constraints:
problem.lb = [dqLow; 0];
problem.ub = [dqUpp; 1];
problem.Aeq = [J, dxGoal];
problem.beq = dxGoal;
problem.Aineq = [];
problem.bineq = [];

% solve the problem
problem.options = optimset('Display','off');
problem.solver = 'quadprog';
[zSoln, ~, exitCode] = quadprog(problem);
dq = zSoln(1:nJnt);
s = 1-zSoln(end);

end
