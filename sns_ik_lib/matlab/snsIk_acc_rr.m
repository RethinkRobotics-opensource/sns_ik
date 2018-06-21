function [ddq, s, exitCode] = snsIk_acc_rr(ddqLow, ddqUpp, ddxGoal, J, dJdq)
% [ddq, s, exitCode] = snsIk_acc_rr(ddqLow, ddqUpp, ddxGoal, J, dJdq)
%
% This function implements the basic version of the SNS-IK acceleration
% solver.
%
% INPUTS:
%   ddqLow (nJnt x 1): lower limit for joint acceleration
%   ddqUpp (nJnt x 1): upper limit for joint acceleration
%   ddxGoal (ndx x 1): task-space acceleration
%   J (ndx x nJnt): jacoabian mapping joint space to task space
%   dJdq (ndx x 1) = dJ * dq
%       dJ = time-derivative of the task jacobian
%       dq = current joint velocity
%
% OUTPUTS:
%   ddq (nJnt x 1): joint acceleration solution with maximum task scale factor
%   s: task scale factor [0, 1]
%   exitCode    (1 == success)
%
%
% NOTES:
%   s * ddxGoal = J * ddq + dJ * dq;
%
% This function is a modification of the basic SNS-IK algorithm that is
% presented in the original SNS-IK papers. It was developed by Andy Park at
% Rethink Robotics in June 2018.

% Copyright 2018 Rethink Robotics
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% TODO: input validation

[nTask, nJnt] = size(J);

W = eye(nJnt);
ddqNull = zeros(nJnt, 1);
s = 1.0;
sStar = 0.0;
exitCode = 1;
tol = 1e-6;
limitExceeded = true;
while limitExceeded == true
    limitExceeded = false;
    ddq = ddqNull + pinv(J*W) * (ddxGoal -dJdq - J*ddqNull);

    if any(ddq < (ddqLow - tol)) || any(ddq > (ddqUpp + tol))
        limitExceeded = true;
    end
    a = pinv(J*W) * ddxGoal;
    b = ddq - a;

    marginL = ddqLow - b;
    marginU = ddqUpp - b;
    sMax = zeros(nJnt, 1);
    for i=1:nJnt
        if W(i,i) == 0
            sMax(i) = inf;
        else
            sMax(i) = FindScaleFactor(marginL(i), marginU(i), a(i));
            % if scale factor is 1 but joint limit is violated,
            % set a scale factor to a small number so that the
            % corresponding joint will be saturated
            if(sMax(i) == 1 && (ddq(i) < (ddqLow(i) - tol) || ddq(i) > (ddqUpp(i) + tol)))
                sMax(i) = 1e-3;
            end
        end
    end

    [~, jntIdx] = min(sMax);
    taskScale = sMax(jntIdx);
    % if the current best scaled solution violates the limit,
    % update sStar to the latest scale factor
    ddqTmp = ddqNull + pinv(J*W)*(sStar * ddxGoal - dJdq - J * ddqNull);
    if ((taskScale > sStar) || (any(ddqTmp > ddqUpp + tol) || any(ddqTmp < ddqLow - tol)))
        sStar = taskScale;
        Wstar = W;
        ddqNullStar = ddqNull;
    end

    W(jntIdx, jntIdx) = 0;
    ddqNull(jntIdx) = min(max(ddqLow(jntIdx), ddq(jntIdx)), ddqUpp(jntIdx));

    if rank(J*W) < nTask
        s = sStar;
        W = Wstar;
        ddqNull = ddqNullStar;
        ddq = ddqNull + pinv(J*W)*(s * ddxGoal - dJdq - J * ddqNull);
        limitExceeded = false;
    end
end
end
