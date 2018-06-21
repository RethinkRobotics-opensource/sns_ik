function [dq, s, exitCode] = snsIk_vel_opt(dqLow, dqUpp, dxGoal, J)
% [dq, s, exitCode] = snsIk_vel_opt(dqLow, dqUpp, dxGoal, J)
%
% This function implements the optimal version of the SNS-IK velocity
% solver for single task.
%
% INPUTS:
%   dqLow (nJnt x 1): lower limit for joint velocities
%   dqUpp (nJnt x 1): upper limit for joint velocities
%   dxGoal (ndx x 1): task-space velocity
%   J (ndx x nJnt): jacoabian mapping joint space to task space
%
% OUTPUTS:
%   dq (nJnt x 1): joint velocity solution with maximum task scale factor
%   s = task scale factor [0, 1]
%   exitCode    (1 == success)
%
%
% NOTES:
%   s * dxGoal = J * dq;

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
dqNull = zeros(nJnt, 1);
s = 1.0;
sStar = 0.0;
exitCode = 1;
tol = 1e-6;
limitExceeded = true;
while limitExceeded == true
    limitExceeded = false;
    P = eye(nJnt) - pinv(J*W)*J;
    dq = pinv(J*W) * dxGoal + P*dqNull;
    if any(dq < (dqLow - tol)) || any(dq > (dqUpp + tol))
        limitExceeded = true;
    end
    a = pinv(J*W) * dxGoal;
    b = dq - a;

    marginL = dqLow - b;
    marginU = dqUpp - b;
    sMax = zeros(nJnt, 1);
    for i=1:nJnt
        if W(i,i) == 0
            sMax(i) = inf;
        else
            sMax(i) = FindScaleFactor(marginL(i), marginU(i), a(i));
        end
    end

    [~, jntIdx] = min(sMax);
    taskScale = sMax(jntIdx);
    if taskScale > sStar
        sStar = taskScale;
        Wstar = W;
        dqNullStar = dqNull;
    end

    W(jntIdx, jntIdx) = 0;
    dqNull(jntIdx) = min(max(dqLow(jntIdx), dq(jntIdx)), dqUpp(jntIdx));

    if rank(J*W) < nTask
        s = sStar;
        W = Wstar;
        dqNull = dqNullStar;
        dq = dqNull + pinv(J*W)*(s * dxGoal - J * dqNull);
        limitExceeded = false;
    end

    % this check is added for optimal version of algorithm
    % do this only if at least two joints are saturated
    if(nJnt > 2)
        if(sum(diag(W)) < nJnt-2)
            mu = -P'*dq;
            for i=1:nJnt
                if W(i,i) == 0
                    if((abs(dq(i)-dqLow(i)) < 1e-6 && mu(i) > 0) || ...
                            (abs(dq(i)-dqUpp(i)) < 1e-6 && mu(i) < 0))
                        % remove i-th joint from saturated joint set
                        W(i,i) = 1;
                        dqNull(i) = 0;
                        limitExceeded = true;
                    end
                end
            end
        end
    end
end
end
