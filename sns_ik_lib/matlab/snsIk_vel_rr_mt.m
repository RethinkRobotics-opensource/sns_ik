function [dq, sData, exitCode] = snsIk_vel_rr_mt(dqLow, dqUpp, dxGoalData, JData)
% [dq, sData, exitCode] = snsIk_vel_rr_mt(dqLow, dqUpp, dxGoalData, JData)
%
% This function implements a generic multi-task version of the SNS-IK
% velocity solver that tries to meet multiple objective terms
% in a prioritized manner.
%
% INPUTS:
%   dqLow (nJnt x 1):  lower limit for joint velocities
%   dqUpp (nJnt x 1):  upper limit for joint velocities
%   dxGoalData (cell array): task-space velocities of all the tasks in the order of priorities
%   JData (cell array): Jacobians of all the task in the order of priorities
%
% OUTPUTS:
%   dq (nJnt x 1): joint velocity solution with maximum task scale factor
%   sData (nTask x 1): task scale factors [0, 1] for all the tasks
%   exitCode    (1 == success)
%
% NOTES:
%
% This function is a modification of the multi-task SNS-IK algorithm that is
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

% get the number of tasks
nTask = size(JData,1);
WData = cell(nTask,1);
sData = zeros(nTask,1);
dqData = cell(nTask,1);
dqNullData = cell(nTask,1);

% initialization
exitCode = 1;
tol = 1e-6;

nJnt = numel(dqLow);
I = eye(nJnt);
Pi = I;
dqi = zeros(nJnt,1);

for iTask = 1:nTask

    % get i-th task jacobian
    Ji = JData{iTask};
    ndxGoal = size(Ji,1);

    % get i-th task velocity
    dxGoali = dxGoalData{iTask};

    % update variables for previous projection matrix and solution
    PiPrev = Pi;
    dqiPrev = dqi;

    % initialize variables for i-th task
    Wi = eye(nJnt);
    dqNulli = zeros(nJnt, 1);
    PiBar = PiPrev;
    si = 1.0;
    siStar = 0.0;

    limitExceeded = true;
    drop_correction_term = false;
    cntLoop = 1;
    while limitExceeded == true
        limitExceeded = false;

        % compute a solution without task scale factor
        PiHat = (I - pinv(Ji*PiBar, tol)*Ji)*pinv((I - Wi)*PiPrev, tol);

        if (drop_correction_term)
            dqi = dqiPrev + pinv(Ji*PiBar, tol) * dxGoali + ...
                PiHat*(dqNulli - dqiPrev);
        else
            dqi = dqiPrev + pinv(Ji*PiBar, tol) * (dxGoali - Ji*dqiPrev) + ...
                PiHat*(dqNulli - dqiPrev);
        end

        % check whether the solution violates the limits
        if any(dqi < (dqLow - tol)) || any(dqi > (dqUpp + tol))
            limitExceeded = true;
        end

        % compute scale factor
        a = pinv(Ji*PiBar, tol) * dxGoali;
        b = dqi - a;

        marginL = dqLow - b;
        marginU = dqUpp - b;
        sMax = zeros(nJnt, 1);
        for iJnt=1:nJnt
            if Wi(iJnt,iJnt) == 0
                sMax(iJnt) = inf;
            else
                sMax(iJnt) = FindScaleFactor(marginL(iJnt), marginU(iJnt), a(iJnt));
            end
        end

        [~, jntIdx] = min(sMax);
        taskScale = sMax(jntIdx);

        if (isinf(taskScale))
            % this means that all joints are saturated, so
            % lower priority tasks become infeasible
            taskScale = 0;
        end

        % do the following only if the task is feasible and the scale
        % factor caculated is correct
        if (iTask == 1 || taskScale > 0)

            % try to store maximum scale factor with associated variables
            if taskScale > siStar
                siStar = taskScale;
                Wistar = Wi;
                dqNulliStar = dqNulli;
                PiBarStar = PiBar;
                PiHatStar = PiHat;
            end

            % saturate the most critical joint
            Wi(jntIdx, jntIdx) = 0;
            dqNulli(jntIdx) = min(max(dqLow(jntIdx), dqi(jntIdx)), dqUpp(jntIdx));

            % update projection matrices
            PiBar = (I - pinv((I - Wi)*PiPrev, tol))*PiPrev;
            PiHat = (I - pinv(Ji*PiBar, tol)*Ji)*pinv((I - Wi)*PiPrev, tol);

            % if rank is below task dimension, terminate the loop and output
            % the current best solution
            if (rank(Ji*PiBar) < ndxGoal)
                si = siStar;
                Wi = Wistar;
                dqNulli = dqNulliStar;
                PiBar = PiBarStar;
                PiHat = PiHatStar;

                if (drop_correction_term)
                    dqi = dqiPrev + pinv(Ji*PiBar, tol) * si*dxGoali + ...
                        PiHat*(dqNulli - dqiPrev);

                    limitExceeded = false;
                else
                    dqi = dqiPrev + pinv(Ji*PiBar, tol) * (si*dxGoali - Ji*dqiPrev) + ...
                        PiHat*(dqNulli - dqiPrev);

                    % if the solution violates the limits, drop the
                    % correction term and run it again
                    if (any(dqi > dqUpp + tol) || any(dqi < dqLow - tol))
                        drop_correction_term = true;
                        cntLoop = 0;
                        Wi = eye(nJnt);
                        dqNulli = zeros(nJnt, 1);
                        PiBar = PiPrev;
                        si = 1.0;
                        siStar = 0.0;
                    else
                        limitExceeded = false;
                    end
                end

            end

        else % if the current task is infeasible
            si = 0;
            Wi = zeros(nJnt, nJnt);
            dqi = dqiPrev;
            limitExceeded = false;
        end

        cntLoop = cntLoop + 1;
    end

    if (si > 0)
        % update nullspace projection
        Pi = PiPrev - pinv(Ji*PiPrev, tol)*(Ji*PiPrev);
    end

    % store data
    WData{iTask} = Wi;
    sData(iTask) = si;
    dqData{iTask} = dqi;
    dqNullData{iTask} = dqNulli;
end

%-- end of algorithm
dq = dqi;

end
