function taskScale = FindScaleFactor(low, upp, a)
% taskScale = FindScaleFactor(low, upp, a)
%
% This function computes task scale factor from upper and lower margins and
% the desired task for a single component. This function is called by all
% of SNS IK algorithms. 
%
% INPUTS:
%   low: lower margin
%   upp: upper margin
%   a: desired task
%
% OUTPUTS:
%   taskScale: task scale factor [0, 1] 
%
% NOTES:
%
% This function is a modification of the SNS-IK algorithm 2 that is
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

if (a < 1e10 && a > -1e10)
    
    if a < 0 && low < 0
        
        if a < low
            taskScale = low / a; % task is feasible with scaling
        else
            taskScale = 1.0; % task is feasible without scaling
        end
        
    elseif a > 0 && upp > 0
        
        if upp < a
            taskScale = upp / a; % task is feasible with scaling
        else
            taskScale = 1.0; % task is feasible without scaling
        end
        
    else
        taskScale = 0.0; % task is infeasible
    end
    
else
    taskScale = 0.0; % task is infeasible
end

end