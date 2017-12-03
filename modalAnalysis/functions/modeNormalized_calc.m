%% Description: modeNormalized_calc (sub function)
% ------------------------------------------------------------------------------------------------------------
% This function calculates the normalized modeshape vectors.
% ------------------------------------------------------------------------------------------------------------
% Input:
%   complexMode             = complex modeshape vectors
%   n_physicalStablePoles   = number of physical stable poles
%   n_outputs               = number of investigated output channels
% Output:
%   complexMode_normalized  = normalized modeshape vetors (complex valued)
%   mode_normalized         = normalized modeshape vetors (real valued)
% ------------------------------------------------------------------------------------------------------------
% Copyright 2017 Raphael Frederik Trumpp
% ------------------------------------------------------------------------------------------------------------
% This file is part of modalAnalysis_pLSCF_main.
%
% modalAnalysis_pLSCF_main is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% modalAnalysis_pLSCF_main is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with modalAnalysis_pLSCF_main.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------------------------------------
function [complexModes_normalized, modes_normalized] = modeNormalized_calc(complexModes)

% Returns maximum value (= maximum magnitude) of each complexMode estimation.
[amplitudeModes_max] = max(complexModes);

% Normalization of the complex modeshape vetors by dividing each modeshape by its maximum entry
complexModes_normalized = complexModes ./ amplitudeModes_max;

% Normalized modes are the real part of complexMode_normalized
modes_normalized = abs(real(complexModes_normalized));

 % Phase
phase_complex = angle(complexModes_normalized);          

for i_physicalStablePoles = 1 : size(complexModes, 2)
    
    isInPhase = (phase_complex(:,i_physicalStablePoles) <= pi/2) & (phase_complex(:,i_physicalStablePoles) >= -pi/2);

    modes_normalized(isInPhase == 0, i_physicalStablePoles) = - modes_normalized(isInPhase == 0, i_physicalStablePoles);

end

end