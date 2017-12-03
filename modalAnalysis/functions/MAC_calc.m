
%% Description: MAC_calc (sub function)
% ------------------------------------------------------------------------------------------------------------
% Input:
%   PSDoFRF                     = PSD or FRF matrix
%   frequencyBand               = discrete frequencies [Hz].
%   naturalFrequencies_pLSCF    = natural frequencies found with the pLSCF
%                                 method
%  
% Output:
%   complexMode_normalized_BFD  = power spectral density matrix
% ------------------------------------------------------------------------------------------------------------
% Programmer: Trumpp, Raphael F.
% Email: raphael.trumpp@web.de
% Last modification: 2017/07/17
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
function [MACvalue] = MAC_calc(complexMode_normalized_pLSCF,complexMode_normalized_BFD)

% Shorter names for the mode shape estiamtes
modeShapeA = complexMode_normalized_pLSCF;
modeShapeB = complexMode_normalized_BFD;

% Number of different poles / mode shape
n_modeShapeA = size(complexMode_normalized_pLSCF, 2);
n_modeShapeB = size(complexMode_normalized_BFD, 2);

% Preallocation
MACvalue = zeros(n_modeShapeA,n_modeShapeB);

for i_modeShapeA = 1 : n_modeShapeA
    for i_modeShapeB = 1 : n_modeShapeB
        
        MACvalue(i_modeShapeA, i_modeShapeB)=...
            (abs(modeShapeA(:, i_modeShapeA)' * modeShapeB(:, i_modeShapeB))) ^ 2 /...
            ((modeShapeA(:, i_modeShapeA)' * modeShapeA(:, i_modeShapeA))...
            * (modeShapeB(:, i_modeShapeB)' * modeShapeB(:, i_modeShapeB)));
        
    end   
end

end