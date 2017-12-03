%% Description: EMA_pLSCF_modeShapes_calc (sub function)
% ------------------------------------------------------------------------------------------------------------
% This function calculates the complex mode shapes, the residues of the
% pole-residue model and the matrix capitalLambda. 
% ------------------------------------------------------------------------------------------------------------
% Input:
%   physicalStablePoles_sort= list which contains all physical stable poles
%   PSDoFRF                 = FRF matrix
%   frequencyBand           = discrete frequencies [Hz]
%
% Output:
%   complexMode             = matrix which contains the complex modeshapes
%                             in each column
%   residue                 = residue matrix
%   capitalLambda           = special form of the PSD matrix
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
function [complexModes, residue, capitalLambda] = ...
    EMA_pLSCF_modeShapes_calc(physicalStablePoles_pLSCF,PSDoFRF, frequencyBand)

% Number of major stable poles
n_physicalStablePoles = size(physicalStablePoles_pLSCF, 1);

n_frequencyBand = size(frequencyBand,1);

% The vector physicalStablePoles_all_con contains in the first
% n_physicalStablePoles entries the complex poles and in second
% n_physicalStablePoles entries the complex conjugate ot the complex poles.
physicalStablePoles_pLSCF_all = cat(1, physicalStablePoles_pLSCF, conj(physicalStablePoles_pLSCF));

% Preallocation
capitalLambda = zeros(n_frequencyBand, 2 * n_physicalStablePoles);

% Computation of gLambda
for i_physicalStablePoles = 1 : 2 * n_physicalStablePoles
    capitalLambda(:,i_physicalStablePoles) =...
        1 ./ (1i * 2 * pi * frequencyBand(:) - physicalStablePoles_pLSCF_all(i_physicalStablePoles));
end

% Calculation of the lower and upper residue
capitalLambda(:, i_physicalStablePoles + 1) = 1;
capitalLambda(:, i_physicalStablePoles + 2) = -1 ./ (1i*2*pi*frequencyBand(:)) .^ 2;

% Avoiding dividing with zero
capitalLambda(1, i_physicalStablePoles + 2) = capitalLambda(2, i_physicalStablePoles + 2);

% Rearrangement of the FRF matrix to an 2-D matrix. First dim is frequency,
% second dim sensors.
FRF_res = squeeze(PSDoFRF);

% Computation of the residues
residue = capitalLambda \ FRF_res;

% Computation of the complex modeshapes
complexModes = residue(1 : n_physicalStablePoles, :).';

end