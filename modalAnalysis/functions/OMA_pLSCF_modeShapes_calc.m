%% Description: OMA_pLSCF_modeShapes_calc (sub function)
% ------------------------------------------------------------------------------------------------------------
% This function calculates the complex mode shapes, the residues of the
% pole-residue model and the matrix capitalLambda.
% ------------------------------------------------------------------------------------------------------------
% Input:
%   physicalStablePoles_sort= list which contains all physical stable poles
%   PSDoFRF                 = PSD matrix
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
    OMA_pLSCF_modeShapes_calc(physicalStablePoles_sort,PSDoFRF, frequencyBand)

% Number of output channels
n_outputs = size(PSDoFRF, 3);

% For multi sensor layout: Reduced dimension (only reference sensors)
n_outputsRed = size(PSDoFRF, 2);

% Number of major stable poles
n_physicalStablePoles = size(physicalStablePoles_sort, 1);

% Number of discrete frequencies
n_frequencyBand = size(frequencyBand,1);

%% 1. Calculation of capitalLambda

% Preallocation
capitalLambda_coeff1 = zeros(n_frequencyBand * n_outputs, n_outputs, n_physicalStablePoles);
capitalLambda_coeff2 = zeros(n_frequencyBand * n_outputs, n_outputs, n_physicalStablePoles);

% Computation of capitalLambda
for i_physicalStablePoles = 1 : n_physicalStablePoles
    for i_frequencyBand = 1 : n_frequencyBand
        
        capitalLambda_coeff1((i_frequencyBand - 1) * n_outputs + 1 : i_frequencyBand * n_outputs, :, i_physicalStablePoles) = ...
            eye(n_outputs) ./ (1i * 2 * pi * frequencyBand(i_frequencyBand) - physicalStablePoles_sort(i_physicalStablePoles));
        
        capitalLambda_coeff2((i_frequencyBand - 1) * n_outputs + 1 : i_frequencyBand * n_outputs, :, i_physicalStablePoles) = ...
            eye(n_outputs) ./ (1i * 2 * pi * frequencyBand(i_frequencyBand) - conj(physicalStablePoles_sort(i_physicalStablePoles)));                 
        
    end    
end

% The matrices capitalLambda_coeff3 and capitalLambda_coeff4 are complex conjugates of
% the previously calculated coefficients.
capitalLambda_coeff3 = conj(capitalLambda_coeff2);
capitalLambda_coeff4 = conj(capitalLambda_coeff1);

% Preallocation
capitalLambda = zeros(n_outputs * n_frequencyBand, 4 * n_outputs * n_physicalStablePoles);

for i_physicalStablePoles = 1 : n_physicalStablePoles
    
    capitalLambda(:, (i_physicalStablePoles - 1) * n_outputs + 1 : i_physicalStablePoles * n_outputs) =...
        capitalLambda_coeff1(:, :, i_physicalStablePoles);
    
    capitalLambda(:, (i_physicalStablePoles - 1) * n_outputs + 1 + n_outputs * n_physicalStablePoles : i_physicalStablePoles * n_outputs + n_outputs * n_physicalStablePoles) =...
        capitalLambda_coeff2(:, :, i_physicalStablePoles);
   
    capitalLambda(:, (i_physicalStablePoles - 1) * n_outputs + 1 + 2 * n_outputs * n_physicalStablePoles : i_physicalStablePoles * n_outputs + 2 * n_outputs * n_physicalStablePoles) =...
        capitalLambda_coeff3(:,:,i_physicalStablePoles);
    
    capitalLambda(:, (i_physicalStablePoles - 1) * n_outputs + 1 + 3 * n_outputs * n_physicalStablePoles : i_physicalStablePoles*n_outputs + 3 * n_outputs * n_physicalStablePoles) =...
        capitalLambda_coeff4(:, :, i_physicalStablePoles);
end

%% 2. Calculation of gLambda

% Preallocation 
gLambda = zeros(i_frequencyBand * n_outputs, n_outputsRed);

% Calculation of gLambda
for i_frequencyBand = 1 : n_frequencyBand
    
    gLambda((i_frequencyBand - 1) * n_outputs + 1 : i_frequencyBand * n_outputs, :) =...
        permute(PSDoFRF(i_frequencyBand, :, :), [3, 2, 1]); 
    
end

%% 3. Estimation of the mode shapes

% The calculation of the residue can suffer from a rank deficient matrix.
% But the warning is off, because this doesn't affect the results much.
% warning('off','MATLAB:rankDeficientMatrix')

% Calculation of the residues
residue = capitalLambda \ gLambda;

% Preallocation
complexModes = zeros(n_outputs, n_physicalStablePoles);

% Following loop calculates the SVD of all the coefficients the residue
% matrix contains . Every first column of the corresponding U Vectors is an
% estimate of the corresponding modeshape vector.

for i_physicalStablePoles = 1 : n_physicalStablePoles
    
    [uVector, ~, ~] = svd(residue((i_physicalStablePoles - 1) * n_outputs + 1 : n_outputs * i_physicalStablePoles, :));

    complexModes(:, i_physicalStablePoles) = uVector(:, 1);
    
end

end