%% Description: BFD_modeShapes_calc (sub function)
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
function [complexMode_normalized_BFD] = BFD_modeShapes_calc(PSDoFRF, frequencyBand, naturalFrequencies)

% Basic values:

% Number of output channels
n_outputs = size(PSDoFRF, 3);

% Reference sensor for mode shape estimation
refPoint = 1;

% Numver of natural frequencies
n_naturalFrequencies_pLSCF = size(naturalFrequencies, 1);

%% 1. Finding the indices of the natural frequencies

%Preallocation
autoPSDoFRF_peaksIndex = zeros(n_naturalFrequencies_pLSCF, 1);

% Identification of the peaks in the set intervals.
for i_naturalFrequencies_pLSCF = 1 : n_naturalFrequencies_pLSCF
    
    % Identifying the frequencies of the frequency band which are
    % nearest to the set intervall.
    [~, autoPSDoFRF_peaksIndex(i_naturalFrequencies_pLSCF)] = ...
        min(abs(frequencyBand-naturalFrequencies(i_naturalFrequencies_pLSCF)));
    
end

%% 2. Estimation of the mode shapes

% Preallocation
complexMode = zeros(n_outputs, n_naturalFrequencies_pLSCF);

% Identifying the mode comples modeshapes corresponding to the identified
% poles respectively natural frequencies
for i_naturalFrequencies_pLSCF = 1 : n_naturalFrequencies_pLSCF
    
    complexMode(:, i_naturalFrequencies_pLSCF) = PSDoFRF(autoPSDoFRF_peaksIndex(i_naturalFrequencies_pLSCF), refPoint, :); 
    
end

%As FDD Method:
for i_naturalFrequencies_pLSCF = 1 : n_naturalFrequencies_pLSCF
    
    PDS_fequLine = zeros(n_outputs);
    
    for i_outputs = 1 : n_outputs
        PDSoFRF_fequLine(i_outputs,:) = PSDoFRF(autoPSDoFRF_peaksIndex(i_naturalFrequencies_pLSCF),:,i_outputs); %#ok<SAGROW>
    end
    
    [uVector, ~, ~] = svd(PDSoFRF_fequLine);
     
    complexMode(:,i_naturalFrequencies_pLSCF) = uVector(:,1);
end

% Calling the function modeNormalized_calc to obtain the normalized mode
% shapes complexMode_normalized and the real valued modeshape
% mode_normalized.
[complexMode_normalized_BFD, ~] = modeNormalized_calc(complexMode);

end