%% Description: OMA_PSD_calc (sub function)
% ------------------------------------------------------------------------------------------------------------
% This function calculates the PSD matrix of the measured dataset.
% ------------------------------------------------------------------------------------------------------------
% Input:
%   records                 = measured dataset
%   n_outputs               = number of investigated output channels
%   window_opt              = specification of the window type for FFT, []
%                             means a Hamming window 
%   noverlap_opt            = number of overlapped samples
%   nfft_opt                = number of points to which the FFT is
%                             calculated, determines frequency resolution
%   fs                      = sampeling frequency [Hz]
% --> for further information about window, noverlap and nfft see 'doc cpsd'.
%
% Output:
%   PSD                     = power spectral density matrix
%   frequencyBand           = discrete frequencies [Hz].
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
function [PSD, frequencyBand] = OMA_PSD_calc(records, window_opt, noverlap_opt, nfft_opt, fs)

% Number of output channels
n_outputs = size(records, 2);

[~, frequencyBand] = ...
        cpsd(records(:, 1), records(:, 1), window_opt, noverlap_opt, nfft_opt, fs, 'onesided');
    
n_frequencyBand = size(frequencyBand,1);    
    
% Preallocation
PSD = zeros(n_frequencyBand, n_outputs, n_outputs);

% Calculation of the PSD matrix
for i_outputs = 1:n_outputs
    for j_outputs = 1:n_outputs
        
    [PSD(:, j_outputs, i_outputs), ~] = ...
        cpsd(records(:, i_outputs), records(:, j_outputs), window_opt, noverlap_opt, nfft_opt, fs, 'onesided');
    
    end
end

end