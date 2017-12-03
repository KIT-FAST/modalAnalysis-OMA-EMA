%% Description: OMA_PSDmultiSensorLayout_calc (sub function)
% ------------------------------------------------------------------------------------------------------------
% This function calculates the PSD matrix of the measured dataset.
% ------------------------------------------------------------------------------------------------------------
% Input:
%   records                 = measured dataset
%   n_outputsRef            = number of reference output channels
%   window_opt              = specification of the window type for FFT, []
%                             means a Hamming window 
%   noverlap-op             = number of overlapped samples
%   nfft_opt                = number of points to which the FFT is
%                             calculated, determines frequency resolution
%   fs                      = sampeling frequency [Hz]
% --> for further information about window, noverlap and nfft see 'doc cpsd'.
%
% Output:
%   PSD                     = power spectral density matrix of the multi
%                             sensor layout.
%   frequencyBand           = identified natural frequencies [Hz].
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
function [PSD_multiSensor, frequencyBand] =...
    OMA_PSDmultiSensor_calc(records, n_outputsRef, window_opt, noverlap_opt, nfft_opt, fs)
%% Basic values

% Number of patches
n_patches = size(records, 3);

% Number of all sensors (reference and moved) per patch
n_outputsPatch = size(records, 2);

% Number of moved sensors per patch
n_outputsMov = n_outputsPatch - n_outputsRef;

%% Calculation of the PSD matrices of each patch

% Calculation of the frequency band of the PSD matrix
[~, frequencyBand] = cpsd(records(:, 1, 1), records(:, 1, 1), window_opt, noverlap_opt, nfft_opt, fs, 'onesided');

% Number of discrete frequency entries
n_frequencyBand = size(frequencyBand, 1);

% Preallocation
PSDref = zeros(n_frequencyBand, n_outputsRef, n_outputsRef, n_patches);

% Calculation of PSDref. The matrix PSDref contains all auto- and
% cross PSD estimates between the reference sensors.
for i_patches =1 : n_patches
    for i_outputsRef = 1 : n_outputsRef
        for j_outputsRef = 1 : n_outputsRef
            
            [PSDref(:, j_outputsRef, i_outputsRef, i_patches), ~] = ...
                cpsd(records(:, i_outputsRef, i_patches), records(:, j_outputsRef, i_patches), window_opt, noverlap_opt, nfft_opt, fs, 'onesided');
            
        end
    end
end

% Preallocation
PSDmov = zeros(n_frequencyBand,n_outputsRef,n_outputsMov,n_patches);

% Calculation of PSDmov. The matrix mov_PSD contains all cross-PSD
% estimates between the moved sensors and the reference sensors. 
for i_patches = 1 : n_patches
    for i_outputsRef = 1 : n_outputsMov
        for j_outputsRef = 1 : n_outputsRef
            
            [PSDmov(:, j_outputsRef, i_outputsRef, i_patches), ~] =...
                cpsd(records(:, n_outputsRef + i_outputsRef, i_patches), records(:, j_outputsRef, i_patches), window_opt, noverlap_opt, nfft_opt, fs, 'onesided');
            
        end
    end
end

%% Re-scaling of the PSD matrices to a single PSD matrix

PSDrescal_zero = zeros(n_outputsRef, n_outputsRef, n_frequencyBand);
PSDrescal_coeff = zeros(n_outputsMov, n_outputsRef, n_frequencyBand, n_patches);
PSDrescal_res = zeros(n_outputsPatch, n_outputsRef, n_frequencyBand);


% Rearrange of the matrices
PSDmov = permute(PSDmov, [3, 2, 1, 4]);
PSDref = permute(PSDref, [3, 2, 1, 4]);

% Re-scaling for every discrete frequency entry
for n_frequencyBand = 1 : size(frequencyBand)
    
    % Calculation of coefficients for the re-scaled data set

    % Re-scaling factor for every patch. Averaged over all ref_PSD from all
    % patches.
    PSDrescal_zero(:, :, n_frequencyBand) = sum(PSDref(:, :, n_frequencyBand, :), 4) / n_patches;
    
    % Calculation of the re-scaled coefficients
    for i_patches = 1 : n_patches     
        
        PSDrescal_coeff(:, :, n_frequencyBand, i_patches) = ...
            PSDmov(:, :, n_frequencyBand, i_patches) / PSDref(:, :, n_frequencyBand, i_patches) * PSDrescal_zero(:, :, n_frequencyBand);
        
    end
    
    % First entry of the re-scaled PSD
    PSDrescal_res(1 : n_outputsRef, :, n_frequencyBand) = PSDrescal_zero(:, :, n_frequencyBand);

    % Construction of the re-scaled PSD with the previous calculated
    % coefficients.
    for i_patches = 1 : n_patches  
        
    PSDrescal_res((i_patches - 1) * n_outputsMov + n_outputsRef + 1 : i_patches * n_outputsMov + n_outputsRef, :, n_frequencyBand)...
        = PSDrescal_coeff(:, :, n_frequencyBand, i_patches);
    
    end

end

%% Rearrangement of the PSD matrix
% Special structure of the PSD matrix is needed for pLSCF_main.m function.
% PSD_multiSensor is the re-scaled matrix, which is needed for other OMA methods

PSD_multiSensor = permute(PSDrescal_res, [3, 2, 1]);

end