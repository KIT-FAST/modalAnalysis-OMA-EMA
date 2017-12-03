
%% Description: pLSCF_findPoles_calc (sub function)
% ------------------------------------------------------------------------------------------------------------
% This function determines the poles and the natural
% frequencies of the system based on the pLSCF method. 
% ------------------------------------------------------------------------------------------------------------
% Input:
%   PSD                     = power spectral density matrix
%   delta_t                 = sampeling time [s]
%   n_modelOrder            = selected model order of the system
%   n_outputs               = number of investigated output channels
%   modelOrder_max          = maximum model order
%   matrix_gamma            = matrix which contains the coefficients of
%                             gamma
%   matrix_yipsilon         = matrix which contains the coefficients of
%                             yipsilon
% Output:
%   poles_tDomain           = identified poles (in frequency domain)  
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
function [poles_fDomain] =...
    pLSCF_findPoles_calc(n_modelOrder, modelOrder_max, n_outputs, n_outputsRed, deltaTime, matrixGamma_max, matrixYpsilon_max)

%% 1. Calculation of matrices [matrixR_o], [matrixS_o] and [matrixT_o]

matrixGamma = matrixGamma_max(:, 1 : n_modelOrder + 1);
    
matrixYpsilon = matrixYpsilon_max(:, 1 : (n_modelOrder + 1) * n_outputsRed, :);

% Computation of matrixR_o
matrixR_o(:,:) = matrixGamma' * matrixGamma;

% Preallocation
matrixS_o = zeros(n_modelOrder + 1, (n_modelOrder + 1) * n_outputsRed, n_outputs);
matrixT_o = zeros((n_modelOrder + 1) * n_outputsRed, (n_modelOrder + 1) * n_outputsRed, n_outputs);

% Computation of matrixS_o and matrixT_o
for i_outputs = 1 : n_outputs
    
    matrixS_o(:, :, i_outputs) = matrixGamma' * matrixYpsilon(:, :, i_outputs);
    
    matrixT_o(:, :, i_outputs) = matrixYpsilon(:, :, i_outputs)' * matrixYpsilon(:, :, i_outputs);

end

% The searched coefficients of the common denominator model are assumpted
% to be real valued. Through that, only the real parts of the matrices
% matrix_R_o, matrix_R_o and matrix_R_ has to be considered.
matrixR_o = real(matrixR_o);
matrixT_o = real(matrixT_o);
matrixS_o = real(matrixS_o);

%% 2. Calculation and solving the reduced normal equation [matrixAlpha]

% Preallocation
matrixM = zeros((n_modelOrder + 1) * n_outputsRed, (n_modelOrder + 1) * n_outputsRed, n_outputs);

% The matrix_M describes the reduced normal equation as follows:
% [matrix_M][matrix_alpha] = [0]. This is the system of equations which
% has to be solved. The result for the matrix [matrix_alpha] contains the
% searched coefficents of the common denominator model.
for i_outputs = 1 : n_outputs
    
    matrixM(:, :, i_outputs) =...
        matrixT_o(:, :, i_outputs) - (matrixS_o(:, :, i_outputs)') / (matrixR_o(:, :)) * matrixS_o(:, :, i_outputs);

end

% Sum over all output sensors
matrixM = sum(matrixM, 3);

matrixAlpha =...
    - (matrixM(n_outputsRed + 1 : (n_modelOrder + 1) * n_outputsRed, n_outputsRed + 1 : (n_modelOrder + 1) * n_outputsRed))...
    \ (matrixM(n_outputsRed + 1 : (n_modelOrder + 1) * n_outputsRed, 1 : n_outputsRed));

% The lowest order coefficent is constrained to be equal to the identity
% matrix of dimension n_outputs x n_outputs. Through that, the parameter
% redundancy is removed. 
matrixAlpha = cat(1,eye(n_outputsRed), matrixAlpha);

%% 3. Calculation of the denominator coefficents

% Preallocation
matrixA_o = zeros(n_outputsRed, n_outputsRed);

% Extracting of the searched coefficients of the common denominator model
% from the matrix [matrix_A_o]. The coefficients are described through
% matrices of dimension n_outputs x n_outputs.

for i_modelOrder = 0 : n_modelOrder
    
    matrixA_o(:, :, i_modelOrder + 1)= ...
        matrixAlpha(i_modelOrder * n_outputsRed + 1 : (i_modelOrder + 1) * n_outputsRed, 1 : n_outputsRed);    
   
end    

%% 4.Calculation of the companion matrix

% The roots of the denominator polynomial are the eigenvalues of the
% corresponding companion matrix. The companion matrix is a square
% n_outputs * n_modelOrder x n_outputs * n_modelOrder matrix and models a
% dynamic system with n_outputs * n_modelOrder / 2 modes.

%Preallocation
companionMatrix_coeff = zeros(n_outputsRed, n_outputsRed, n_modelOrder);

for i_modelOrder = 1 : n_modelOrder  
    
    companionMatrix_coeff(:, :, i_modelOrder) =...
        - (matrixA_o(:, :, n_modelOrder + 1)) \ matrixA_o(:, :, n_modelOrder + 1 - i_modelOrder);
    
end

%Preallocation
companionMatrix = zeros(n_outputsRed * n_modelOrder, n_outputsRed * n_modelOrder);

% Construction of the companion matrix of the denominator polynomial
for i_modelOrder = 1 : n_modelOrder 
    
    companionMatrix(1 : n_outputsRed, (i_modelOrder - 1) * n_outputsRed + 1: i_modelOrder * n_outputsRed) =...
        companionMatrix_coeff(:, :, i_modelOrder); 
    
end

% Filling the matrix with identity matrices
companionMatrix(n_outputsRed + 1 : n_modelOrder * n_outputsRed, 1 : n_outputsRed * n_modelOrder - n_outputsRed) =...
    eye(n_modelOrder * n_outputsRed - n_outputsRed);

%% 5. Estimation of the eigenfrequencies
% n_modelOrder-th order matrix fraction model yields n_outputs *
% n_modelOrder poles. The poles are obtained trough the eigenvalue
% decomposition of the companion matrix. They are descriebed in z-domain.

% Poles in z-domain
poles_zDomain = eig(companionMatrix);

% Filling with zeros so that returned vectors have always the same length
poles_fDomain = zeros(modelOrder_max * n_outputsRed, 1);

% Transformation of the poles into frequency-domain
poles_fDomain(1 : size(poles_zDomain)) = log(poles_zDomain) / deltaTime;