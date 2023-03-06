%% Description: modalAnalysis_pLSCF_main (main function)
% ------------------------------------------------------------------------------------------------------------
% This code provides you the p-LSCF (poly-reference least square complex
% frequency) method for Operational Modal Analysis (OMA) and Experimental
% Modal Analysis (EMA). This programm will compute the natural frequencies,
% damping ratios and mode shapes of the studied structure.  
% For detailed information about this method see: C. Rainieri, G. Fabbrocino, Operational Modal
% Analysis of Civil Engineering Structures, Springer Science+Business, DOI
% 10.1007/978-1-4939-0767-0.
% ------------------------------------------------------------------------------------------------------------
% This programm was developed as a part of a Bachlor Thesis at the
% Karlsruher Institut of Technology (Germany).
%
% Programmer: Trumpp, Raphael F.
% Email: raphael.trumpp@web.de
% Last modification: 2017/07/17
% ------------------------------------------------------------------------------------------------------------
% Copyright 2017 Raphael Frederik Trumpp
% ------------------------------------------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------------------------------------
%% 1. Clean-Up
clc, clf, clearvars 

%% 2. User input

% Name of the file containing the measurement data/records
fileNameRecords = 'EMA_records_demo.mat';

% OMA or EMA data?
dataType = 'EMA';

% Comparison function with BFD method true or false?
isCompareWithBFD = true;

% Multi sensor layout function true or false?
isMultiSensorLayout = false;
n_outputsRef = 2;

% Setting of the minimum and maxium model order 
modelOrder_min = 5;
modelOrder_max = 30;

% Specify the examined frequency band. frequencyMax must be smaller than
% fs/2.

frequencyBand_min = 20;
frequencyBand_max = 150;

% Sampling frequency [Hz] of the sampeling records
fs = 600;

% Settings for the FFT of the cpsd (welch method) function
window_opt = [];
noverlap_opt = [];
nfft_opt = [];

%% 3. Initialization

disp('Initialization...')
disp('----------------------------------------------------------------------')
disp('User inputs:')

switch dataType
    
    case 'OMA'
    
        disp('--> Data typ: Operational Modal Analysis (OMA).')

        if isMultiSensorLayout == true
            disp('--> Multi sensor layout function "on".') 
        else
            disp('--> Multi sensor layout function "off".')
        end
    
    case 'EMA'
        
        disp('--> Data typ: Eyperimental Modal Analysis (EMA).')   
        
    otherwise
        
        error(['Unknown data type: ', dataType, ' --> either ''OMA'' or ''EMA'' possible'])
end

if isCompareWithBFD == true
    disp('--> Comparison Function is "on".')
else
    disp('--> Comparison Function is "off".') 
end

disp(['--> Minimum model order is set to: "', num2str(modelOrder_min),...
    '" and maximum model order is set to: "' num2str(modelOrder_max),'".']) 
disp(['--> Smallest examined frequency is set to: "', num2str(frequencyBand_min),...
    ' Hz" and highest examined frequency is set to: "' num2str(frequencyBand_max),' Hz".'])

disp('----------------------------------------------------------------------')

% Adding path to sub folders
addpath('./functions');
addpath('./records');

disp('Importing the sampling records.')

% Importing the sample records. Imported Data must be named 'records'
load(fileNameRecords);

disp('--> The recordings were successfully imported.')
disp('----------------------------------------------------------------------')

%% 4. Estimation of the PSD matrix or loading of the FRF matrix

% Different procedures, depending on user input

switch dataType
    
    case 'OMA'

        disp('Calculating the PSD matrix')

        if isMultiSensorLayout == true

            % Calling the function PSD_calc which calculates the PSD of the records
            % (auto and cross spectra between all output channels).
            [PSD_multiSensor, frequencyBand] = ...
            OMA_PSDmultiSensor_calc(records, n_outputsRef, window_opt, noverlap_opt, nfft_opt, fs);

            % Renaming the PSD matrix to an common name (used by OMA and EMA
            PSDoFRF = PSD_multiSensor;

        else

            % Importing the PSD matrix calculated by the function for multi sensor
            % layout.  
            [PSD, frequencyBand] = OMA_PSD_calc(records, window_opt, noverlap_opt, nfft_opt, fs);

            % Renaming the PSD matrix to an common name (used by OMA and EMA
            PSDoFRF = PSD;

        end  

        disp('--> The PSD matrix was successfully calculated.')
        disp('----------------------------------------------------------------------')

    case 'EMA'
    
        disp('Loading the FRF matrix')

        % Renaming the FRF matrix to an common name (used by OMA and EMA)
        % and permuting the dimension to the required structure
        PSDoFRF = permute(FRF, [1 3 2]);

        % Sampling frequency [Hz] of the sampeling records
        fs = 2 * frequencyBand_max;      
        disp('--> The FRF matrix was successfully loaded.')
        disp('----------------------------------------------------------------------')

end

% Basic values:

% Number of output channels
n_outputs = size(PSDoFRF, 3);

% For multi sensor layout: Reduced dimension (only reference sensors)
n_outputsRed = size(PSDoFRF, 2);

% Sampling intervall [s]
deltaTime = 1 / fs;

% Identifying the frequencies of the frequency band which are
% nearest to the set frequency intervall.
[~, frequencyBand_minIndex] = min(abs(frequencyBand - frequencyBand_min));
[~, frequencyBand_maxIndex] = min(abs(frequencyBand - frequencyBand_max));

% Reduced frequency band and PSDoFRF matrix
frequencyBand = frequencyBand(frequencyBand_minIndex : frequencyBand_maxIndex);
PSDoFRF = PSDoFRF(frequencyBand_minIndex : frequencyBand_maxIndex, :, :);

%% 5. Calculation of the poles (lambdaR) for different model orders

% Index of the highest frequency
n_frequencyBand = length(frequencyBand);

disp('Calculating the poles and the natural frequencies.')

% The chosen polynomial basis function is exp(i*2*Pi*frequency*delta_t)=Omega_f=z_f
% in the z-domain. This complex polynomial basis function ensures good numerical
% conditioning. The matrix matrix_gamma has dimension n_frequencyBand x (n_modelOrder + 1).
matrixGamma_max = exp(1i * 2 * pi * frequencyBand * (0:modelOrder_max) * deltaTime);

% Preallocation
matrixYpsilon_max = zeros(n_frequencyBand - 1, (modelOrder_max + 1) * n_outputsRed, n_outputs);

% The operator "kron" is the Kronecker product.The matrix matrix_yipsilon has
% dimension n_frequencyBand x (n_modelOrder + 1) * n_outputs. 
for i_outputs = 1 : n_outputs
    for i_frequencyBand = 1 : n_frequencyBand
        
        matrixYpsilon_max(i_frequencyBand, :, i_outputs) =...
            -kron(matrixGamma_max(i_frequencyBand, :), PSDoFRF(i_frequencyBand, :, i_outputs));
        
    end
end

% Preallocation
lambdaR = zeros(modelOrder_max * n_outputsRed, modelOrder_max);

% Estimation of lambda_r for different modal orders
for i_modelOrder = modelOrder_min : modelOrder_max
             
    [lambdaR(:, i_modelOrder)]...
        = pLSCF_findPoles_calc(i_modelOrder, modelOrder_max, n_outputs, n_outputsRed, deltaTime,...
        matrixGamma_max, matrixYpsilon_max);

end

disp('--> The poles and the natural frequencies were successfully calculated.')
disp('----------------------------------------------------------------------')

%% 6. Identification of the stable poles

disp('Identifying the physical poles.')

% Following statement finds all stable poles from every modal order n. Criterium for
% stable poles: real(lambda_r) < 0. Result is the vector stablePoles_all
% which contains all stable poles (redundancy of the poles is removed by
% just picking the poles of the complex conjugate pairs with positiv imag part).
for i_modelOrder = modelOrder_max : -1 : modelOrder_min
    
    stablePoles_index = (real(lambdaR(:, i_modelOrder)) < 0) & (imag(lambdaR(:, i_modelOrder)) > 0);
    
    stablePoles_unsorted(1:size(find(stablePoles_index)), i_modelOrder) = ...
        lambdaR(stablePoles_index, i_modelOrder);
    
end

stablePoles_unsorted = fliplr(stablePoles_unsorted);

% Sorting the list in ascending order (compares the abs values)
stablePoles = sort(stablePoles_unsorted, 'ComparisonMethod', 'abs');

%% 7. Construction of the stabilization diagram

% Different procedures, depending on user input
switch dataType
    
    case 'OMA'
        % Preallocation
        PSDoFRF_plot = zeros(n_frequencyBand, 1);

        % Calculating PSD_auto_all, which is the sum of all auto-PSDs
        for i_outputs = 1 : n_outputsRed

            PSDoFRF_plot = PSDoFRF_plot + PSDoFRF(:, i_outputs, i_outputs);

        end

        % Averaging to a single representative sensor
        PSDoFRF_plot = PSDoFRF_plot / n_outputsRed;
    
    case 'EMA'

        % Calculating FRF_auto_all, which is the sum of all FRFs and averaging to a
        % single representative sensor. 
        PSDoFRF_plot = sum(abs(PSDoFRF(:, 1, :)), 3) / n_outputs;
  
end

% Plot of the auto-PSD or FRF
subplot(1, 2, 1)
yyaxis left
plot(frequencyBand, PSDoFRF_plot, 'lineWidth', 1.5);

hold on

for i_modelOrder = 1 : modelOrder_max - modelOrder_min

    % Conversion to frequenzy domain 
    frequenciesPlot = abs(stablePoles(:, i_modelOrder)) ./ (2 * pi);
    
    % No displaying of the zero entries
    frequenciesPlot = frequenciesPlot(frequenciesPlot > 0);
    
    if ~isempty(frequenciesPlot)
        
        % Displaying points for each pole at different model orders.
        yyaxis right
        plot (frequenciesPlot,modelOrder_max + 1 - i_modelOrder,...
            '+', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 8);  
        
    end
end

% Title and axis labels of the stabilization diagramm
% Caption of the left y axis
yyaxis left
title('Stabilization diagram', 'FontSize', 26, 'Color', 'k')
xlabel({'Frequency [Hz]'}, 'FontSize', 20);
ylabel({'Magnitude'}, 'FontSize', 20);

% Caption of the right y axis
yyaxis right
ylabel({'Model order'}, 'FontSize', 20)

% Different procedures, depending on user input
switch dataType
    
    case 'OMA'

        % Legend
        legend({'Auto PSD', 'Stable poles'}, 'FontSize', 14, 'TextColor', 'k')
    
    case 'EMA'
    
        % Legend
        legend({'FRF', 'Stable poles'}, 'FontSize', 14, 'TextColor', 'k')
    
end

% Specific settings 
grid on

hold off

%% 8. Identification of the the "physical" stable poles

% Conversion to frequency domain
stablePoles_abs = abs(stablePoles) ./ (2*pi);

% Preallocation
stablePoles_maxOrder = zeros(2, size(stablePoles, 1));
 
% Assuming that all physical stable poles occures at the maximum model order,
% the identification of the physical stable poles starts from the poles,
% which occures at maximum model order. 
stablePoles_maxOrder(1, :) = stablePoles(:, 1); 

% Number of stable poles
n_stablePoles = size(stablePoles, 1);

% Criterum for "major stable pole":
% A major stable pole occuries at 20% of the examined model orders or more. Also its abs value changes
% less than +-1%.
for i_stablePoles = 1 : n_stablePoles
    
    counterStablePoles =...
        0.99 * stablePoles_abs(i_stablePoles, 1) < stablePoles_abs(:) & ...
        stablePoles_abs(:) < 1.01 * stablePoles_abs(i_stablePoles, 1);

    % Number of times a stable pole is occurring over all model orders.
    stablePoles_maxOrder(2, i_stablePoles) = size(find(counterStablePoles),1); 

end

% Number of times a stable pole must occur until it is declared as a physical stable pole.
n_StablePolesRequired = round((modelOrder_max - modelOrder_min) * 0.2);

% Results:

% List with all physical stable poles
physicalStablePoles_pLSCF = stablePoles_maxOrder(1, stablePoles_maxOrder(2, :) >= n_StablePolesRequired).'; % Checken
physicalStablePoles_pLSCF = sort(physicalStablePoles_pLSCF, 'ComparisonMethod', 'abs');

% Number of major stable poles
n_physicalStablePoles_pLSCF = size(physicalStablePoles_pLSCF, 1);

% List with the corresponding natural frequencies to the physical stable poles
naturalFrequencies_pLSCF = abs(physicalStablePoles_pLSCF(:, 1)) / (2 * pi);

% List with the corresponding damping ratios to the physical stable poles
dampingRatio_pLSCF = - real(physicalStablePoles_pLSCF(:, 1)) ./ abs(physicalStablePoles_pLSCF(:, 1)); 

disp('--> The physical stable poles were succesfully found.')
disp('----------------------------------------------------------------------')

disp('The corresponding natural frequencies [Hz] are: ')
disp(naturalFrequencies_pLSCF.')
disp('----------------------------------------------------------------------')

disp('The corresponding damping ratios [%] are: ')
disp(dampingRatio_pLSCF.')
disp('----------------------------------------------------------------------')

%% 9. Calculation of the mode shapes

disp('Determining the associated mode shapes to the identified physical poles.')

switch dataType
    
    case 'OMA'
    
        disp('The following warning doesn´t affect the results much. Rank deficient can´t be avoided.')

        [complexModes_pLSCF, residue, capitalLambda] = ...
            OMA_pLSCF_modeShapes_calc(physicalStablePoles_pLSCF, PSDoFRF, frequencyBand);

    case 'EMA'
        disp('The following warning doesn´t affect the results much. Rank deficient can´t be avoided.')

        [complexModes_pLSCF, residue, capitalLambda] = ...
            EMA_pLSCF_modeShapes_calc(physicalStablePoles_pLSCF, PSDoFRF, frequencyBand);

end

[complexModes_normalized_pLSCF, modes_normalized_pLSCF] = modeNormalized_calc(complexModes_pLSCF);

disp('--> The mode shapes of the system were succesfully identified.')
disp('----------------------------------------------------------------------')

%disp('The complex mode shapes are:')
%disp(complexModes_normalized_pLSCF)
%disp('----------------------------------------------------------------------')

%% 11. PSD/FRF reconstruction 

switch dataType
    
    case 'OMA'
    
        disp('Reconstruction of the PSD matrix.')

        % Reconstruction of the PSD matrix by the parametric model of the LSFD step
        PSD_rec = capitalLambda * residue;

        PSDoFRF_plot_rec = zeros(n_frequencyBand, 1);

        % The trace of the PSD matrix are the auto PSD
        for i_frequencyBand = 1 : n_frequencyBand    
            for i_outputs_red = 1 : n_outputsRed

                PSDoFRF_plot_rec(i_frequencyBand, 1) = PSDoFRF_plot_rec(i_frequencyBand, 1)...
                    + abs(PSD_rec(n_outputs * (i_frequencyBand - 1) + i_outputs_red, i_outputs_red));


            end    
        end

        % Averaging to a single representative sensor
        PSDoFRF_plot_rec = PSDoFRF_plot_rec / n_outputsRed;

        disp('--> Succesfully reconstructed the PSD matrix.')
        disp('----------------------------------------------------------------------')

    case 'EMA'
    
        disp('Reconstruction of the FRF matrix.')

        % Reconstruction of the FRF matrix by the paramteric model of the LSFD step
        FRF_rec = capitalLambda * residue;

        % Sum of all sensors
        PSDoFRF_plot_rec = sum(abs(FRF_rec), 2) / n_outputs;

        disp('--> Succesfully reconstructed the FRF matrix.')
        disp('----------------------------------------------------------------------')

end

% Plot of the measured PSD and reconstructed PSD to compare them with each
% other. 

subplot(1, 2, 2) 
yyaxis left
plot(frequencyBand, PSDoFRF_plot, 'lineWidth', 1.5);

hold on
yyaxis right
plot(frequencyBand, PSDoFRF_plot_rec, 'lineWidth', 1.5);
ylabel({'Magnitude'}, 'FontSize', 20);

% Title and axis labels of the reconstructed FRF curve
yyaxis left

xlabel({'Frequency [Hz]'}, 'FontSize', 20);
ylabel({'Magnitude'}, 'FontSize', 20);

if isequal(dataType, 'OMA') 
    title('PSD Curve', 'FontSize', 26, 'Color', 'k')
else
    title('FRF Curve', 'FontSize', 26, 'Color', 'k')
end

% Legend
legend({'Measured data', 'Parametric model'}, 'FontSize', 14, 'TextColor', 'k')

% Specific settings 
grid on

hold off

%% 12. Comparison of mode shapes with BFD

if isCompareWithBFD == true
    
    disp('Comparison of the mode shape estimates with the BFD method.')

    % Mode shape estimation of the BFD method
    [complexMode_normalized_BFD] = BFD_modeShapes_calc(PSDoFRF, frequencyBand, naturalFrequencies_pLSCF);

    % MAC (modal assurance criterion) compares and assesses the various
    % estimates of the mode shape vectors. MAC values over 0.9 indicates good accuracy of the estimates. MAC values
    % under 0.1 should be avoided and their corresponding mode shape vectors
    % used with caution. 
    [MACvalue] = MAC_calc(complexModes_normalized_pLSCF,complexMode_normalized_BFD);

    % Columns are the estimates of the pLSCF Method, rows from the BFD method
    disp('Columns are the estimates of the pLSCF Method, rows from the BFD method.')

    disp(MACvalue)

    disp('----------------------------------------------------------------------')

end

%% 13. Saving results

% Current date and time is added to the file name.
% Format:YearMonthDayTHourMinuteSecond (ISO 8601)
currentDate = datestr(datetime('now'), 30);

% Path of the folder "results"
pathFolderResults = [pwd filesep 'results' filesep];

switch dataType
    
    case 'OMA'

        save([pathFolderResults, 'OMA_pLSCF_results_' currentDate], 'physicalStablePoles_pLSCF', ...
            'naturalFrequencies_pLSCF', 'dampingRatio_pLSCF', 'complexModes_normalized_pLSCF','modes_normalized_pLSCF', ...
            'modelOrder_min', 'modelOrder_max', 'frequencyBand_min', 'frequencyBand_max')

        disp(['The results are saved as: "ResultsOMA_plscf_', currentDate, '.mat"'])
        disp('----------------------------------------------------------------------')
    
    case 'EMA'
   
        save([pathFolderResults, 'EMA_pLSCF_results_' currentDate], 'physicalStablePoles_pLSCF', ...
            'naturalFrequencies_pLSCF', 'dampingRatio_pLSCF', 'complexModes_normalized_pLSCF','modes_normalized_pLSCF', ...
            'modelOrder_min', 'modelOrder_max', 'frequencyBand_min', 'frequencyBand_max')

        disp(['The results are saved as: "ResultsEMA_plscf_', currentDate, '.mat"'])
        disp('----------------------------------------------------------------------')
    
end

disp('Program completed successfully.')
disp('----------------------------------------------------------------------')