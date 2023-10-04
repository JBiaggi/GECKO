% Script to determine fluxes to be loaded to ecModel

%initCobraToolbox(false)

% Set up the adapter and load the model
adapterLocation = fullfile(findGECKOroot, 'ecYaliGEM', 'ecYaliGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
ecModel = loadEcModel('ecYaliGEM.yml');

rxnID = 'xMAINTENANCE';

% Set the maintenance objective and calculate the maximum NGAM
ecModel = setParam(ecModel, 'obj', rxnID, 1);
sol = solveLP(ecModel);
maxNGAM = -sol.f;

% Set the biomass objective and calculate the initial growth rate
ecModel = setParam(ecModel, 'obj', 'xBIOMASS', 1);
sol = solveLP(ecModel);
initial_growth = -sol.f;

% Define the number of steps and NGAM values
nSteps = 49;
NGAMs = linspace(0, maxNGAM, nSteps + 1);
growthRates = zeros(1, nSteps + 1);
growthRates(1) = initial_growth;

% Iterate over NGAM constraints and calculate growth rates
for i = 1:nSteps
    ecModel = setParam(ecModel, 'lb', rxnID, NGAMs(i + 1));
    sol = solveLP(ecModel);
    
    growthRates(i + 1) = -sol.f;
end

% Plotting NGAM vs Growth Rates
figure;
plot(NGAMs, growthRates, 'o-', 'LineWidth', 2);
xlabel('NGAM Constraint');
ylabel('Growth Rate');
title('NGAM vs Growth Rate');
grid on;

% Adding labels to some key points
text(NGAMs(1), initial_growth, 'Initial Point', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(NGAMs(end), growthRates(end), 'Final Point', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Display the plot
xlim([min(NGAMs), max(NGAMs)]);
ylim([min(growthRates), max(growthRates)]);
title('NGAM vs Growth Rate');
xlabel('NGAM Constraint');
ylabel('Growth Rate');
grid on;

% Save the plot as an SVG file
outputFileName = 'NGAM_vs_GrowthRate.svg';
saveas(gcf, outputFileName, 'svg');

disp(['Plot saved as ' outputFileName]);