% Script to determine fluxes to be loaded to ecModel

%initCobraToolbox(false)

% Set up the adapter and load the model
adapterLocation = fullfile(findGECKOroot, 'ecYaliGEM', 'ecYaliGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

ecModel = loadEcModel('ecYaliGEM_stage4.yml');

rxnID = 'prot_pool_exchange';

% Set the biomass objective and calculate the initial growth rate
ecModel = setParam(ecModel, 'obj', 'xBIOMASS', 1);
sol = solveLP(ecModel);
initial_growth = -sol.f;

% Define the number of steps and prot_pool values
nSteps = 49;
poolSize = linspace(ecModel.lb(strcmp(ecModel.rxns, 'prot_pool_exchange')), 0, nSteps + 1);
growthRates = zeros(1, nSteps + 1);
growthRates(1) = initial_growth;

% Iterate over prot_pool constraints and calculate growth rates
for i = 1:nSteps
    ecModel = setParam(ecModel, 'lb', rxnID, poolSize(i + 1));
    sol = solveLP(ecModel);
    growthRates(i + 1) = -sol.f;
end

% Plotting prot_pool vs Growth Rates
figure;
plot(poolSize, growthRates, 'o-', 'LineWidth', 2);
xlabel('Protein pool size');
ylabel('Growth Rate');
title('Protein pool size vs Growth Rate');
grid on;

% Adding labels to some key points
text(poolSize(1), initial_growth, 'Initial Point', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(poolSize(end), growthRates(end), 'Final Point', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% Display the plot
xlim([min(poolSize), max(poolSize)]);
ylim([min(growthRates), max(growthRates)]);
title('Protein pool size vs Growth Rate');
xlabel('Protein pool size');
ylabel('Growth Rate');
grid on;

% Save the plot as an SVG file
outputFileName = 'protPool_vs_GrowthRate.svg';
saveas(gcf, outputFileName, 'svg');

disp(['Plot saved as ' outputFileName]);