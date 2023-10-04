% Define the adapter location
adapterLocation = fullfile(findGECKOroot, 'ecYaliGEM', 'ecYaliGEMAdapter.m');

% Set the model adapter and get parameters
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

% Load the ecModel
ecModel = loadEcModel('ecYaliGEM.yml');
Yxs = 0.463;
exprGrowth = 0.03;
CsourceMM = 92.05;
uptakeCsource = (exprGrowth/Yxs)*1000/CsourceMM;

% Get relevant rxn indexes
poolIdx      = strcmpi(ecModel.rxns, 'prot_pool_exchange');

% Fix experimental growth and carbon source uptake
ecModel = setParam(ecModel,'eq',{'xBIOMASS', 'y001808'}, [exprGrowth -uptakeCsource]);

% Minimize protein pool
ecModel = setParam(ecModel, 'obj', 'prot_pool_exchange', 1);
sol   = solveLP(ecModel);
ecModel = setParam(ecModel, 'lb', 'prot_pool_exchange', sol.x(poolIdx) * 1.01);

% Maximize target, store original ecModel
ecModel = setParam(ecModel, 'obj', 'EXC_OUT_m1640', 1);
ecModel_OG = ecModel;
sol = solveLP(ecModel);

% Calculate lipOG
lipOG = -sol.f;

% Set Qa
Qa = 1.01;

% Initialize FCCs
FCCs = zeros(numel(ecModel.ec.kcat), 1);

progressbar('Flux control coefficients calculation')
for i = 1:numel(ecModel.ec.kcat)
    % Calculate new kcat
    ecModel.ec.kcat(i) = ecModel.ec.kcat(i) * Qa;
    
    % Apply kcat constraints and solve
    ecModel = applyKcatConstraints(ecModel);
    sol = solveLP(ecModel, 1);
    Qlip = -sol.f;
    
    % Calculate FCC
    kcat_i = ecModel_OG.ec.kcat(i);
    FCCs(i) = (Qlip - lipOG) * kcat_i / (lipOG * (kcat_i * Qa - kcat_i));
    
    % Restore the original ecModel
    ecModel = ecModel_OG;
    progressbar(i / numel(ecModel.ec.kcat))
end
progressbar(1)

% Create output variable 
result.FCCs = FCCs;
result.kcat = ecModel_OG.ec.kcat;
result.rxns = ecModel_OG.ec.rxns;

% Sort the results in descending order based on FCCs
[result.FCCs, sortOrder] = sort(result.FCCs,'descend');
result.kcat = result.kcat(sortOrder);
result.rxns = result.rxns(sortOrder);
result.rxnEnzMat = ecModel.ec.rxnEnzMat(sortOrder, :);

% Initialize rxnNames and enzymes as cell arrays
result.rxnNames = cell(length(result.rxns), 1);
result.enzymes = cell(length(result.rxns), 1);

% Populate rxnNames and enzymes
for i = 1:length(result.rxns)
    result.rxnNames(i,1) = ecModel.rxnNames(strcmp(ecModel.rxns,result.rxns(i)));
    if i < length(ecModel.ec.enzymes)
    result.enzymes{i} = ecModel.ec.enzymes(logical(result.rxnEnzMat(i,:))');
    if length(result.enzymes{i}) > 1
        result.enzymes{i} = strjoin(result.enzymes{i},', ');
    end
    end
end

%% Make chart
% Sample data (replace with your data)
FCCs = result.FCCs(1:10); % Top 10 FCCs
rxnNames = result.rxnNames(1:10); % Top 10 rxnNames

% Define a set of pleasant colors
pleasantColors = [0.2 0.4 0.8; 0.8 0.2 0.4; 0.4 0.6 0.2; 0.6 0.2 0.6; 0.2 0.6 0.6;
                  0.8 0.4 0.2; 0.4 0.2 0.8; 0.6 0.6 0.2; 0.2 0.8 0.4; 0.4 0.4 0.4];

% Create a horizontal bar chart
h = barh(FCCs);

% Customize the chart
xlabel('Flux Control Coefficients')

% Set custom labels for the y-axis ticks
yticks(1:10) % Set the number of ticks to match the number of bars
yticklabels(rxnNames) % Set the rxnNames as labels

% Set the chart area size to 5 in. tall by 6.75 in. wide
set(gcf, 'Position', [100, 100, 2 * 6.75 * 100, 5 * 100]) % [left, bottom, width, height]

% Invert the y-axis to display the highest value at the top
set(gca, 'YDir', 'reverse')

% Remove grid lines
grid off

% Set the font size for the axes labels and ticks
set(gca, 'FontSize', 14) % Set font size to 14

% Change axes line color to black and set to 1 pt thick
ax = gca;
ax.XAxis.Color = 'black';
ax.YAxis.Color = 'black';
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;

% Set major tick marks to "cross" and minor tick marks to "outside"
ax.XAxis.TickDirection = 'both';
ax.YAxis.TickDirection = 'both';

% Set the plot area to have a solid black border (1 pt thick) and no fill
set(gca, 'Box', 'on', 'LineWidth', 1, 'Color', 'none')

% Format the data series with pleasant colors
for i = 1:numel(h)
    set(h(i), 'FaceColor', pleasantColors(i, :), 'EdgeColor', 'none') % Use pleasant colors
end

% Save the chart as an SVG file
saveas(gcf, 'FCC_ecYali_lipids.svg', 'svg')

%% Make output file
% Create a table from the result variable
resultTable = table(result.FCCs, result.kcat, result.rxns, result.rxnNames, result.enzymes,  'VariableNames', {'FCCs', 'kcat', 'rxns', 'rxnNames', 'enzymes'});

% Define the file path and name
filePath = fullfile(params.path, 'output', 'ecYali_FCCs_lipids.csv');

% Save the result as a CSV file to the specified folder
writetable(resultTable, filePath);