
checkInstallation; % Confirm that RAVEN is functional, should be 2.8.3 or later.

%   - GECKO can be installed via cloning or direct download of ZIP file.
%     See installation instructions in the README.md:
%     https://github.com/SysBioChalmers/GECKO/tree/main#installation
%   - Alternatively, GECKO can be installed as MATLAB Add-On:
%     https://se.mathworks.com/help/matlab/matlab_env/get-add-ons.html
%   - Add the appropriate GECKO (sub)folders to MATLAB path:
GECKOInstaller.install

%% STAGE 1: expansion from a starting metabolic model to an ecModel structure
% STEP 1 Set modelAdapter - WORKS
adapterLocation = fullfile(findGECKOroot,'ecYali4GEM','ecYali4Adapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

% STEP 2 Load conventional iYali - works
model = loadConventionalGEM();
       
% STEP 3-4 Prepare ecModel - I generated a custom uniprot.tsv file where I
% got the KEGG crossref and exchanged it for the gene_oln
[ecModel, noUniprot] = makeEcModel(model,false);

% STEP 5 Store model in YAML format
saveEcModel(ecModel,'ecYali4GEM_stage1.yml');

%% STAGE 2: integration of kcat into the ecModel structure
%ecModel=loadEcModel('ecYali4GEM_stage1.yml'); 

% STEP 6 Fuzzy matching with BRENDA
% Requires EC numbers, which are here first taken from the starting model,
% with the missing ones taken from Uniprot & KEGG databases.
ecModel         = getECfromDatabase(ecModel);

% Do the actual fuzzy matching with BRENDA.
kcatList_fuzzy  = fuzzyKcatMatching(ecModel);

% STEP 7 DLKcat prediction through machine learning
% Requires metabolite SMILES, which are gathered from PubChem.
[ecModel, noSmiles] = findMetSmiles(ecModel);

% DLKcat runs in Python. An input file is written, which is then used by
% DLKcat, while the output file is read back into MATLAB.
writeDLKcatInput(ecModel,[],[],[],[],true);

% runDLKcat will run the DLKcat algorithm via a Docker image. If the
% DLKcat.tsv file already has kcat values, these will all be overwritten.
runDLKcat();
kcatList_DLKcat = readDLKcatOutput(ecModel);

% STEP 8 Combine kcat from BRENDA and DLKcat
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy, 2); % JSB: Changed DLKcat priority

% STEP 9 Take kcatList and populate edModel.ec.kcat
ecModel  = selectKcatValue(ecModel, kcatList_merged);
ecModel  = selectKcatValue(ecModel, kcatList_DLKcat,'max','ifHigher'); % JSB: another way of prioritizing DLKcat

% STEP 10 Apply custom kcat values
[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel);

% STEP 11 Get kcat values across isozymes
ecModel = getKcatAcrossIsozymes(ecModel);

% STEP 12 Get standard kcat
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);

% STEP 13 Apply kcat constraints from ecModel.ec.kcat to ecModel.S
ecModel = applyKcatConstraints(ecModel);

% STEP 14 Set upper bound of protein pool
% The protein pool exchange is constrained by the total protein content
% (Ptot), multiplied by the f-factor (ratio of enzymes/proteins) and the
% sigma-factor (how saturated enzymes are on average: how close to their
% Vmax to they function based on e.g. metabolite concentrations). In 
% modelAdapter Ptot, f- and sigma-factors can all be specified (as rough
% estimates, 0.5 for each of the three parameters is reasonable).
Ptot  = params.Ptot;
f     = params.f;
sigma = params.sigma;

ecModel = setProtPoolSize(ecModel,Ptot,f,sigma);

saveEcModel(ecModel,'ecYali4GEM_stage2.yml');

%% STAGE 3: model tuning
%ecModel=loadEcModel('ecYali4GEM_stage2.yml'); % Uncomment if you want to
%reload model.

% STEP 15 Test maximum growth rate
% Test whether the model is able to reach maximum growth if glycerol uptake
% is unlimited. First set glycerol uptake unconstraint.
ecModel = setParam(ecModel,'eq','1714',0);
ecModel = setParam(ecModel,'lb','1808',-1000);
% And set growth maximization as the objective function.
ecModel = setParam(ecModel,'obj','xBIOMASS',1);
% Run FBA.
sol = solveLP(ecModel,1);
bioRxnIdx = getIndexes(ecModel,params.bioRxn,'rxns');
fprintf('Growth rate: %f /hour.\n', sol.x(bioRxnIdx))

% STEP 17 Sensitivity tuning
ecModel = setProtPoolSize(ecModel);
[ecModel, tunedKcats] = sensitivityTuning(ecModel);

% Inspect the tunedKcats structure in table format.
struct2table(tunedKcats)

% STEP 18 Curate kcat values based on kcat tuning

saveEcModel(ecModel,'ecYali4GEM_stage3.yml');

% This functional ecModel will also be kept in the GECKO GitHub.
saveEcModel(ecModel,'ecYali4GEM.yml');
saveEcModel(ecModel,'ecYali4GEM.xml');

%% STAGE 4 integration of proteomics data into the ecModel.
%ecModel=loadEcModel('ecYali4GEM_stage3.yml'); % Uncomment if you want to
%reload model.

% STEP 19 Load proteomics data and constrain ecModel
protData = loadProtData(1); %Number of replicates, only one experiment.
ecModel = fillEnzConcs(ecModel,protData);
ecModel = constrainEnzConcs(ecModel);

% STEP 20 Update protein pool
% The protein pool reaction will be constraint by the remaining, unmeasured
% enzyme content. This is calculated by subtracting the sum of 
% ecModel.ec.concs from the condition-specific total protein content. The
% latter is stored together with the flux data that otherwise will be used
% in Step 21.
fluxData = loadFluxData();
ecModel = updateProtPool(ecModel,fluxData.Ptot(1));

% STEP 21 Load flux data
% Matching the proteomics sample(s), condition-specific flux data needs to
% be loaded to constrain the model. This was already loaded in Step 20 for
% gathering Ptot, but is repeated here nonetheless. Flux data is read from
% /data/fluxData.tsv.
fluxData = loadFluxData();
ecModel = constrainFluxData(ecModel,fluxData,1,'max','loose'); % Use first condition.
sol = solveLP(ecModel); % To observe if growth was reached.
fprintf('Growth rate that is reached: %f /hour.\n', abs(sol.f))
% Growth rate of 0.1 is by far not reached, flexibilize protein
% concentrations.

% STEP 22 Protein concentrations are flexibilized (increased), until the
% intended growth rate is reached. This is condition-specific, so the
% intended growth rate is gathered from the fluxData structure.
[ecModel, flexProt] = flexibilizeEnzConcs(ecModel,fluxData.grRate(1),10);

% Neither individual protein levels nor total protein pool are limiting
% growth. Test whether the starting model is able to reach 0.1.
model = constrainFluxData(model,fluxData);
sol = solveLP(model)

% It also only reaches 0.0889! So the metabolic network would not be able
% to adhere to all measured constraints. Perhaps there is something
% incorrect with the measurements? Regardless, the ecModel is now able to
% reach about 0.0889, which will be fine for now.
sol = solveLP(ecModel)

% Inspect the flexibilized proteins.
struct2table(flexProt)

% Growth is reached! Let's make sure we store this functional model.
saveEcModel(ecModel,'ecYaliGEM_stage4.yml');

%% STAGE 5: simulation and analysis
% STEP 23 Example of various useful RAVEN functions
% % Set the upper bound of reaction r_0001 to 10.
% ecModel = setParam(ecModel,'ub','r_0001',10);
% % Set the lower bound of reaction r_0001 to 0.
% ecModel = setParam(ecModel,'lb','r_0001',0);
% % Set the objective function to maximize reaction 'r_4041'.
% ecModel = setParam(ecModel,'obj','r_4041',1);
% % Set the objective function to minimize protein usage.
% ecModel = setParam(ecModel,'obj','prot_pool_exchange',1);
% % Perform flux balance analysis (FBA).
% sol = solveLP(ecModel);
% % Perform parsimonious FBA (minimum total flux).
% sol = solveLP(ecModel,1);
% % Inspect exchange fluxes from FBA solution.
% printFluxes(ecModel,sol.x)
% % Inspect all (non-zero) fluxes from FBA solution.
% printFluxes(ecModel,sol.x,false)
% % Export to Excel file (will not contain potential model.ec content).
% exportToExcelFormat(ecModel,'filename.xlsx');

% STEP 24 Simulate Crabtree effect with protein pool

% (Re)load the ecModel without proteomics integration.
ecModel = loadEcModel('ecYaliGEM.yml');

% We will soon run a custom plotCrabtree function that is kept in the code
% subfolder. To run this function we will need to navigate into the folder
% where it is stored, but we will navigate back to the current folder
% afterwards.
currentFolder = pwd;
cd(fullfile(params.path,'code'))
[fluxes, gRate] = plotCrabtree(ecModel);
% fluxes has all the predicted fluxes, while gRate is a vector with the
% corresponding growth rates that were simulated, as visualized on the
% x-axis in the graph.
% The plot will also be saved in the output subfolder.
saveas(gcf,fullfile(params.path,'output','crabtree.pdf'))

% The two graphs show (left:) exchange fluxes from simulations (lines) and
% experiments (circles, from doi:10.1128/AEM.64.11.4226-4233.1998); and
% (right:) the fraction of protein pool that is used by enzymes. In the
% left graph, the y-axis indicates absolute fluxes, so that glucose uptake
% and CO2 excretion both have positive numbers. The model simulation
% demonstrates the Crabtree-effect: at increasing growth rates yeast
% switches from respiration to fermentation, and this occurs when the
% protein pool becomes fully used and thereby limiting. The shift away from
% respiration is most clearly shown by reduced oxygen uptake and increased
% ethanol excretion.

% For comparison, make a similar Crabtree plot for a conventional GEM.
% Set protein pool to infinite, to mimic a conventional GEM.
ecModel_infProt=setProtPoolSize(ecModel,Inf);
plotCrabtree(ecModel_infProt);
saveas(gcf,fullfile(params.path,'output','crabtree_infProt.pdf'))
% It is obvious that no total protein constraint is reached, and Crabtree
% effect is not observed.

% Perform the Crabtree simulation on the pre-Step 17 ecModel (where kcat
% sensitivity tuning has not yet been applied).
ecModel_preTuning = loadEcModel('ecYeastGEM_stage2.yml');
ecModel_preTuning = setParam(ecModel_preTuning,'lb','r_1714',-1000);
plotCrabtree(ecModel_preTuning);
saveas(gcf,fullfile(params.path,'output','crabtree_preStep17.pdf'))
% Without kcat tuning, the model gets constrained too early (at too low
% growth rates), which means that no solutions exist at high growth rates.

% STEP 25 Selecting objective functions
ecModel = setParam(ecModel,'obj',params.bioRxn,1);
sol = solveLP(ecModel)
fprintf('Growth rate that is reached: %f /hour.\n', abs(sol.f))
% Set growth lower bound to 99% of the previous value.
ecModel = setParam(ecModel,'lb',params.bioRxn,0.99*abs(sol.f));
% Minimize protein pool usage. As protein pool exchange is defined in the
% reverse direction (with negative flux), minimization of protein pool
% usage is computationally represented by maximizing the prot_pool_exchange
% reaction.
ecModel = setParam(ecModel,'obj','prot_pool_exchange',1);
sol = solveLP(ecModel)
fprintf('Minimum protein pool usage: %.2f mg/gDCW.\n', abs(sol.f))

% STEP 26 Inspect enzyme usage
% Show the result from the earlier simulation, without mapping to
% non-ecModel.
ecModel = setParam(ecModel, 'obj', 'r_1714', 1);
ecModel = setParam(ecModel, 'lb', params.bioRxn, 0.25);
sol = solveLP(ecModel, 1);
usageData = enzymeUsage(ecModel, sol.x);
usageReport = reportEnzymeUsage(ecModel,usageData);
usageReport.topAbsUsage

% STEP 27 Compare fluxes from ecModel and conventional GEM
sol = solveLP(ecModel);
% Map the ecModel fluxes back to the conventional GEM.
[mappedFlux, enzUsageFlux, usageEnz] = mapRxnsToConv(ecModel, model, sol.x);
% Confirm that mappedFlux is of the same length as model.rxns.
numel(mappedFlux)
numel(model.rxns)

% STEP 28 Perform (ec)FVA
% Perform FVA on a conventional GEM, ecModel, and ecModel plus proteomics
% integration, all under similar exchange flux constraints.
% First make sure that the correct models are loaded.
model = loadConventionalGEM();
ecModel = loadEcModel('ecYaliGEM.yml');
ecModelProt = loadEcModel('ecYaliGEM_stage4.yml');

% As protein model can maximum reach 0.4529, also set this as constrain for
% all models.
fluxData.grRate(1) = 0.4529;

% Apply same constraints on exchange fluxes
model = constrainFluxData(model,fluxData,1,'max','loose');
ecModel = constrainFluxData(ecModel,fluxData,1,'max','loose');
ecModelProt = constrainFluxData(ecModelProt,fluxData,1,'max','loose');

solveLP(model)
solveLP(ecModel)
solveLP(ecModelProt)

% Prepare output structure.
minFlux = zeros(numel(model.rxns),3);
maxFlux = minFlux;

% Run ecFVA for each model.
[minFlux(:,1), maxFlux(:,1)] = ecFVA(model, model);
[minFlux(:,2), maxFlux(:,2)] = ecFVA(ecModel, model);
[minFlux(:,3), maxFlux(:,3)] = ecFVA(ecModelProt, model);

% Write results to output file.
output = [model.rxns, model.rxnNames, num2cell([minFlux(:,1), maxFlux(:,1), ...
    minFlux(:,2), maxFlux(:,2), minFlux(:,3), maxFlux(:,3)])]';
fID = fopen(fullfile(params.path,'output','ecFVA.tsv'),'w');
fprintf(fID,'%s %s %s %s %s %s %s %s\n','rxnIDs', 'rxnNames', 'minFlux', ...
            'maxFlux', 'ec-minFlux', 'ec-maxFlux', 'ecP-minFlux', 'ecP-maxFlux');
fprintf(fID,'%s %s %g %g %g %g %g %g\n',output{:});
fclose(fID);

% Plot ecFVA results and store in output/.
plotEcFVA(minFlux, maxFlux);
saveas(gca, fullfile(params.path,'output','ecFVA.pdf'))

% STEP 29 Compare light and full ecModels
% For a fair comparison of the two types of ecModels, the custom 
% plotlightVSfull function makes a light and full ecModel for yeast-GEM and
% compares their flux distributions at maximum growth rate.
cd(fullfile(findGECKOroot,'tutorials','full_ecModel','code'))
[fluxLight, fluxFull] = plotlightVSfull();

% The ratio between the two ecModel results indicates which reactions have
% a flux that deviates > 0.1% across light vs. full ecModels.
fluxRatio = fluxFull ./ fluxLight;
changedFlux = abs(fluxRatio-1) > 0.001;
model.rxnNames(changedFlux)

% STEP 30 Contextualize ecModels
% This step is exemplified in tutorials/light_ecModel/protocol.m.
