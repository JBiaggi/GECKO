% Script to update the iYali_corr GEM for GECKO usability

% Get model adapter
adapterLocation = fullfile(findGECKOroot,'ecCollab','ecCollabAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

cd ../models/

% Loads model (in RAVEN) to edit the biomass equations
% Make sure to get the latest iYali_corr from the ecModel folder in GECKO
% The stoichiometric coefficients are from biomass_C equation

model = importModel('iYali4_corr_OG.xml');

%% Remove rxns
% There is no evidence that Yarrowia encodes a acyl
% dihydroxyacetonephosphate reductase. There was no GPR for this rxn, and
% this enzyme cannot be found in Yarrowia's Uniprot.
model = removeReactions(model,'iYL0336',false,false,false);

% Fix rxn bounds
% DAG acyltransferase rxn (lb) is not reversible
model = setParam(model,'lb',{'336u'},0);

% DAG acyltransferase rxn (mm) was blocked for no apparent reason
model = setParam(model,'ub',{'336u_1','337u_2'},[1000,1000]);

% CDP-diacylglycerol synthase in this compartment was resulting in TAGs
% being produced in a weird way. Blocking this reaction fixes the issue
% without hampering the production of CDP-diacylglycerol
model = setParam(model,'eq','258',0);


%% Add rxns
% Malic enzyme (NADPH) was improperly removed. Restore the original
% reaction and fix the name of malic enzyme (NADH). For more details, see
% (https://doi.org/10.1007/s10529-013-1302-7)
idx = find(strcmp(model.rxns,'718'));
model.rxnNames(idx) = {'malic enzyme (NAD)'};
model = addReaction(model, '719', {'m538','m178','m6','m44','m176'}, [-1, -1, 1, 1, 1], 'malic enzyme (NADP)');
model = changeGeneAssoc(model,'719','YALI0E18634g');

% Add carbohydrate pseudoreaction
model = addMetabolites(model, 'm1836', 'carbohydrate', 'c', 'Pseudometabolite for biomass equation. Necessary for sumBiomMass function of RAVEN based GEMs');
model = addReaction(model, 'xCARBOHYDRATE', {'m294','m401','m1123','m1324','m1836'}, [-0.00686, -0.868358, -0.943397, -0.235849, 1], 'carbohydrate pseudoreaction');

% Add DNA pseudoreaction
model = addMetabolites(model, 'm1837', 'DNA', 'c', 'Pseudometabolite for biomass equation. Necessary for sumBiomMass function of RAVEN based GEMs');
model = addReaction(model, 'xDNA', {'m89','m459','m465','m505','m1837'}, [-0.01007, -0.010436, -0.009226, -0.010353, 1], 'DNA pseudoreaction');

% Add RNA pseudoreaction
model = addMetabolites(model, 'm1838', 'RNA', 'c', 'Pseudometabolite for biomass equation. Necessary for sumBiomMass function of RAVEN based GEMs');
model = addReaction(model, 'xRNA', {'m86','m93','m95','m149','m1838'}, [-0.055401, -0.050871, -0.05721, -0.059363, 1], 'RNA pseudoreaction');

% Add ion pseudoreaction
model = addMetabolites(model, 'm1839', 'ion', 'c', 'Pseudometabolite for biomass equation. Necessary for sumBiomMass function of RAVEN based GEMs');
model = addReaction(model, 'xION', {'m964','m1839'}, [-0.02, 1], 'ion pseudoreaction');

% Modify protein pseudoreaction
rxns = {'xAMINOACID'};
equations.mets = {'m50','m74','m114','m130','m267','m272','m310','m319','m441','m443',...
                  'm615','m743','m765','m770','m772','m775','m793','m859','m992','m1008','m1726'};
equations.stoichCoeffs = {[-0.243904, -0.044232, -0.567939, -0.243866, -0.186531, -0.509943, -0.125563,...
                            -0.186498, -0.284699, -0.003632, -0.053881, -0.090566, -0.206372, -0.210543,...
                            -0.002154, -0.19455, -0.275005, -0.082571, -0.041282, -0.172773, 1]};

model = changeRxns(model, rxns, equations, 1);

% Change the name to protein
    % get xAMINOACID index
    idx = find(strcmp(model.rxns,rxns));
    
    model.rxns(idx) = {'xPROTEIN'};
    
    model.rxnNames(idx) = {'protein pseudoreaction'};

% Modify lipid pseudoreaction
rxns = {'xLIPID'};
equations.mets = {'m359','m1000','m1631','m1640','m1648','m1700','m1701','m1705','m1727'};
equations.stoichCoeffs = {[-0.003029, -0.035251, -0.000035, -0.000234, -0.000075, -0.000237, -0.000346, -0.000058, 1]};

model = changeRxns(model, rxns, equations, 1);

% Change the name to lipid
    % get xLIPID index
    idx = find(strcmp(model.rxns,rxns));
    
    model.rxnNames(idx) = {'lipid pseudoreaction'};

% Remove old lipid pseudoreaction

model = removeReactions(model,'2108',false,false,false);

% Modify xBIOMASS
rxns = {'xBIOMASS'};
equations.mets = {'m32', 'm141', 'm1836', 'm1837', 'm1838', 'm1839', 'm1726', 'm1727',...
                   'm10', 'm35', 'm143', 'm1401'};
equations.stoichCoeffs = {[-23.09, -23.09, -1, -1, -1, -1, -1, -1, 23.09, 23.09, 23.09, 1]};

model = changeRxns(model, rxns, equations, 1);

% Change the name to biomass
% get xBIOMASS index
idx = find(strcmp(model.rxns,rxns));

model.rxnNames(idx) = {'biomass pseudoreaction'};

% Test changes
sol = solveLP(model, 1);

model = setParam(model, 'obj', 'xBIOMASS', 1);

sol2 = solveLP(model, 1);

% Calculate GAMnonPol

GAMnonPol = calculateGAM(model);

% Fix biomass composition
[X,~,C,~,~,~,~] = sumBioMassYali4(model, false);
delta = X - 1;  % difference to balance
fC = (C - delta) / C;
model = rescalePseudoReaction(model, 'carbohydrate', fC);

% Recalculate and set GAM
[~,~,model] = calculateGAM(model, GAMnonPol,true);

% Save changes
exportModel(model,'iYali.xml');
exportToExcelFormat(model,'iYali.xlsx');

%% Helper function to add metabolites
function model = addMetabolites(model, metID, metName, compartment, metNotes)
    metsToAdd.mets = {metID};
    metsToAdd.metNames = {metName};
    metsToAdd.compartments = {compartment};
    metsToAdd.metNotes = {metNotes};

    model = addMets(model, metsToAdd);
end

% Helper function to add reactions
function model = addReaction(model, rxnID, mets, stoichCoeffs, rxnName)
    rxnsToAdd.rxns = rxnID;
    rxnsToAdd.mets = mets;
    rxnsToAdd.stoichCoeffs = stoichCoeffs;
    rxnsToAdd.rxnNames = {rxnName};
    rxnsToAdd.lb = 0;
    rxnsToAdd.ub = 1000;
    rxnsToAdd.subSystems = {''};

    model = addRxns(model, rxnsToAdd);
end
