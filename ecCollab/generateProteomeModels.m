%initCobraToolbox(false)

% Get model adapter
adapterLocation = fullfile(findGECKOroot,'ecCollab','ecCollabAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

% Load pooled model
ecModel = loadEcModel('eciYali.yml');

% Load cultivation and lipids-related data
fluxData = loadFluxData();
lipidNchainData = loadLipidNchainData();

% Adjust model to experimental data
%TODO: model name -> ecModel_StrainName
%TODO: Loop through fluxData should start here

% Calculate GANnonPol for later
GAMnonPol = calculateGAM(ecModel);

%% Update the biomass equation for each condition

% Change protein content in the biomass equation
ecModel_1 = scaleBioMassYali4(ecModel, 'protein', fluxData.Ptot(1));

% Change lipid content in the biomass equation
ecModel_1 = scaleBioMassYali4(ecModel_1, 'lipid', lipidNchainData.Ltot(1));

% Adjust fatty acid distribution
ecModel_1 = updateAcylPool(ecModel_1, lipidNchainData, ModelAdapter);

% Balance out mass with carbohydrate content
[X,~,C,~,~,~,~] = sumBioMassYali4(ecModel_1, false);

delta = X - 1;  % difference to balance
fC = (C - delta)/C;
ecModel_1 = rescalePseudoReaction(ecModel_1, 'carbohydrate', fC);

% Recalculate GAM
[~,~,ecModel_1] = calculateGAM(ecModel_1, GAMnonPol,true);