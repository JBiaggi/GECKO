%% This script will add a lipids exchange reaction to iYali before it us turned into an ecModel so it can be used for FSEOF

adapterLocation = fullfile(findGECKOroot,'ecYaliGEM','ecYaliGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();
model = loadConventionalGEM();
model = setParam(model,'eq','y001714',0);

[model, addedRxns] = addExchangeRxns(model,'out','m1727');