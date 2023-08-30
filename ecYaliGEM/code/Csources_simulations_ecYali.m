adapterLocation = fullfile(findGECKOroot,'ecYaliGEM','ecYaliGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();

ecModel = loadEcModel('ecYaliGEM.yml');

%% Change media batch part

flux = +1000;

c_source = 'glycerol exchange';

pos = find(strcmp(ecModel.rxnNames,c_source));

% First block any uptake

[rxnIDs,exchange] = getExchangeRxns(ecModel);
%Exclude protein pool from exchange reactions list
protIndex = find(contains(ecModel.rxnNames,'prot_'));
exchange  = setdiff(exchange,protIndex);
%First allow any exchange (uptakes and secretions)
ecModel.ub(exchange) = +1000;
%Then block all uptakes
uptakes            = exchange(find(contains(rxnIDs,'_REV')));
model.ub(uptakes)  = 0;