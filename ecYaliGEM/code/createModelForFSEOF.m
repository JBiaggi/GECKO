%% This script will add a TAG exchange reaction to iYali before it us turned into an ecModel so it can be used for FSEOF

adapterLocation = fullfile(findGECKOroot,'ecYaliGEM','ecYaliGEMAdapter.m');
ModelAdapter = ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapter.getParameters();
model = loadConventionalGEM();
model = setParam(model,'eq','y001714',0);

[model, addedRxns] = addExchangeRxns(model,'out','m1640');

% rxns with wrong name in the original iYali - JSB
model.rxnNames(find(strcmp(model.rxns,'y000027'))) = {'homoaconitase'};
model.rxnNames(find(strcmp(model.rxns,'y000117'))) = {'2-methylcitrate dehydratase'};
model = setParam(model,'rev','y001808',1);
model = changeGeneAssoc(model,'y000958','YALI0C24101g'); % had another GPR associated to an enzyme with a different function

%Make irreversible with RAVEN function:
model = convertToIrrev(model);

%Analysis will be around the experimental biomass yield:
Ysx = 0.456;    %gDW/g(gluc)
Ysx = Ysx*92.05;  %gDW/mol(gluc)
Ysx = Ysx/1000; %gDW/mmol(gluc)

%Simulation:
results.model   = model;
results.glycerol = compare_substrate(model,'triglyceride exchange (OUT)','glycerol',Ysx);

% Initialize a cell array to store gene names
results.glycerol.geneNames = cell(size(results.glycerol.genes));

% Loop through each gene and find the corresponding name in results.glycerol.rxns
for i = 1:numel(results.glycerol.genes)
    gene = results.glycerol.genes{i};
    
    % Find rows in results.glycerol.rxns where the gene identifier partially matches
    idx = cellfun(@(x) contains(x, gene), results.glycerol.rxns(:, 3));
    
    if any(idx)
        % Retrieve the corresponding name from column 2 for the first partial match found
        results.glycerol.geneNames{i} = results.glycerol.rxns{idx, 2};
    else
        results.glycerol.geneNames{i} = 'Gene not found';  % Set a placeholder if gene is not found
    end
end

%% Functions
function FC = compare_substrate(model,product,substrate,Ysx)

%Simulate WT (100% growth):
FC.flux_WT = simulateGrowth(model,product,substrate,1);

%Simulate forced (X% growth and the rest towards product) based on yield:
posX     = strcmp(model.rxnNames,'Biomass production');
alphaExp = Ysx/FC.flux_WT(posX);
%alpha    = (alphaExp/2):(alphaExp/10):(2*alphaExp);
alpha    = (alphaExp/2):(alphaExp/30):(alphaExp);
FC.alpha = alpha;
v_matrix = zeros(length(model.rxns),length(alpha));
k_matrix = zeros(length(model.rxns),length(alpha));
for i = 1:length(alpha)
    FC.flux_MAX   = simulateGrowth(model,product,substrate,alpha(i));
    v_matrix(:,i) = FC.flux_MAX;
    k_matrix(:,i) = FC.flux_MAX./FC.flux_WT;
end

%Generate rxn equations:
rxnEqs = printRxnFormula(model,model.rxns,false,true,true);

%Take out rxns with no grRule:
withGR   = ~cellfun(@isempty,model.grRules);
v_matrix = v_matrix(withGR,:);
k_matrix = k_matrix(withGR,:);
gene_rxn = model.rxnGeneMat(withGR,:);
FC.rxns  = [model.rxns(withGR) model.rxnNames(withGR) model.grRules(withGR) rxnEqs(withGR)];

%Filter out rxns that are always zero -> k=0/0=NaN:
non_nan  = sum(~isnan(k_matrix),2) > 0;
v_matrix = v_matrix(non_nan,:);
k_matrix = k_matrix(non_nan,:);
gene_rxn = gene_rxn(non_nan,:);
FC.rxns  = FC.rxns(non_nan,:);

%Replace remaining NaNs with 1s:
k_matrix(isnan(k_matrix)) = 1;

%Replace any Inf value with 1000 (maximum value is ~700):
k_matrix(isinf(k_matrix)) = 1000;

%Filter out values that are inconsistent at different alphas:
always_down  = sum(k_matrix <= 1,2) == length(alpha);
always_up    = sum(k_matrix >= 1,2) == length(alpha);
incons_rxns  = always_down + always_up == 0;
incons_genes = sum(gene_rxn(incons_rxns,:),1) > 0;
incons_rxns  = sum(gene_rxn(:,incons_genes),2) > 0;
v_matrix     = v_matrix(~incons_rxns,:);
k_matrix     = k_matrix(~incons_rxns,:);
gene_rxn     = gene_rxn(~incons_rxns,:);
FC.rxns      = FC.rxns(~incons_rxns,:);

%Order from highest to lowest k:
FC.k_rxns   = mean(k_matrix,2);
[~,order]   = sort(FC.k_rxns,'descend');
FC.k_rxns   = FC.k_rxns(order,:);
FC.v_matrix = v_matrix(order,:);
FC.k_matrix = k_matrix(order,:);
gene_rxn    = gene_rxn(order,:);
FC.rxns     = FC.rxns(order,:);

%Create list of remaining genes and filter out any inconsistent score:
FC.genes     = model.genes(sum(gene_rxn,1) > 0);
FC.geneNames = model.geneShortNames(sum(gene_rxn,1) > 0);
FC.k_genes   = zeros(size(FC.genes));
gene_rxn     = gene_rxn(:,sum(gene_rxn,1) > 0);
cons_genes   = false(size(FC.genes));
for i = 1:length(FC.genes)
    k_set         = FC.k_rxns(gene_rxn(:,i) > 0);
    always_down   = sum(k_set <= 1) == length(k_set);
    always_up     = sum(k_set >= 1) == length(k_set);
    cons_genes(i) = always_down + always_up == 1;
    FC.k_genes(i) = mean(k_set);
end
FC.genes     = FC.genes(cons_genes);
FC.geneNames = FC.geneNames(cons_genes);
FC.k_genes   = FC.k_genes(cons_genes);

%Filter any value between mean(alpha) and 1:
unchanged    = (FC.k_genes >= mean(alpha) - 1e-3) + (FC.k_genes <= 1 + 1e-3) == 2;
FC.genes     = FC.genes(~unchanged);
FC.geneNames = FC.geneNames(~unchanged);
FC.k_genes   = FC.k_genes(~unchanged);

%Order from highest to lowest k:
[~,order]    = sort(FC.k_genes,'descend');
FC.genes     = FC.genes(order,:);
FC.geneNames = FC.geneNames(order,:);
FC.k_genes   = FC.k_genes(order,:);

% %% Create cell array for genes names
% % Initialize a cell array to store gene names
% FC.geneNames = cell(size(genes));
% 
% % Loop through each gene and find the first corresponding name in results.glycerol.rxns
% for i = 1:numel(genes)
%     gene = genes{i};
% 
%     % Find the first row in results.glycerol.rxns where the gene identifier matches
%     idx = find(strcmp(results.glycerol.rxns(:, 3), gene), 1);
% 
%     if ~isempty(idx)
%         % Retrieve the corresponding name from column 2
%         geneNames{i} = results.glycerol.rxns{idx, 2};
%     else
%         geneNames{i} = 'Gene not found';  % Set a placeholder if gene is not found
%     end
% end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flux = simulateGrowth(model,rxn,substrate,alpha)

if strcmp(substrate,'glycerol')
    model.ub(strcmp(model.rxnNames,'glycerol exchange'))              = 0;
    model.lb(strcmp(model.rxnNames,'glycerol exchange'))              = 0;
    model.ub(strcmp(model.rxnNames,'glycerol exchange (reversible)')) = 1;
end

%Positions of biomass & target rxn:
posX = strcmp(model.rxnNames,'Biomass production');
posP = strcmp(model.rxnNames,rxn);

%Max growth:
sol = optModel(model,posX,+1);

%Fix growth suboptimal and and minimize fluxes
model.lb(posX) = sol.x(posX)*0.999*alpha;
sol            = optModel(model,1:length(model.rxns),-1);

%Max product:
model.lb(1:length(model.rxns)) = sol.x(1:length(model.rxns));
sol            = optModel(model,posP,+1);
flux           = sol.x;

%Fix also product and minimize fluxes:
model.lb(posP) = sol.x(posP)*0.999;
sol            = optModel(model,1:length(model.rxns),-1);
flux           = sol.x;

% %Fix growth suboptimal and then max product:
% model.lb(posX) = sol.x(posX)*0.999*alpha;
% sol            = optModel(model,posP,+1);
% 
% %Fix also product and minimize fluxes:
% model.lb(posP) = sol.x(posP)*0.999;
% sol            = optModel(model,1:length(model.rxns),-1);
% flux           = sol.x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sol = optModel(model,pos,c,base_sol)

if length(pos) == 1
    rxn_code = model.rxns{pos};
    if strcmp(rxn_code(end-3:end),'_REV')
        rev_pos = strcmp(model.rxns,rxn_code(1:end-4));
    else
        rev_pos = strcmp(model.rxns,[rxn_code '_REV']);
    end
    model.lb(rev_pos) = 0;
    model.ub(rev_pos) = 0;
end

model.c      = zeros(size(model.rxns));
model.c(pos) = c;
sol          = optimizeCbModel(model);

if isempty(sol.x) && length(pos) == 1
    model.lb(rev_pos) = base_sol.x(rev_pos);
    model.ub(rev_pos) = base_sol.x(rev_pos);
    sol               = optimizeCbModel(model);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%