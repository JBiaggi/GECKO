modelOG = readCbModel('iYali.xml');

%%
load('iYali4_corr.mat')

% Correct the data type of model.subSystems
var = modelOG.subSystems;

for i = 1:length(model.subSystems)
    var{i} = model.subSystems(i);
end

model.subSystems = var;

% Loop to fix how metabolites are names
for i = 1:numel(model.mets)
    model.mets{i} = [model.mets{i}, model.Comps{i}];
end

% Creates compart
model.comps = modelOG.comps;
model.compNames = modelOG.compNames;

model.comps{1} = 'cy';
model.comps{2} = 'en';
model.comps{3} = 'ex';

model.comps{5} = 'em'; 
model.comps{6} = 'go';


model.comps{9} = 'mi';

model.comps{11} = 'nu';
model.comps{12} = 'pe';
model.comps{13} = 'va';

% Turns model into RAVEN format and fixes some comps names for better GECKO
% compatibility

model = ravenCobraWrapper(model);

model.comps{1} = 'c';
model.comps{3} = 'e';
model.comps{9} = 'm';


% Adds lipid exchange reaction so it can be maximized by FSEOF. susSystem
% must be fixed manually
[model, addedRxns] = addExchangeRxns(model,'out','m1640');
model.subSystems(1923) = model.subSystems(1922);

% Reaction names fixes
model.rxnNames(find(strcmp(model.rxns,'27'))) = {'homoaconitase'};
model.rxnNames(find(strcmp(model.rxns,'117'))) = {'2-methylcitrate dehydratase'};
model.rxnNames(find(strcmp(model.rxns,'iYL0459'))) = {'phosphoribosylglycinamide formyltransferase 1'};

% GPR fixes
    model = changeGeneAssoc(model,'958','YALI0C24101g'); % had another GPR associated to an enzyme with a different function
    model = changeGeneAssoc(model,'336u','YALI0E32769g or YALI0D07986g'); % GPR for DGA2, which is missing in the original model
    model = changeGeneAssoc(model,'336u_1','YALI0E32769g or YALI0D07986g'); % GPR for DGA2, which is missing in the original model
    model = changeGeneAssoc(model,'yli0053','YALI0E32769g or YALI0D07986g'); % GPR for DGA2, which is missing in the original model
    % Ethanol dehydrogenases
    % The reverse rxn (acetaldehyde to ethanol is not as ubiquitous as the
    % forward. GPR kept only for the protein with the greatest homology to
    % S. cerevisae ADH1, shown experimentally to be the only ADH of yeast
    % to metabolize acetaldehyde
    % (https://doi.org/10.1111/j.1567-1364.2011.00760.x)
    model = changeGeneAssoc(model,'2115','YALI0A16379g');

    % Mitochondrial alcohol dehydrogenases
    % Some genes (YALI0D01738g, YALI0C06171g, YALI0E12463g, YALI0A15147g, YALI0F25003g, YALI0F08129g, YALI0B08404g, YALI0E19921g) were removed due to low homology according to Uniprot
    model = changeGeneAssoc(model,'165','YALI0F29623g');

    % Aldehyde dehydrogenases
rxns = {'172';...
    '173';...
    '174';...
    '175';...
    '201';...
    '2116'};

grRules = {'YALI0D07942g or YALI0F04444g';...
    'YALI0C03025g';...
    'YALI0E00264g';...
    'YALI0E00264g';...
    'YALI0D07942g or YALI0F04444g';...
    'YALI0D07942g or YALI0F04444g';...
    };

model = changeGrRules(model,rxns,grRules);

    % Phosphoribosylglycinamide formyltransferase 1 - iYL0459
    % A0A1D8ND08 must be found in the uniprot.tsv file and the gene
    % YALI1D03865g added to it.
    model = changeGeneAssoc(model,'iYL0459','YALI1D03865g');

% BOUNDS FIXES
    % Block glucose uptake because we are working with glycerol
    model = setParam(model,'eq','1714',0);
    % Allow glycerol uptake
    model = setParam(model,'lb','1808', -1000);
    model = setParam(model,'ub','1808', 1000);
    % Glycerol dehydrogenase. It is not native to Yarrowia.
    model = setParam(model,'eq','487', 0);
    % Aldehyde dehydrogenases. The cytosolic ones are propoably unused.
    model = setParam(model,'eq',{'173','2116'}, 0);
    % Isocitratre dehydrogenase is only mitochondrial in Yali
    % (https://doi.org/10.1007/s12010-013-0373-1). The cytosolic and
    % peroxissomal reactions are YeastGEM artifacts.
    model = setParam(model,'eq',{'659','661'}, 0);
    % See reversibility tab
    model = setParam(model,'lb',{'163','165','658'}, 0);

% Reversibility changes
    % The reversibility of alcohol dehydrogenases must be dealt with
    % separate reactions
    model = setParam(model,'rev',{'163','165'},0);
    % Fixing the reversibility of isocitrate dehydrogenase (NAD+)
    model = setParam(model,'rev','658',0);

% Set growth as objective funcion
model = setParam(model,'obj',{'biomass_C'},1);

% Create model ID

%model.id = {'iYali4_corr'};

% Save model
exportModel(model,'iYali4_corr.xml')
