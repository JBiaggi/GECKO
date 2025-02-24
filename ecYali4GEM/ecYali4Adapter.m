classdef ecYali4Adapter < ModelAdapter
    methods
        function obj = ecYali4Adapter()
            % Set initial values of the obj.params - they can be changed by the user
            
            % Directory where all model-specific files and scripts are kept.
            % Is assumed to follow the GECKO-defined folder structure.
            obj.params.path = fullfile('C:\Users\jusabe\Documents\GECKO\ecYali4GEM');

			% Path to the conventional GEM that this ecModel will be based on.
			obj.params.convGEM = fullfile(obj.params.path,'models','iYali4_corr.xml');

			% Average enzyme saturation factor
			obj.params.sigma = 0.49;

			% Total protein content in the cell [g protein/gDw]
			obj.params.Ptot = 0.202;

			% Fraction of enzymes in the model [g enzyme/g protein]
			obj.params.f = 0.5;
            
            % Growth rate the model should be able to reach when not
            % constraint by nutrient uptake (e.g. max growth rate) [1/h]
			obj.params.gR_exp = 0.455;

			% Provide your organism scientific name
			obj.params.org_name = 'yarrowia lipolytica';
            
            % Taxonomic identifier for Complex Portal
            obj.params.complex.taxonomicID = 'Yarrowia lipolytica';

			% Provide your organism KEGG ID, selected at
			% https://www.genome.jp/kegg/catalog/org_list.html
			obj.params.kegg.ID = 'yli';
            % Field for KEGG gene identifier; should match the gene
            % identifiers used in the model. With 'kegg', it takes the
            % default KEGG Entry identifier (for example YER023W here:
            % https://www.genome.jp/dbget-bin/www_bget?sce:YER023W).
            % Alternatively, gene identifiers from the "Other DBs" section
            % of the KEGG page can be selected. For example "NCBI-GeneID",
            % "UniProt", or "Ensembl". Not all DB entries are available for
            % all organisms and/or genes.
            obj.params.kegg.geneID = 'kegg';

			% Provide what identifier should be used to query UniProt.
            % Select proteome IDs at https://www.uniprot.org/proteomes/
            % or taxonomy IDs at https://www.uniprot.org/taxonomy.
            obj.params.uniprot.type = 'taxonomy'; % 'proteome' or 'taxonomy'
			obj.params.uniprot.ID = '4952'; % should match the ID type
            % Field for Uniprot gene ID - should match the gene ids used in the 
            % model. It should be one of the "Returned Field" entries under
            % "Names & Taxonomy" at this page: https://www.uniprot.org/help/return_fields
            obj.params.uniprot.geneIDfield = 'gene_oln';
            % Whether only reviewed data from UniProt should be considered.
            % Reviewed data has highest confidence, but coverage might be (very)
            % low for non-model organisms
            obj.params.uniprot.reviewed = false;

			% Reaction ID for glucose exchange reaction (or other preferred carbon source)
			obj.params.c_source = '1808'; 

			% Reaction ID for biomass pseudoreaction
			obj.params.bioRxn = 'xBIOMASS';

			% Compartment name in which the added enzymes should be located
			obj.params.enzyme_comp = 'cytoplasm';
        end
		
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
            % Indicates how spontaneous reactions are identified. Here it
            % is done by the reaction have 'spontaneous' in its name.
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
	end
end