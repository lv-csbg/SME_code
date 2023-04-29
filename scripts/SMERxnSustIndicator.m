function rxnImpactArray = SMERxnSustIndicator(model, KeggID, metIndicator)
%SMERxnToImpact Matches metabolite sustainability indicator 
%               to the corresponding exchange reaction
%
%   Input:
%         model = model with KeggID for metabolites
%         KeggID = array of KeggIDs for metabolites which have an indicator 
%                  value 
%         metIndicator = array of metabolite sustainability indicators,
%                   corresponds with KeggID array      
%   Output:
%         rxnEconArray = array of sustainability 
%                        indicators for each reaction of a model

%% Array of metabolite economic impact
% matches metabolites from the array to metabolites in the model using
% KeggIDs
metabolite_indicator_array = zeros(size(model.metKEGGID));

for z = 1:length(metabolite_indicator_array)
    
    if isempty(model.metKEGGID{z}) == 0
        index = find(strcmp(KeggID, model.metKEGGID{z}));
        if index >0
            metabolite_indicator_array(z) = metIndicator(index); 
        end
    end
    
end
%% Make impact of each reaction
clc
exc_reactions = findExcRxns(model);

reaction_indicator_array =zeros(size(model.rxns));
for z = 1:length(exc_reactions)
    
    if exc_reactions(z) == 1
        
        met_temp = findMetsFromRxns(model, model.rxns(z));
        met_temp = met_temp{1};
        index_temp = find(strcmp(model.mets, met_temp));
        
        reaction_indicator_array(z) = metabolite_indicator_array(index_temp);
    end    
end

rxnImpactArray = reaction_indicator_array;
end

