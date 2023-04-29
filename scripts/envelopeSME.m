function [outputx outputmin outputmax] = envelopeSME(model, obj_func, biomass_reaction , n)
%Outputs the points of an envelope. The product is calculated by a
%given objective function. Can be used to draw sustainability envelopes
%with a sustainability objective function

%   Inputs:
%       model = COBRA model structure
%       obj_func = objective function to use for optimization, where a
%                  value for each model reaction is given
%       biomass_reaction = vector the same length as model.rxns, where the
%                          biomass reaction has a value 1
%       n = number of points to plot
%   Outputs:
%       outputx = vector of biomass production values
%       outputmin = vector of minimum objective values at different
%                   growth-rates
%       outputmax = vector of maximum objective values at different
%                   growth-rates
%%
    model.c = biomass_reaction;
    max_biomass = optimizeCbModel(model, 'max').f;
    min_biomass = optimizeCbModel(model,'min').f;
    model.c = obj_func;       
    x =[];
    maxList = [];
    minList = [];
    biomass_space = linspace(min_biomass, max_biomass,n);
    for m = 1:n
        
        biomass_bound = biomass_space(m);
        model.lb(biomass_reaction==1) = biomass_bound;    
                
        max = optimizeCbModel(model, 'max'); 
        min = optimizeCbModel(model, 'min');             
        
        minList = [minList;min.f];
        maxList = [maxList; max.f];
        
        x = [x;biomass_bound];
               
        model.c = obj_func;
        
    end
    
    outputx = x;
    outputmin = minList;
    outputmax = maxList;
end

