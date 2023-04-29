clear
clc
model = readCbModel('iML1515.xml');
price_table = readtable('metab_econ.xlsx'); 
env_table = readtable('metab_env.xlsx');
soc_table = readtable('metab_soc.xlsx');
%% Make the model anaeorobic
oxygen_bound = 0;

oxygen = ismember(model.rxns,'EX_o2_e','rows' );
model.ub(oxygen==1) = oxygen_bound;
model.lb(oxygen ==1) =-oxygen_bound;

%% Essential genes
% determine the essential genes from Baba et al. 2006
clc
essential = readtable('essential.xlsx');
essential_gene_loc = ismember(model.proteins, essential.Gene);
geneList =[];
for z = 1:length(essential_gene_loc)
    if essential_gene_loc(z) == 0
        geneList = [geneList; model.geneNames(z)];
    end
end

%% Create arrays of metabolite environmental and economic indicators
%The indicators are previously calculated in an excel file
clc
metabolite_econ_array = price_table.EconImpact1e_3USD_mmol;
metabolite_env_array = env_table.EnvImpact1e_3_USD_mmol_;
metabolite_soc_array = soc_table.SocialIndicator1e_3USD_mmol;


%% Generate objective function
clc
objFuncModel = model;

%match metabolite sustainability indicators with their corresponding
%exchange reactions
reaction_economic_array = SMERxnSustIndicator(objFuncModel, price_table.KeggID, metabolite_econ_array);
reaction_env_array = SMERxnSustIndicator(objFuncModel, env_table.KeggID, metabolite_env_array);
reaction_soc_array = SMERxnSustIndicator(objFuncModel, soc_table.KeggID, metabolite_soc_array);

envObj = -reaction_env_array;
%replace Nan values with 0
for z = 1:length(envObj)
    if isnan(envObj(z))
        envObj(z) = 0;
    end
end

econObj = reaction_economic_array;
%replace Nan values with 0
for z = 1:length(econObj)
    if isnan(econObj(z))
        econObj(z) = 0;
    end
end

socObj = reaction_soc_array;
%replace Nan values with 0
for z = 1:length(socObj)
    if isnan(socObj(z))
        socObj(z) = 0;
    end
end

objFunc = (econObj+ envObj + socObj);
%replace Nan values with 0
for z = 1:length(objFunc)
    if isnan(objFunc(z))
        objFunc(z) = 0;
    end
end

%% Show metabolites that are in the objective function

clc
tableEc = [];
tableEnv = [];
tableSoc = [];
tableObj = [];
tableMet = [];
tableKegg =[];
for z = 1:length(model.c)
    if objFunc(z) ~=0
        tableMet = [tableMet;model.rxns(z)];
        %tableKegg = [tableKegg ; model.metKEGGID(x)];
        tableEc = [tableEc;econObj(z)];
        tableEnv = [tableEnv;envObj(z)];
        tableObj = [tableObj;objFunc(z)];
        tableSoc = [tableSoc;socObj(z)];
    end
end
% shows metabolites and their corresponding economic, environmental, social
% and sustainability indicators
objective_table = table(tableMet,  tableEc, tableEnv,tableSoc, tableObj)



%% Genetic algorithms
clc
% objFunc is inserted in the modified optGeneFitnessTilt function using a
% global variable
global objFuncGA
objFuncGA = objFunc;

targetRxn = 'EX_succ_e'; %not taken into account in the modified optGene
substrateRxn = 'EX_glc__D_e';

% uses modified version of optGene -> optGene1
% uses geneList where only non-essential genes are allowed to delete
[z, population, scores, optGeneSol, hashtable] = optGene1(model, ...
    targetRxn, substrateRxn, geneList,'MaxKOs',7)

%% Best design envelope
clc
knockModel = deleteModelGenes(model, optGeneSol.geneList);

%knockModel = deleteModelGenes(model,{'b2415', 'b3737', 'b1380'}) %design #1

[z,minProd,maxProd]=envelopeSME(knockModel, objFunc, knockModel.c , 20);
figure(2)
hold on
plot([z;flip(z)],[minProd;flip(maxProd)],'r','LineWidth',2);
[z,minProd,maxProd]=envelopeSME(model, objFunc, model.c , 20);
plot([z;flip(z)],[minProd;flip(maxProd)],'b','LineWidth',2);

xlabel('Biomass (1/h)')
ylabel(['ISS 1e-4*USD/gDW/h'])
legend({'Design','Wild-type'})
hold off
%% Show exchange reactions active of a design at max growth-rate
clc
model1= knockModel;
%model1 = model;
model1.c = model.c;
max_growth_rate = optimizeCbModel(model1).f
model1.lb(model1.c ==1 )=optimizeCbModel(model1).f;
model1.c = objFunc;
s=optimizeCbModel(model1,'min');
exchange_reactions= findExcRxns(knockModel);
c_array = [];% values in the objective coefficient array
rxn_array = [];% reactions
v_array = [];% flux 
obj_array = [];% sustainability indicator coefficient value
econ_array = [];% economic indicator
env_array = [];% environemntal indicator
soc_array =[];% social indicator
for z = 1:length(objFunc)
    if exchange_reactions(z) ~= 0 && abs(s.x(z)) >1e-3 && objFunc(z) ~=0
    
        c_array = [c_array; objFunc(z)];
        rxn_array = [rxn_array; model1.rxns(z)];
        v_array = [v_array; s.x(z)];
        obj_array = [obj_array; s.x(z).*objFunc(z)];
        econ_array = [econ_array; s.x(z).*econObj(z)];
        env_array = [env_array; s.x(z).*envObj(z)];
        soc_array = [soc_array; s.x(z).*socObj(z)];
    end
end

exc_reaction_table = table(rxn_array, v_array, obj_array, c_array, econ_array, env_array, soc_array)
%% Analysis of top designs from optGene
clc
% takes the hashtable returned from optGene and ranks the designs
clc
allKeys = arrayfun(@char, hashtable.keySet.toArray, 'UniformOutput', false);
allValues = cellfun(@(x) hashtable.get(x), allKeys, 'UniformOutput', false);

design_start = 1;%start from which rank
design_end = 1000; %end with which rank
designs_n = design_end- design_start;

allValues = cell2mat(allValues);
[~,p] = sort(allValues,'ascend');
r = 1:length(allValues);
r(p) = r;

gene_list_design_rank = [];
value_array = [];
deletion_nr_array = [];
substrate = ismember(model.rxns,'EX_glc__D_e');
for n= design_start:design_end
    loc = ismember(r, n);
    temp_design = allKeys{loc == 1};
    temp_gene_list = [];
    temp_deletion_nr = 0;
    for y = 1:length(temp_design)
        
        if temp_design(y) == '1'
            temp_deletion_nr = temp_deletion_nr + 1;
            temp_gene_list = [temp_gene_list; geneList(y)];
        end
    end
    gene_list_design_rank = [gene_list_design_rank; {temp_gene_list}];
    value_array = [value_array; allValues(loc == 1);];
    deletion_nr_array = [deletion_nr_array; temp_deletion_nr];
end

product_v_array = zeros(designs_n, 1);
product_obj_array = zeros(designs_n, 1);
product_array = cell(designs_n, 1);
substrate_v_array = zeros(designs_n, 1);
rank_array = zeros(designs_n, 1);
obj_value_array = zeros(designs_n, 1);
obj_zeroGR_array = zeros(designs_n, 1);
env_array = zeros(designs_n, 1);
econ_array = zeros(designs_n, 1);
soc_array = zeros(designs_n, 1);
deletion_array = cell(designs_n, 1);
GR_array = zeros(designs_n, 1);
rxn_array =cell(designs_n, 1);

for n=1:designs_n+1
    [temp_model, hasEffect, constrRxnNames, deletedGenes] = deleteModelGenes(model, gene_list_design_rank{n});
    constrRxnNames = sort(constrRxnNames);
    temp_optCb = optimizeCbModel(temp_model);
    
    if string(temp_optCb.origStat) == 'INFEASIBLE'
        continue
    end
    
    maxGR = temp_optCb.f;
    temp_model.lb(model.c == 1) = maxGR;
    temp_model.ub(model.c == 1) = maxGR;
    
    temp_model.c = objFunc;
    temp_opt = optimizeCbModel(temp_model,'min');
    obj_result = objFunc.*temp_opt.x;
    
    temp_model.lb(model.c == 1) = 0;
    temp_model.ub(model.c == 1) = 0;
    temp_opt_zero = optimizeCbModel(temp_model, 'min');
    obj_at_zero_GR = temp_opt_zero.f;
    econ_result = sum(econObj.*temp_opt.x);
    env_result = sum(envObj.*temp_opt.x);
    soc_result = sum(socObj.*temp_opt.x);
    
    [temp_obj_max, temp_product_index] = max(obj_result);
    temp_product = model.rxns(temp_product_index);
    temp_product_v = temp_opt.x(temp_product_index);
    substrate_v = temp_opt.x(substrate==1);
    
    GR_array(n) = maxGR;
    product_array (n)= temp_product;
    product_v_array(n) = temp_product_v;
    rank_array(n) =  n+design_start-1;
    obj_value_array(n) = temp_opt.f;
    product_obj_array(n) = temp_obj_max;
    substrate_v_array(n) = substrate_v;
    env_array(n) = env_result;
    econ_array(n) = econ_result;
    soc_array(n) = soc_result;
    obj_zeroGR_array(n) = obj_at_zero_GR;
    deletion_array{n} = string(strjoin(gene_list_design_rank{n}, ', '));
    rxn_array{n} = string(strjoin(constrRxnNames,', '));
    
end
ranking = table(rank_array,product_array, product_v_array, substrate_v_array,GR_array,...
    deletion_nr_array, obj_value_array,obj_zeroGR_array,econ_array, env_array, soc_array, deletion_array,rxn_array)


