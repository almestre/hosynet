%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  NULL MODELS
%% Creation of the null models for the analyses of network modularity
%% Methods:
% Number of null models = 999
% The models will be randomised versions of the model of certain associations in America
% Type of null model = fixed null model (see Krasnov2012)
%% Description of the fixed null model:
% In this model, each random network has the same number of connections per species as the real network. Random networks (pmatrices) were constructed using an independent swap algorithm in which the original matrix is reshuffled by repeatedly swapping 2 # 2 submatrices that preserve the row and column totals (Roberts and Stone 1990). The null model implemented in our study maintained the bipartite structure of the networks by allowing connections between hosts and parasites and not allowing connections among hosts or among parasites.
% MATLAB function swap.m written by B. Semmens
% 999 random matrices for each real network
% We calculate the Q for each of these matrices
% Statistical significance of the modularity of each observed network (p value) is the fraction of the 999 random matrices with a modularity value equal to or higher than the observed one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To see the matrix before importing it
% type interaction_matrix.csv

% Set the random seed for repeatability 
seed = 5485798;
rng(seed);

% Import the matrix labels to matlab
host_labels = readtable ('host_labels.txt','TextType','string');
symbiont_labels = readtable ('symbiont_labels.txt','TextType','string');

% Import the matrix to matlab
A = readmatrix('interaction_matrix.txt');

% Generate the null models and store them in the cell array null_models
A_null = cell(999,1);
for k = 1 : 999
    A_null{k} = swap(A,10000);
end
