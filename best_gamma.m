% Selection of the best gamma parameter (i.e. resolution) for the modularity
    % 1. Create a set of partitions covering a representative range of gamma values
    % 2. Identify the partition with the largest domain of optimisation (using the CHAMP package)

% Set the random seed for repeatability 
seed = 6429819;
rng(seed);

% Vector with the gamma values
gamma_values = 0:0.01:6;

% Create the containers to store the data generated in the for loop
B_gtest = cell(size(gamma_values,2),1);
twom_gtest = zeros(size(gamma_values,2),1);
Q_gtest = zeros(size(gamma_values,2),1);
S_gtest = cell(size(gamma_values,2),1);
n_it_gtest =  zeros(size(gamma_values,2),1);

for k1 = 1 : size(gamma_values,2)*100
    if k1 == 1
        k2 = 1;
    end
    
    % Build modularity matrix and null model from the bipartite network
    [B_gtest_temp,twom_gtest_temp]=bipartite(A,gamma_values(k2));
    % Apply the genlouvain algorithm to the Barber modularity matrix
    [S_gtest_temp,Q_gtest_temp,n_it_gtest_temp]=iterated_genlouvain(B_gtest_temp);
    Q_gtest_temp=Q_gtest_temp/twom_gtest_temp;
    
    if Q_gtest(k2) < Q_gtest_temp
        B_gtest{k2} = B_gtest_temp;
        twom_gtest(k2) = twom_gtest_temp;
        Q_gtest(k2) = Q_gtest_temp;
        S_gtest{k2} = S_gtest_temp;
        n_it_gtest(k2) =  n_it_gtest_temp;
    end
    
    if mod(k1,100) == 0
        k2 = k2 + 1;
    end
  
end

% Number of modules with >5 nodes in each partition
Nmod5_gtest = zeros(size(gamma_values,2),1);
Nmod6_gtest = zeros(size(gamma_values,2),1);
Nmod7_gtest = zeros(size(gamma_values,2),1);
Nmod8_gtest = zeros(size(gamma_values,2),1);
Nmod9_gtest = zeros(size(gamma_values,2),1);
Nmod10_gtest = zeros(size(gamma_values,2),1);
Nmod20_gtest = zeros(size(gamma_values,2),1);

for k = 1 :  size(gamma_values,2)
    Nmod5_gtest(k) = sum(groupcounts(S_gtest{k})>4);
    Nmod6_gtest(k) = sum(groupcounts(S_gtest{k})>5);
    Nmod7_gtest(k) = sum(groupcounts(S_gtest{k})>6);
    Nmod8_gtest(k) = sum(groupcounts(S_gtest{k})>7);
    Nmod9_gtest(k) = sum(groupcounts(S_gtest{k})>8);
    Nmod10_gtest(k) = sum(groupcounts(S_gtest{k})>9);
    Nmod20_gtest(k) = sum(groupcounts(S_gtest{k})>19);
end

% Plot gamma versus Q, and gamma versus number of modules.
plot(gamma_values,Q_gtest)
hold on; % to retain the current plot when adding new plot
plot(gamma_values,Nmod5_gtest/max(Nmod5_gtest))
hold on; % to retain the current plot when adding new plot
plot(gamma_values,Nmod20_gtest/max(Nmod20_gtest))

% Create partition_array
partition_array = zeros(size(S_gtest,1),size(S_gtest{1},1));

for k=1 : size(S_gtest,1)
    partition_array(k,:) = S_gtest{k};
end