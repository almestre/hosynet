% Modularity of the original bipartite network
%    Lay out the B matrix which will be fed to the genlouvain algorithm to calculate Q
%    The input (A) is an MxN adjacency matrix for an undirected bipartite network
%    The function returns monolayer Barber modularity matrix for undirected bipartite networks

% Set the random seed for repeatability 
seed = 547857;
rng(seed);

gamma = 1.3;
[B,twom,A_mat,P_mat]=bipartite2(A,gamma);

% Create the containers to store the data generated in the for loop
Q_rep = zeros(100,1);
S_rep = cell(100,1);
n_it_rep =  zeros(100,1);

for k = 1 : 100
    % Genlouvain algorithm
    [S_rep{k},Q_rep(k),n_it_rep(k)]=iterated_genlouvain(B);
    Q_rep(k)=Q_rep(k)/twom;
end

% Number of modules with >4 nodes in each partition
Nmod20_rep = zeros(100,1);

% Create an array to store partitions
part_g1 = zeros(size(S_rep,1),size(S_rep{1},1));

for k=1 : size(S_rep,1)
    part_g1(k,:) = S_rep{k};
end

for k = 1 :  100
    Nmod20_rep(k) = sum(groupcounts(S_rep{k})>19);
end

% Identify the replicate with the max Q
[Q_best,idx_best] = max(Q_rep);
S_best = S_rep{idx_best};
Nmod20_best = Nmod20_rep(idx_best);

% Modularity of the 999 null models
B_null = cell(size(A_null,1),1);
twom_null = zeros(size(A_null,1),1);
S_null = cell(size(A_null,1),1);
Q_null = zeros(size(A_null,1),1);
n_it_null = zeros(size(A_null,1),1);

for k1 = 1 : size(A_null,1)*100
    if k1 == 1
        k2 = 1;
    end
    
    [B_null_temp,twom_null_temp] = bipartite(A_null{k2},gamma);
    [S_null_temp,Q_null_temp,n_it_null_temp]=iterated_genlouvain(B_null_temp);
    Q_null_temp=Q_null_temp/twom_null_temp;
    
    if Q_null(k2) < Q_null_temp
        B_null{k2} = B_null_temp;
        twom_null(k2) = twom_null_temp;
        Q_null(k2) = Q_null_temp;
        S_null{k2} = S_null_temp;
        n_it_null(k2) =  n_it_null_temp;
    end
    
    if mod(k1,100) == 0
        k2 = k2 + 1;
    end
end

% Calculating the z-score (standard score or normal score) of Q:
zQ = (Q_best - mean(Q_null))/std(Q_null);

% Historgam with a normal distribution fit:
h_modtest = histfit(Q_null);

% plottools('on')

% Set papersize
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [16 16]);

% Set figure size
set(gcf, 'Windowstyle', 'normal');
x0=1;
y0=1;
width=14;
height=12;
set(gcf,'units','centimeters','position',[x0,y0,width,height]);

% Set tick values for the axes
xticks([.48 .5 .52 .54, .56, .58, .6]);
xticklabels({'0.48','0.50','0.52','0.54','0.56','0.58','0.60'});
yticks([0 20 40 60 80 100]);
yticklabels({'0','20','40','60','80','100'});

% Set axes labels
xlabel('Q');
ylabel('Density');

% Set colors
h_modtest(1).FaceColor = [173/255 216/255 230/255]; % color of hist
h_modtest(2).Color = [0 0 0]; % color of fitted line

% Set figure axes
ax = gca;
ax.FontSize = 8;
ax.XLim = [.475 .6];
ax.YLim = [0 110];

% Add line with the observed value of Q
hold on; % to retain the current plot when adding new plot
line([Q_best, Q_best], ylim, 'LineWidth', 2, 'Color', 'r');

% Save figure as jpg
%saveas(gcf,'h_modtest.jpg'); % adjusts image size to figure size

% Print figure to pdf
print(gcf,'-dpdf','htest_g1.pdf'); % Papersize independent of figure size and defined in the script (see "Set papersize")