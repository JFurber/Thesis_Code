% Read data from CSV file
datalocation = readmatrix("F2_2015_data_location.csv");

% Adjust the 7th column to be 1-based instead of 0-based
datalocation(:,7) = datalocation(:,6)+1; %move to 1 base

Site = '2015';


% Plot the clustered data to identify each cluster to change below if needed (so numbers are next to each other)
figure();
Coords_color1 = datalocation(:,7);
k = 5; %Amount of clusters plotting
% Initialize arrays to store centroid positions
centroid_x = zeros(k, 1);
centroid_y = zeros(k, 1);
% Calculate the centroid for each cluster
for i = 1:k
    cluster_points = datalocation(Coords_color1 == i, :);
    centroid_x(i) = mean(cluster_points(:, 3));
    centroid_y(i) = mean(cluster_points(:, 4));
end
CM = summer(k);
% CM = jet(k);
for i = 1:k
    hold on
    scatter(datalocation(Coords_color1== i,3),datalocation(Coords_color1== i,4),15,CM(i,:))
    % scatter(centroid_x(i),centroid_y(i))
    text(centroid_x(i),centroid_y(i),num2str(i));
end


%%
% Woodchester

% Map the cluster identifiers to new values
for i = 1:length(datalocation)
    if datalocation(i,7)==1
        datalocation(i,8)=2;
    end
    if datalocation(i,7)==2
        datalocation(i,8)=6;
    end
    if datalocation(i,7)==3
        datalocation(i,8)=4;
    end
    if datalocation(i,7)==4
        datalocation(i,8)=7;
    end
    if datalocation(i,7)==5
        datalocation(i,8)=3;
    end
    if datalocation(i,7)==6
        datalocation(i,8)=1;
    end
    if datalocation(i,7)==7
        datalocation(i,8)=5;
    end
    % if datalocation(i,7)==8
    %     datalocation(i,8)=10;
    % end
    % if datalocation(i,7)==9
    %     datalocation(i,8)=6;
    % end
    % if datalocation(i,7)==10
    %     datalocation(i,8)=9;
    % end
    % if datalocation(i,7)==11
    %     datalocation(i,8)=4;
    % end
    % if datalocation(i,7)==12
    %     datalocation(i,8)=8;
    % end
  
end
datalocation(:,9) = datalocation(:,8);
%%
% Ireland
% Map the cluster identifiers to new values
for i = 1:length(datalocation)
    if datalocation(i,7)==1
        datalocation(i,8)=12;
    end
    if datalocation(i,7)==2
        datalocation(i,8)=5;
    end
    if datalocation(i,7)==3
        datalocation(i,8)=17;
    end
    if datalocation(i,7)==4
        datalocation(i,8)=22;
    end
    if datalocation(i,7)==5
        datalocation(i,8)=21;
    end
    if datalocation(i,7)==6
        datalocation(i,8)=16;
    end
    if datalocation(i,7)==7
        datalocation(i,8)=7;
    end
    if datalocation(i,7)==8
        datalocation(i,8)=1;
    end
    if datalocation(i,7)==9
        datalocation(i,8)=9;
    end
    if datalocation(i,7)==10
        datalocation(i,8)=25;
    end
    if datalocation(i,7)==11
        datalocation(i,8)=24;
    end
    if datalocation(i,7)==12
        datalocation(i,8)=3;
    end
    if datalocation(i,7)==13
        datalocation(i,8)=20;
    end
    if datalocation(i,7)==14
        datalocation(i,8)=11;
    end
    if datalocation(i,7)==15
        datalocation(i,8)=6;
    end
    if datalocation(i,7)==16
        datalocation(i,8)=4;
    end
    if datalocation(i,7)==17
        datalocation(i,8)=10;
    end
    if datalocation(i,7)==18
        datalocation(i,8)=13;
    end
    if datalocation(i,7)==19
        datalocation(i,8)=8;
    end
    if datalocation(i,7)==20
        datalocation(i,8)=2;
    end
    if datalocation(i,7)==21
        datalocation(i,8)=15;
    end
    if datalocation(i,7)==22
        datalocation(i,8)=23;
    end
    if datalocation(i,7)==23
        datalocation(i,8)=19;
    end
    if datalocation(i,7)==24
        datalocation(i,8)=14;
    end
    if datalocation(i,7)==25
        datalocation(i,8)=18;
    end
    % if datalocation(i,7)==26
    %     datalocation(i,8)=6;
    % end
    % if datalocation(i,7)==27
    %     datalocation(i,8)=21;
    % end
    % if datalocation(i,7)==28
    %     datalocation(i,8)=24;
    % end
end
datalocation(:,9) = datalocation(:,8);
%% C2
% Map the cluster identifiers to new values
for i = 1:length(datalocation)
    if datalocation(i,7)==1
        datalocation(i,8)=4;
    end
    if datalocation(i,7)==2
        datalocation(i,8)=1;
    end
    if datalocation(i,7)==3
        datalocation(i,8)=2;
    end
    if datalocation(i,7)==4
        datalocation(i,8)=3;
    end
    % if datalocation(i,7)==5
    %     datalocation(i,8)=4;
    % end
    % if datalocation(i,7)==6
    %     datalocation(i,8)=2;
    % end
    % if datalocation(i,7)==7
    %     datalocation(i,8)=7;
    % end
end
datalocation(:,9) = datalocation(:,8);
%% C4
% Map the cluster identifiers to new values
for i = 1:length(datalocation)
    if datalocation(i,7)==1
        datalocation(i,8)=1;
    end
    if datalocation(i,7)==2
        datalocation(i,8)=4;
    end
    if datalocation(i,7)==3
        datalocation(i,8)=3;
    end
    if datalocation(i,7)==4
        datalocation(i,8)=2;
    end
    % if datalocation(i,7)==5
    %     datalocation(i,8)=2;
    % end
    % if datalocation(i,7)==6
    %     datalocation(i,8)=6;
    % end
end
datalocation(:,9) = datalocation(:,8);
%% F1
% Map the cluster identifiers to new values
for i = 1:length(datalocation)

    if datalocation(i,7)==1
        datalocation(i,8)=1;
    end
    if datalocation(i,7)==2
        datalocation(i,8)=3;
    end
    if datalocation(i,7)==3
        datalocation(i,8)=2;
    end
    if datalocation(i,7)==4
        datalocation(i,8)=4;
    end
    % if datalocation(i,7)==5
    %     datalocation(i,8)=5;
    % end
end
datalocation(:,9) = datalocation(:,8);
%% F2
% Map the cluster identifiers to new values
for i = 1:length(datalocation)
    if datalocation(i,7)==1
        datalocation(i,8)=1;
    end
    if datalocation(i,7)==2
        datalocation(i,8)=4;
    end
    if datalocation(i,7)==3
        datalocation(i,8)=5;
    end
    if datalocation(i,7)==4
        datalocation(i,8)=3;
    end
    if datalocation(i,7)==5
        datalocation(i,8)=2;
    end
    % if datalocation(i,7)==6
    %     datalocation(i,8)=8;
    % end
    % if datalocation(i,7)==7
    %     datalocation(i,8)=4;
    % end
    % if datalocation(i,7)==8
    %     datalocation(i,8)=3;
    % end
    % if datalocation(i,7)==9
    %     datalocation(i,8)=5;
    % end
end
datalocation(:,9) = datalocation(:,8);

%%

datalocation(:,8) = datalocation(:,7);
datalocation(:,9) = datalocation(:,8);


%%

% Split data into groups based on unique badger IDs
[~,~,gId] = unique(datalocation(:,1));
groups = splitapply(@(x){x}, datalocation, gId);
groups = transpose(groups);

%Calculate the length of each group
len = zeros(1,length(groups));
for i = 1:length(groups)
    len(i) = length(groups{i});
end

% %Pad the groups to have the same length
% for j = 1:length(groups)
%     if length(groups{j}) ~= max(len)
%     Badger_zero{j} = padarray(groups{j},max(len)-length(groups{j}),0,'post');
%     else
%         Badger_zero{j} = groups{j};
%     end
% endc




%
% For Markov Chain

% Define the column index for specific layout and interpolation time
a = 9; %7 is not changed excel to specific layout %9 if have

interp = 20; %time that you have interpolated for

% Pad the groups to have the same length and add interpolation time
for j = 1:length(groups)
    if length(groups{j}) ~= max(len)
    Badger_zero1{j} = padarray(groups{j},max(len)-length(groups{j}),0,'post'); %pad the array with 0s
    Badger_zero1{j}(:,1) = j; %give first row a number
    for k = len(j)+1:max(len)
        Badger_zero1{j}(k,2) = Badger_zero1{j}(k-1,2)+interp; %+interpolation time
        Badger_zero1{j}(k,a) = 30; %'removed badgers'
    end
    else
        Badger_zero1{j} = groups{j};
    end
end

% Concatenate all padded groups into one array
Badger_M = vertcat(Badger_zero1{:}); %put into one array 

% Split into same time (unique by time)
[~,~,gId] = unique(Badger_M(:,2));
groups_TIME = splitapply(@(x){x}, Badger_M, gId);
groups_TIME = transpose(groups_TIME);

%

Badgers_Markov = groups_TIME{1}(:,a);
for p = 2:length(groups_TIME)
    Badgers_Markov = [Badgers_Markov, groups_TIME{p}(:,a)];
end

clear Badger_zero Badger_zero1 p len k i j groups groups_TIME gId

%% Find how many are in each box at each time point - not necessary

Badgers_Markov_SUM = zeros(9,length(Badgers_Markov));
for i=1:length(Badgers_Markov_SUM)
    Badgers_Markov_SUM(1,i) = sum(Badgers_Markov(:,i)==0);
    Badgers_Markov_SUM(2,i) = sum(Badgers_Markov(:,i)==1);
    Badgers_Markov_SUM(3,i) = sum(Badgers_Markov(:,i)==2);
    Badgers_Markov_SUM(4,i) = sum(Badgers_Markov(:,i)==3);
    Badgers_Markov_SUM(5,i) = sum(Badgers_Markov(:,i)==4);
    Badgers_Markov_SUM(6,i) = sum(Badgers_Markov(:,i)==5);
    Badgers_Markov_SUM(7,i) = sum(Badgers_Markov(:,i)==6);
    Badgers_Markov_SUM(8,i) = sum(Badgers_Markov(:,i)==7);
    Badgers_Markov_SUM(9,i) = sum(Badgers_Markov(:,i)==8);
end

%%
% Define the number of clusters
clusters = 5;

% Initialize the transition matrix
matrix = zeros(clusters);

% Loop through each row of the Badgers_Markov matrix
for k = 1:size(Badgers_Markov,1)
    % Initialize a temporary matrix for the current row
    matrix1 = zeros(clusters);
    % Get unique states and their indices
    [u,~,n] = unique(Badgers_Markov(k,:));
    % NumU is the number of unique states
    NumU = length(u);     
    % Calculate the counts of transitions between states
    Counts = accumarray([n(1:end-1),n(2:end)],1,[NumU,NumU]);
    % If the last state is 'removed', ignore it in the transition matrix
    if u(:,length(u))==30 
        for i = 1:length(u)-1
            for j = 1:length(u)-1
             matrix1(u(i),u(j)) = Counts(i,j);
            end
        end
    else 
        for i = 1:length(u)
            for j = 1:length(u)
             matrix1(u(i),u(j)) = Counts(i,j);
            end
        end
    end
    % Add the temporary matrix to the main transition matrix
    matrix = matrix + matrix1;
end

% Normalize the transition matrix to get the probability transition matrix
P1 = matrix./sum(matrix,2);  

P2 = transpose(P1); % sum of columns is 1

clear u n NumU matrix1 k j i

%% Plot graph

% Create a discrete-time Markov chain object
mc = dtmc(P2);

% Get the number of states in the Markov chain
numstates = mc.NumStates

figure()

% Plot the Markov chain graph with colored edges
h = graphplot(mc, 'ColorEdges',true)

% coordinates of the nodes to space them out

% C2
% h.XData = [3 2 1 1.5 2.2];
% h.YData = [1 1 1 1.5 2];

% % C4
% h.XData = [3.7 3 2.5 2 1.5 1];
% h.YData = [1.5 1.5 2 3 2.5 3];

% F1
% h.XData = [1 3.25 2 3 4];
% h.YData = [1 1 1.5 2 2.5];

% F2
% h.XData = [0.5 1.5 2.5 2.5 2.5 2.5 1 0.5];
% h.YData = [4 4 4 3 2 1 2 2.5];

% Woodchester
% h.XData = [1 2 2 3 4 3 4 5];
% h.YData = [1.5 2 1 2 2 1 1 1.5];

% % Ireland
% h.XData = [-0.5 -0.5 -0.5 -0.5 1 2 3 3 3 4 3 3 1 5 5 8 7 6.5 9 8 9 8 7 7 6 3 4 5];
% h.YData = [4 5.5 7 8.5 7 4 6 7 8 8 10 9 10 9 10.5 10.5 9 7.5 7.5 6 4 3 4 2 0.5 2 3 6];


% title('Transition Graph (Ireland)',Interpreter='latex')
% latex_fig(15, 3.85, 5)

%% 

numSteps = 50;
X = simulate(mc,numSteps);
figure;
simplot(mc,X);
%% 

x0 = 1*ones(1,mc.NumStates);
numSteps = 10;
X = simulate(mc,numSteps,'X0',x0);

%% Plot clusters

mc2 = dtmc(P2);
Coords_color1 = datalocation(:,9);
CM = summer(clusters);
% CM = jet(k);
% Number of clusters
k = max(Coords_color1);

% Initialize arrays to store centroid positions
centroid_x = zeros(k, 1);
centroid_y = zeros(k, 1);

% Calculate the centroid for each cluster
for i = 1:k
    cluster_points = datalocation(Coords_color1 == i, :);
    centroid_x(i) = mean(cluster_points(:, 3));
    centroid_y(i) = mean(cluster_points(:, 4));
end

figure();
hold on;

for i = 1:clusters
    hold on
    scatter(datalocation(Coords_color1== i,3),datalocation(Coords_color1== i,4),15,CM(i,:))
end

h = graphplot(mc2, 'ColorEdges',true)

h.XData = centroid_x ;
h.YData = centroid_y;
h.MarkerSize = 3;
h.NodeColor = [0 0 0];

% h.NodeLabel = [];
colorbar.label = 'Transition Probability';
colormap cool
title(Site, Interpreter='latex');
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('Easting', Interpreter='latex')
ylabel('Northing', Interpreter='latex')
set(gca,'TickLabelInterpreter','latex')
% colorbar(TickLabelInterpreter='latex')
hold off
latex_fig(15, 2.6, 2.4)

% Cluster Assignments and Transition Probabilities

%%

print -depsc F215_MC_26_24.eps
% print -dpng spec_compare.png
% savefig('transition.fig')
savefig('F2_MC_2015.fig')

%% Permeable Barrier

transitionMatrix = P2; % Your kxk transition matrix
[permMatrix, countOffDiagonalNonZero] = Permability(transitionMatrix);
% disp(permMatrix);
% disp(['Number of off-diagonal non-zero entries: ', num2str(countOffDiagonalNonZero)]);

% Site = 'Cornwall (F2)';

figure();
% Permeable = P2;
heatmap(permMatrix,'CellLabelColor','none')

% Define the number of points for negative and positive values
total_points = countOffDiagonalNonZero;

%Define min and max values
min_value = min(permMatrix(:));
max_value = max(permMatrix(:));

% Calculate the proportion of negative and positive ranges
negative_range = -min_value;
positive_range = max_value;
total_range = negative_range + positive_range;

% Calculate the number of points for negative and positive ranges
negative_points = round(total_points * (negative_range / total_range));
positive_points = total_points - negative_points;
% 
% % Generate the colormap for negative values (red to white)
% red_negative = [linspace(1, 1, negative_points)', linspace(0, 1, negative_points)', linspace(0, 1, negative_points)'];
% 
% % Generate the colormap for positive values (white to green)
% green_positive = [linspace(1, 0, positive_points)', linspace(1, 1, positive_points)', linspace(1, 0, positive_points)'];
% 
% % Combine the negative, zero, and positive colormaps
% custom_colormap = [red_negative; [1, 1, 1]; green_positive];

% Create a transition from blue to white
blue_to_white = [linspace(0, 1, negative_points)', linspace(0, 1, negative_points)', linspace(1, 1, negative_points)'];

% Create a transition from white to yellow
white_to_pink = [linspace(1, 1, positive_points)', linspace(1, 0.08, positive_points)', linspace(1, 0.58, positive_points)'];

custom_colormap = [blue_to_white;[1,1,1]; white_to_pink];

% Apply the custom colormap
colormap(custom_colormap);

title(Site);


latex_fig(15, 2.6, 2.4)

%%

print -depsc F215_perm_26_24.eps
savefig('F215_permeable.fig')


%% Put MC into 1 Figure

%  Load saved figures
aa = openfig('WCT_MC.fig', 'reuse', 'invisible');
a = openfig('WC18_MC.fig', 'reuse', 'invisible');
b = openfig('WC19_MC.fig', 'reuse', 'invisible');
c = openfig('WC20_MC.fig', 'reuse', 'invisible');
d = openfig('WC21_MC.fig', 'reuse', 'invisible');
e = openfig('WC22_MC.fig', 'reuse', 'invisible');

% Prepare subplots
figure;

h(1) = subplot('Position',[0.097,0.5503,0.2134,0.3188])
ylabel("Northing",'Interpreter','latex')

colormap cool
set(gca,'XTick',[])
set(gca,'YTick',[])

h(2) = subplot('Position',[0.3778,0.5503,0.2134,0.3188])
colormap cool
set(gca,'XTick',[])
set(gca,'YTick',[])

h(3) = subplot('Position',[0.6586,0.5503,0.2134,0.3188])
colorbar.label = 'Transition Probability';
colormap cool
set(gca,'XTick',[])
set(gca,'YTick',[])

h(4) = subplot('Position',[0.097,0.1037, 0.2134,0.3188])
ylabel("Northing",'Interpreter','latex')
colormap cool
set(gca,'XTick',[])
set(gca,'YTick',[])

h(5) = subplot('Position',[0.3778,0.1037,0.2134,0.3188])
xlabel("Easting",'Interpreter','latex')
colormap cool
set(gca,'XTick',[])
set(gca,'YTick',[])

h(6) = subplot('Position',[0.6586,0.1037,0.2134,0.3188])
colorbar.label = 'Transition Probability';
colormap cool
set(gca,'XTick',[])
set(gca,'YTick',[])

% Copy figures to the subplots
copyobj(allchild(get(aa, 'CurrentAxes')), h(1));
copyobj(allchild(get(a, 'CurrentAxes')), h(2));
copyobj(allchild(get(b, 'CurrentAxes')), h(3));
copyobj(allchild(get(c, 'CurrentAxes')), h(4));
copyobj(allchild(get(d, 'CurrentAxes')), h(5));
copyobj(allchild(get(e, 'CurrentAxes')), h(6));

% Set titles for subplots if needed
title(h(1), 'Total','Interpreter','latex');
title(h(2), '2018','Interpreter','latex');
title(h(3), '2019','Interpreter','latex');
title(h(4), '2020','Interpreter','latex');
title(h(5), '2021','Interpreter','latex');
title(h(6), '2022','Interpreter','latex');

% Adjust other properties if needed
set(h(1), 'TickLabelInterpreter', 'latex');
set(h(2), 'TickLabelInterpreter', 'latex');
set(h(3), 'TickLabelInterpreter', 'latex');
set(h(4), 'TickLabelInterpreter', 'latex');
set(h(5), 'TickLabelInterpreter', 'latex');
set(h(6), 'TickLabelInterpreter', 'latex');

sgtitle('Woodchester','Interpreter','latex')


latex_fig(15, 6, 3.5)

% print -depsc WC_MC_6_35.eps
savefig('WC_MC.fig')

% close all hidden

%%

figFiles = {'WCT_permeable.fig', 'WC18_permeable.fig', 'WC19_permeable.fig'};

for i = 1:length(figFiles)
    % Open the figure file
    aa = openfig(figFiles{i}, 'reuse', 'invisible');
    
    % Check if the figure was opened successfully
    if isempty(aa)
        error('Failed to open figure file: %s', figFiles{i});
    end
    
    % Find the HeatmapChart object using findall
    heatmapObj = findall(aa, 'Type', 'HeatmapChart');
    
    % Check if the HeatmapChart object was found
    if isempty(heatmapObj)
        error('No HeatmapChart object found in figure file: %s', figFiles{i});
    end
    
    % Create a subplot
    h(i) = subplot(3, 2, i);
    
    % Copy the HeatmapChart object to the subplot
    newHeatmap = copyobj(heatmapObj, gcf);
    
    % Adjust the position of the new heatmap to fit the subplot
    newHeatmap.Position = h(i).Position;
    
    % Delete the original subplot axes
    delete(h(i));
    
    % Close the original figure
    close(aa);
end

% Refresh the figure window to ensure all graphics are updated
drawnow;
