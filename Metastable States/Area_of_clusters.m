%% Import Data
clear
clc

Site = '2017';

% Read the CSV file for data location
data = readmatrix('Ireland_2017_data_location.csv');

% Extract the coordinates and labels
coordinates = data(:, 3:4);
labels = data(:, 6);

% Find unique labels
unique_labels = unique(labels);

cluster_areas_km2 = zeros(size(unique_labels));
for i = 1:numel(unique_labels)
    cluster_coordinates = coordinates(labels == unique_labels(i), :);
    k_clust = convhull(cluster_coordinates(:, 1), cluster_coordinates(:, 2));
    cluster_areas_km2(i) = polyarea(cluster_coordinates(k_clust, 1), cluster_coordinates(k_clust, 2))/(1000*1000);
end

clear data

% Read CSV file for SG information
data = readmatrix('Ireland2017_AUGUST2024.csv');

% Cornwall
% coordinates_data = data(:, 17:18);
% labels_SG = data(:, 6); %C2
% labels_SG = data(:, 7); %C4, F1,F2

% Ireland
coordinates_data = data(:, 21:22);
labels_SG = data(:, 7);

%Yearly WC
% coordinates_data = data(:, 21:22);
% labels_SG = data(:, 5);

%For Total WC:
% coordinates_data = data(:, 25:26);
% labels_SG = data(:, 6);

% Find unique labels
unique_labels_SG = unique(labels_SG);

% Calculate the area of each cluster
cluster_areas_km2_SG = zeros(size(unique_labels_SG));
for i = 1:numel(unique_labels_SG)
    cluster_coordinates_SG = coordinates_data(labels_SG == unique_labels_SG(i), :);
    k = convhull(cluster_coordinates_SG(:, 1), cluster_coordinates_SG(:, 2));
    cluster_areas_km2_SG(i) = polyarea(cluster_coordinates_SG(k, 1), cluster_coordinates_SG(k, 2))/(1000*1000);
end


%%

% print -depsc EDMD_fig.eps
% print -dpng spec_compare.png
savefig('17NIEDMD_area1.fig')



%% Metastable cluster and Social Group 


figure;
hold on
%Plot metastable clusters

colors = parula(numel(unique_labels));

% Generate a custom grayscale colormap excluding white
num_clusters = numel(unique_labels);
% gray_colors = linspace(0, 0.8, num_clusters)' * [1 1 1]; % Shades from black to light gray

for i = 1:num_clusters
    cluster_coordinates = coordinates(labels == unique_labels(i), :);
    
    % Plot the points (optional, can be removed if only convex hulls are needed)
    plot(cluster_coordinates(:, 1), cluster_coordinates(:, 2), 'o', 'Color', colors(i, :));
    
    % Calculate and plot the convex hull
    k_clust = convhull(cluster_coordinates(:, 1), cluster_coordinates(:, 2));
    a = plot(cluster_coordinates(k_clust, 1), cluster_coordinates(k_clust, 2), 'Color', '#5D3A9B','LineWidth',3);

    cluster_areas = polyarea(cluster_coordinates(k_clust, 1), cluster_coordinates(k_clust, 2))/(1000*1000);
    
    % Plot the filled the convex hull
    % fill(cluster_coordinates(k_clust, 1), cluster_coordinates(k_clust, 2), colors(i, :), 'FaceAlpha', 0.5);
    
    % Label the cluster with its area
    % text(mean(cluster_coordinates(:, 1)), mean(cluster_coordinates(:, 2)), sprintf('%.2f km$^2$', cluster_areas), 'Color', 'k',Interpreter='latex');
end

hold off;

hold on;

% colors_SG = jet(numel(unique_labels_SG)); % Generate distinct colors for each cluster

for i = 1:numel(unique_labels_SG)
    cluster_coordinates_SG = coordinates_data(labels_SG == unique_labels_SG(i), :);
    % plot(cluster_coordinates_SG(:, 1), cluster_coordinates_SG(:, 2), 'o','Color', colors_SG(i, :));
    
    % Calculate and plot the convex hull
    k = convhull(cluster_coordinates_SG(:, 1), cluster_coordinates_SG(:, 2));
    % k = boundary(cluster_coordinates_SG(:, 1), cluster_coordinates_SG(:, 2));
    % plot(cluster_coordinates_SG(k, 1), cluster_coordinates_SG(k, 2), 'r');
    

    % Plot and fill the convex hull
    b = plot(cluster_coordinates_SG(k, 1), cluster_coordinates_SG(k, 2), 'Color', '#E66100','LineWidth',3);
    
    % fill(cluster_coordinates_SG(k, 1), cluster_coordinates_SG(k, 2), colors(i, :), 'FaceAlpha', 0.5);

    % Label the cluster with its area
    % cluster_area_km2_SG = cluster_areas_SG(i) / (1000 * 1000);
    % text(mean(cluster_coordinates_SG(:, 1)), mean(cluster_coordinates_SG(:, 2)), sprintf('%.2f km^2', cluster_area_km2_SG));
    
    % Add legend for clusters
    %legend('Social Groups');
end
hold off;

xlabel('x (meters)', Interpreter='latex')
ylabel('y (meters)', Interpreter='latex')
% set(gca,'YTick',[])
set(gca,'TickLabelInterpreter','latex')

title(Site,Interpreter="latex")

legend([a,b],'Metastable Cluster','Social Group','Interpreter','latex');

latex_fig(15, 6, 5)


%% Plot Metastable Clusters (for Total Data)

% Plot the polyareas and convex hulls

figure;

hold on
%Plot metastable clusters

colors = parula(numel(unique_labels));

% Generate a custom grayscale colormap excluding white
num_clusters = numel(unique_labels);
% gray_colors = linspace(0, 0.8, num_clusters)' * [1 1 1]; % Shades from black to light gray

for i = 1:num_clusters
    cluster_coordinates = coordinates(labels == unique_labels(i), :);
    
    % Plot the points (optional, can be removed if only convex hulls are needed)
    a = plot(cluster_coordinates(:, 1), cluster_coordinates(:, 2), 'o', 'Color', colors(i, :));
    
    % Calculate and plot the convex hull
    k_clust = convhull(cluster_coordinates(:, 1), cluster_coordinates(:, 2));
    a = plot(cluster_coordinates(k_clust, 1), cluster_coordinates(k_clust, 2), 'Color', '#5D3A9B','LineWidth',3);

    cluster_areas = polyarea(cluster_coordinates(k_clust, 1), cluster_coordinates(k_clust, 2))/(1000*1000);
    
    % Plot the filled the convex hull
    % fill(cluster_coordinates(k_clust, 1), cluster_coordinates(k_clust, 2), colors(i, :), 'FaceAlpha', 0.5);
    
    % Label the cluster with its area
    % text(mean(cluster_coordinates(:, 1)), mean(cluster_coordinates(:, 2)), sprintf('%.2f km$^2$', cluster_areas), 'Color', 'k','Interpreter','latex');
end

hold off;

xlabel('x (meters)', Interpreter='latex')
ylabel('y (meters)', Interpreter='latex')
% set(gca,'YTick',[])
set(gca,'TickLabelInterpreter','latex')

title(Site,Interpreter="latex")

latex_fig(15, 6, 5)

%%

% print -depsc EDMD_fig.eps
% print -dpng spec_compare.png
savefig('14NI_EDMD_meta.fig')

%% Plot Metastable Clusters (for 1 Metastable Cluster)

% Plot the polyareas and convex hulls

figure;

hold on
%Plot metastable clusters

colors = parula(numel(unique_labels_SG));

% Generate a custom grayscale colormap excluding white
num_clusters = numel(unique_labels_SG);
% gray_colors = linspace(0, 0.8, num_clusters)' * [1 1 1]; % Shades from black to light gray

for i = 1:num_clusters
    cluster_coordinates_SG = coordinates_data(labels_SG == unique_labels_SG(i), :);
    
    % Plot the points (optional, can be removed if only convex hulls are needed)
    a = plot(cluster_coordinates_SG(:, 1), cluster_coordinates_SG(:, 2),'o', 'Color', colors(i, :));
    
    % Calculate and plot the convex hull
    k_clust = convhull(cluster_coordinates_SG(:, 1), cluster_coordinates_SG(:, 2));
    a = plot(cluster_coordinates_SG(k_clust, 1), cluster_coordinates_SG(k_clust, 2), 'Color', '#5D3A9B','LineWidth',3);

    cluster_areas = polyarea(cluster_coordinates_SG(k_clust, 1), cluster_coordinates_SG(k_clust, 2))/(1000*1000);
    
    % Plot the filled the convex hull
    % fill(cluster_coordinates(k_clust, 1), cluster_coordinates(k_clust, 2), colors(i, :), 'FaceAlpha', 0.5);
    
    % Label the cluster with its area
    text(mean(cluster_coordinates_SG(:, 1)), mean(cluster_coordinates_SG(:, 2)), sprintf('%.2f km$^2$', cluster_areas), 'Color', 'k','Interpreter','latex');
end

hold off;

xlabel('x (meters)', Interpreter='latex')
ylabel('y (meters)', Interpreter='latex')
% set(gca,'YTick',[])
set(gca,'TickLabelInterpreter','latex')

title(Site,Interpreter="latex")

latex_fig(15, 6, 5)

%%

figure;
hold on;
for i = 1:numel(unique_labels)
    cluster_coordinates = coordinates(labels == unique_labels(i), :);
    plot(cluster_coordinates(:, 1), cluster_coordinates(:, 2), 'o');

    % Calculate and plot the convex hull
    k = convhull(cluster_coordinates(:, 1), cluster_coordinates(:, 2));
    plot(cluster_coordinates(k, 1), cluster_coordinates(k, 2), 'r');

    % Label the cluster with its area
    cluster_area_km2 = cluster_areas(i) / (1000 * 1000);
    text(mean(cluster_coordinates(:, 1)), mean(cluster_coordinates(:, 2)), sprintf('%.2f km^2', cluster_area_km2));

    % Add legend for clusters
    %legend('Clusters');
end
hold off;
% 
cluster_areas_km2 = cluster_areas/(1000*1000) ;

% Display the cluster areas in command window
% disp('Cluster Areas (km^2):');
% disp(cluster_areas_km2);
% 
% clear data

%% Plot Area of Social Groups

% Plot the polyareas and convex hulls
figure;
hold on;
for i = 1:numel(unique_labels_SG)
    cluster_coordinates_SG = coordinates_data(labels_SG == unique_labels_SG(i), :);
    plot(cluster_coordinates_SG(:, 1), cluster_coordinates_SG(:, 2), 'o');
    
    % Calculate and plot the convex hull
    k = convhull(cluster_coordinates_SG(:, 1), cluster_coordinates_SG(:, 2));
    plot(cluster_coordinates_SG(k, 1), cluster_coordinates_SG(k, 2), 'r');
    
    % Label the cluster with its area
    cluster_area_km2_SG = cluster_areas_SG(i) / (1000 * 1000);
    text(mean(cluster_coordinates_SG(:, 1)), mean(cluster_coordinates_SG(:, 2)), sprintf('%.2f km^2', cluster_area_km2_SG));
    
    % Add legend for clusters
    %legend('Social Groups');
end
hold off;

cluster_areas_km2_SG = cluster_areas_SG/(1000*1000) ;

% Plot both together

% Plot the polyareas and convex hulls in grayscale without white

figure;


hold on;

colors = jet(numel(unique_labels));

% Generate a custom grayscale colormap excluding white
num_clusters = numel(unique_labels);
% gray_colors = linspace(0, 0.8, num_clusters)' * [1 1 1]; % Shades from black to light gray

for i = 1:num_clusters
    cluster_coordinates = coordinates(labels == unique_labels(i), :);
    
    % Plot the points (optional, can be removed if only convex hulls are needed)
    % plot(cluster_coordinates(:, 1), cluster_coordinates(:, 2), 'o', 'Color', gray_colors(i, :));
    
    % Calculate and plot the convex hull
    k_clust = convhull(cluster_coordinates(:, 1), cluster_coordinates(:, 2));
    % k_clust = boundary(cluster_coordinates(:, 1), cluster_coordinates(:, 2));
    a = plot(cluster_coordinates(k_clust, 1), cluster_coordinates(k_clust, 2), 'Color', 'b','LineWidth',3);

    % Plot and fill the convex hull
    fill(cluster_coordinates(k_clust, 1), cluster_coordinates(k_clust, 2), colors(i, :), 'FaceAlpha', 0.5);
    
    % Label the cluster with its area
    cluster_area_km2 = cluster_areas(i) / (1000 * 1000);
    text(mean(cluster_coordinates(:, 1)), mean(cluster_coordinates(:, 2)), sprintf('%.2f km^2', cluster_area_km2), 'Color', 'k');
end


hold off;

hold on;

colors_SG = jet(numel(unique_labels_SG)); % Generate distinct colors for each cluster

for i = 1:numel(unique_labels_SG)
    cluster_coordinates_SG = coordinates_data(labels_SG == unique_labels_SG(i), :);
    % plot(cluster_coordinates_SG(:, 1), cluster_coordinates_SG(:, 2), 'o','Color', colors_SG(i, :));
    
    % Calculate and plot the convex hull
    k = convhull(cluster_coordinates_SG(:, 1), cluster_coordinates_SG(:, 2));
    % k = boundary(cluster_coordinates_SG(:, 1), cluster_coordinates_SG(:, 2));
    % plot(cluster_coordinates_SG(k, 1), cluster_coordinates_SG(k, 2), 'r');
    

    % Plot and fill the convex hull
    b = plot(cluster_coordinates_SG(k, 1), cluster_coordinates_SG(k, 2), 'Color', 'r','LineWidth',3);
    
    % fill(cluster_coordinates_SG(k, 1), cluster_coordinates_SG(k, 2), colors(i, :), 'FaceAlpha', 0.5);

    % Label the cluster with its area
    cluster_area_km2_SG = cluster_areas_SG(i) / (1000 * 1000);
    % text(mean(cluster_coordinates_SG(:, 1)), mean(cluster_coordinates_SG(:, 2)), sprintf('%.2f km^2', cluster_area_km2_SG));
    
    % Add legend for clusters
    %legend('Social Groups');
end
hold off;

title(Site)

legend([a,b], 'Metastable Cluster','Social Group');



%%


%  Load saved figures
aa = openfig('TNIEDMD_area.fig', 'reuse', 'invisible');
a = openfig('14NIEDMD_area1.fig', 'reuse', 'invisible');
b = openfig('15NIEDMD_area1.fig', 'reuse', 'invisible');
c = openfig('16NIEDMD_area1.fig', 'reuse', 'invisible');
d = openfig('17NIEDMD_area1.fig', 'reuse', 'invisible');
% e = openfig('17F1EDMD_area.fig', 'reuse', 'invisible');

% aa = openfig('TNI_EDMD_meta.fig', 'reuse', 'invisible');
% a = openfig('14NI_EDMD_meta.fig', 'reuse', 'invisible');
% b = openfig('15NI_EDMD_meta.fig', 'reuse', 'invisible');
% c = openfig('16NI_EDMD_meta.fig', 'reuse', 'invisible');
% d = openfig('17NI_EDMD_meta.fig', 'reuse', 'invisible');




% Prepare subplots
figure;

h(1) = subplot(3, 2, 1);
ylabel("y (meters)",'Interpreter','latex')

h(2) = subplot(3, 2, 2);
h(3) = subplot(3, 2, 3);
ylabel("y (meters)",'Interpreter','latex')
% xlabel("x (meters)",'Interpreter','latex')
h(4) = subplot(3, 2, 4);
% xlabel("x (meters)",'Interpreter','latex')
h(5) = subplot(3, 2, 5);
ylabel("y (meters)",'Interpreter','latex')
xlabel("x (meters)",'Interpreter','latex')
% h(6) = subplot(3, 2, 6);
% xlabel("x (meters)",'Interpreter','latex')

% Copy figures to the subplots
copyobj(allchild(get(aa, 'CurrentAxes')), h(1));
copyobj(allchild(get(a, 'CurrentAxes')), h(2));
copyobj(allchild(get(b, 'CurrentAxes')), h(3));
copyobj(allchild(get(c, 'CurrentAxes')), h(4));
copyobj(allchild(get(d, 'CurrentAxes')), h(5));
% copyobj(allchild(get(e, 'CurrentAxes')), h(6));

% % Copy legends to the subplots
% legends_c = findobj(c, 'Type', 'Legend');
% legends_k = findobj(k, 'Type', 'Legend');
% 
% % Create new legends in the subplots
% new_legend_c = legend(h(1), get(legends_c, 'String'));
% new_legend_k = legend(h(2), get(legends_k, 'String'));

% Set titles for subplots if needed
title(h(1), 'Total','Interpreter','latex');
title(h(2), '2014','Interpreter','latex');
title(h(3), '2015','Interpreter','latex');
title(h(4), '2016','Interpreter','latex');
title(h(5), '2017','Interpreter','latex');
% title(h(6), '2017','Interpreter','latex');

% Adjust other properties if needed
set(h(1), 'TickLabelInterpreter', 'latex');
set(h(2), 'TickLabelInterpreter', 'latex');
set(h(3), 'TickLabelInterpreter', 'latex');
set(h(4), 'TickLabelInterpreter', 'latex');
set(h(5), 'TickLabelInterpreter', 'latex');
% set(h(6), 'TickLabelInterpreter', 'latex');



sgtitle('Northern Ireland','Interpreter','latex')

latex_fig(15, 6, 7)

% close all hidden

%% Add in legends to already saved figures

% Open the saved .fig file
fig = openfig('WC_EDMD_area_6_7.fig');

% % Get the handles of the axes
% ax = gca;
% 
% % Use findobj to get handles of the lines
% lines = findobj(ax, 'Type', 'Line');

% Get the handles of all axes in the figure
allAxes = findobj(fig, 'Type', 'axes');

% Identify the third subplot
% Assuming subplots are created in a 2x2 grid, for example
thirdSubplot = allAxes(2); % This assumes the third subplot is the third in the list

% Alternatively, you can use the subplot function if you know the layout
% thirdSubplot = subplot(2, 2, 3); % For a 2x2 grid, third subplot

% Use findobj to get handles of the lines in the third subplot
% lines = findobj(thirdSubplot, 'Type', 'Line');

% Assuming you know the properties of the lines you want to include in the legend
% For example, if you know the color or other properties of the lines
Meta = findobj(thirdSubplot, 'Type', 'Line', 'Color', '#5D3A9B');
Social = findobj(thirdSubplot, 'Type', 'Line', 'Color', '#E66100');

h1 = Meta(1);
h2 = Social(1);

% Add the legend for the two lines
legend(thirdSubplot,[h1, h2], {'Metastable Cluster', 'Social Group'});
legend(thirdSubplot,'Interpreter','latex',"NumColumns",2,'Position',[0.13,0.025,0.776,0.0234])


%%

print -depsc WC_EDMD_area_6_7.eps
% print -dpng spec_compare.png
savefig('WC_EDMD_area_6_7.fig')