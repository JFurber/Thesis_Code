% Read the data from the Excel file 'Noise_GLMM_data.xlsx' into a table
data = readtable('Noise_GLMM_data.xlsx');

% Convert 'Site', 'Capture_Year', and 'Sex' columns to categorical data type if not already
data.Site = categorical(data.Site);
data.Capture_Year = categorical(data.Capture_Year);
data.Sex = categorical(data.Sex);
data.Month = categorical(data.Month);

% Reorder the categories for 'Site' to set a specific reference order
data.Site = reordercats(data.Site, {'Ireland','Woodchester','C2', 'C4', 'F1', 'F2'});

% Reorder the categories for 'Capture_Year' to set a specific reference order
data.Capture_Year = reordercats(data.Capture_Year, {'Year 1','Year 2', 'Year 3', 'Year 4', 'Year 5'});

% Reorder the categories for 'Sex' to set a specific reference order
data.Sex = reordercats(data.Sex, {'Male','Female'});

% Fit a generalized linear mixed-effects model (GLMM) to the data
% The model predicts 'L' using 'Capture_Year', 'Sex', 'Month', and 'Site' as fixed effects
% and includes a random intercept for 'Animal'
% Distribution is normal (default for continuous data)
% Binomial distribution for binary data, poisson for count 
glme = fitglme(data,'L ~ 1 + Capture_Year + Sex + Month + Site + (1|Animal)','Distribution','Normal','Link','log');


glme_fixed = fitglme(data,'L ~ 1 + Capture_Year + Sex + Month + Site','Distribution','Normal','Link','log');


%%
% Fit a reduced GLMM excluding 'Capture_Year' as a fixed effect
glme_reduced_CY = fitglme(data,'L ~ 1 + Sex + Month + Site + (1|Animal)','Distribution','Normal','Link','log');
% Compare the reduced model with the full model to test the significance of 'Capture_Year'
result_CY = compare(glme_reduced_CY, glme)

% Fit a reduced GLMM excluding 'Sex' as a fixed effect
glme_reduced_Sex = fitglme(data,'L ~ 1 + Capture_Year + Month + Site + (1|Animal)','Distribution','Normal','Link','log');
% Compare the reduced model with the full model to test the significance of 'Sex'
result_Sex = compare(glme_reduced_Sex, glme)

% Fit a reduced GLMM excluding 'Month' as a fixed effect
glme_reduced_Month = fitglme(data,'L ~ 1 + Capture_Year + Sex + Site + (1|Animal)','Distribution','Normal','Link','log');
% Compare the reduced model with the full model to test the significance of 'Month'
result_Month = compare(glme_reduced_Month, glme)

% Fit a reduced GLMM excluding 'Site' as a fixed effect
glme_reduced_Site = fitglme(data,'L ~ 1 + Capture_Year + Sex + Month + (1|Animal)','Distribution','Normal','Link','log');
% Compare the reduced model with the full model to test the significance of 'Site'
result_Site = compare(glme_reduced_Site, glme)

%%

% In MATLAB, you can obtain the estimated marginal means (also known as "least-squares means") 
% for all variables and categories in a generalized linear mixed-effects model (GLME) using the predict function. 
% Calculate marginal means manually for categorical predictors by creating new datasets where you hold 
% other variables constant and vary the one of interest.


% Extract unique levels
siteLevels = unique(data.Site);
yearLevels = unique(data.Capture_Year);
monthLevels = unique(data.Month);
sexLevels = unique(data.Sex);

% Initialize arrays to store the marginal means
marginalMeans_Site = zeros(length(siteLevels), 1);
marginalMeans_Year = zeros(length(yearLevels), 1);
marginalMeans_Month = zeros(length(monthLevels), 1);
marginalMeans_Sex = zeros(length(sexLevels), 1);

% Calculate marginal means for 'Site'
for i = 1:length(siteLevels)
    newdata = data;
    newdata.Site(:) = siteLevels(i);  % Set 'Site' to each level
    marginalMeans_Site(i) = mean(predict(glme, newdata));
end

% Calculate marginal means for 'Capture_Year'
for j = 1:length(yearLevels)
    newdata = data;
    newdata.Capture_Year(:) = yearLevels(j);  % Set 'Capture_Year' to each level
    marginalMeans_Year(j) = mean(predict(glme, newdata));
end

% Calculate marginal means for 'Month'
for j = 1:length(monthLevels)
    newdata = data;
    newdata.Month(:) = monthLevels(j);  % Set 'Capture_Year' to each level
    marginalMeans_Month(j) = mean(predict(glme, newdata));
end

% Calculate marginal means for 'Sex'
for j = 1:length(sexLevels)
    newdata = data;
    newdata.Sex(:) = sexLevels(j);  % Set 'Capture_Year' to each level
    marginalMeans_Sex(j) = mean(predict(glme, newdata));
end

% Display the marginal means
disp(table(siteLevels, marginalMeans_Site, 'VariableNames', {'Site', 'MarginalMean'}));
disp(table(yearLevels, marginalMeans_Year, 'VariableNames', {'Capture_Year', 'MarginalMean'}));
disp(table(monthLevels, marginalMeans_Month, 'VariableNames', {'Month', 'MarginalMean'}));
disp(table(sexLevels, marginalMeans_Sex, 'VariableNames', {'Sex', 'MarginalMean'}));

%%

glme_nosite_mpl = fitglme(data,'L ~ 1 + Capture_Year + Sex + Month  + (1|Animal)','Distribution','Normal','Link','log','FitMethod','MPL');

glme_nosite_rempl = fitglme(data,'L ~ 1 + Capture_Year + Sex + Month  + (1|Animal)','Distribution','Normal','Link','log','FitMethod','REMPL');


%%

model_site = fitglme(data,'L ~ 1 + Capture_Year + Sex + Month  + Site + (1|Animal)','Distribution','Normal','Link','log','FitMethod','REMPL');

res = residuals(model_site); % Get residuals
months = data.Month;    % Grouping variable
CY = data.Capture_Year;
sex = data.Sex;
site = data.Site;

figure();
boxplot(res, months);
xlabel('Month');
ylabel('Residual');
title('Residuals by Month');


figure();
boxplot(res, CY);
xlabel('Capture Year');
ylabel('Residual');
title('Residuals by Capture Year');

figure();
boxplot(res, sex);
xlabel('Sex');
ylabel('Residual');
title('Residuals by Sex');

figure();
boxplot(res, site);
xlabel('Site');
ylabel('Residual');
title('Residuals by Site');

%%

model = fitglme(data,'L ~ 1 + Capture_Year + Sex + Month  + (1|Animal)','Distribution','Normal','Link','log','FitMethod','REMPL');

res = residuals(model); % Get residuals
months = data.Month;    % Grouping variable
CY = data.Capture_Year;
sex = data.Sex;
site = data.Site;

figure();
boxplot(res, months);
xlabel('Month');
ylabel('Residual');
title('Residuals by Month');


figure();
boxplot(res, CY);
xlabel('Capture Year');
ylabel('Residual');
title('Residuals by Capture Year');

figure();
boxplot(res, sex);
xlabel('Sex');
ylabel('Residual');
title('Residuals by Sex');

figure();
boxplot(res, site);
xlabel('Site');
ylabel('Residual');
title('Residuals by Site');

%% 

data_downsampled = data(1:2:end,:);

glme_downsample = fitglme(data_downsampled,'L ~ 1 + Capture_Year + Sex + Month  + (1|Animal)','Distribution','Normal','Link','log','FitMethod','REMPL');

%%

n = height(data);
idx = randperm(n, round(0.8 * n)); % Randomly select 80% of indices
data_downsampled_80 = data(idx, :);

glme_downsample_80 = fitglme(data_downsampled_80,'L ~ 1 + Capture_Year + Sex + Month  + (1|Animal)','Distribution','Normal','Link','log','FitMethod','REMPL');

%%

n = height(data);
idx = randperm(n, round(0.9 * n)); % Randomly select 90% of indices
data_downsampled_90 = data(idx, :);

glme_downsample_90 = fitglme(data_downsampled_90,'L ~ 1 + Capture_Year + Sex + Month  + (1|Animal)','Distribution','Normal','Link','log','FitMethod','REMPL');


%%

n = height(data);
idx = randperm(n, round(0.92 * n)); % Randomly select 99% of indices
data_downsampled_92 = data(idx, :);

glme_downsample_92 = fitglme(data_downsampled_92,'L ~ 1 + Capture_Year + Sex + Month  + (1|Animal)','Distribution','Normal','Link','log','FitMethod','REMPL');


%%

% The correct model to use (renamed model below)
glme_site_animal = fitglme(data,'L ~ 1 + Capture_Year + Sex + Month + (1|Site:Animal)','Distribution','Normal','Link','log');


glme_site_siteanimal = fitglme(data,'L ~ 1 + Capture_Year + Sex + Month + (1|Site) + (1|Site:Animal)','Distribution','Normal','Link','log');

glme_model_site = fitglme(data,'L ~ 1 + Capture_Year + Sex + Month + (1|Site)','Distribution','Normal','Link','log');

%%

% The correct model to use:
model = fitglme(data,'L ~ 1 + Capture_Year + Sex + Month + (1|Site:Animal)','Distribution','Normal','Link','log');

%% Compute variances

% This is the "Variance explained by fixed-effect" and "Variance explained by random effect".
% MATLAB does not seem to compute this for the output.
% You can think of it this way: V_total = V_f + V_r + sigma^2  
%
% V_f: fixed-effect variance
% Intuitively, it captures how much of the variability in the fitted log-response is explained by the predictors Capture_Year, Month, and Sex.
% 
% V_r: random-effect variance
% This quantifies the contribution of the random intercepts (Animal-within-Site, or Site and Animal if you use nesting) to the overall variability in log-L.
% 
% sigma^2:   This is just model.Dispersion
% It accounts for the remaining unexplained variability after fitting both fixed and random effects.


yhat_fixed = predict(model, data, 'Conditional', false);
% Variance of the fixed‐only predictor
V_f = var(yhat_fixed);

% “Conditional” predictions (fixed + random effects)
yhat_full  = predict(model, data, 'Conditional', true);
% Variance due to random effects: difference between full and fixed‐only
V_r = var(yhat_full - yhat_fixed);

V_total = V_f + V_r + model.Dispersion;

% Marginal R2 is concerned with variance explained by fixed factors
R_2_marginal = V_f/V_total;

% conditional R2 is concerned with variance explained by both fixed and random factors.
R_2_conditional  =  (V_f+V_r)/V_total;

%% Plot residuals

res = residuals(model); % Get residuals

figure(); 
plotResiduals(model,'fitted');

figure();
plotResiduals(model,'histogram');

figure();
boxplot(res, data.Month);
xlabel('Month');
ylabel('Residual');
title('Residuals by Month');


figure();
boxplot(res, data.Capture_Year);
xlabel('Capture Year');
ylabel('Residual');
title('Residuals by Capture Year');

figure();
boxplot(res, data.Sex);
xlabel('Sex');
ylabel('Residual');
title('Residuals by Sex');

% figure();
% boxplot(res, site);
% xlabel('Site');
% ylabel('Residual');
% title('Residuals by Site');
%%
figure()

subplot(2,2,1)
plotResiduals(model,'fitted');
set(gca,'TickLabelInterpreter','latex')

subplot(2,2,2)
plotResiduals(model,'histogram');
set(gca,'TickLabelInterpreter','latex')

% subplot(3,2,3)
% boxplot(res, data.Capture_Year);
% xlabel('Capture Year',Interpreter='latex');
% ylabel('Residual',Interpreter='latex');
% title('Residuals by Capture Year',Interpreter='latex');
% set(gca,'TickLabelInterpreter','latex')

subplot(2,2,3)
boxplot(res, data.Month);
xlabel('Month',Interpreter='latex');
ylabel('Residual',Interpreter='latex');
title('Residuals by Month',Interpreter='latex');
set(gca,'TickLabelInterpreter','latex')

subplot(2,2,4)
boxplot(res, data.Sex);
xlabel('Sex',Interpreter='latex');
ylabel('Residual',Interpreter='latex');
title('Residuals by Sex',Interpreter='latex');
set(gca,'TickLabelInterpreter','latex')

% title(t,'Mean: (2688, 1863)','fontweight','bold',fontsize = 15, Interpreter='latex')
% title(t,'Residuals','fontweight','bold',fontsize = 15, Interpreter='latex')
latex_fig(15, 5, 4)

%%

% The correct model to use:
model = fitglme(data,'L ~ 1 + Month + Sex + (Capture_Year|Site) + (1|Site:Animal)','Distribution','Normal','Link','log')

% model = fitglme(data,'L ~ 1 + Sex + Capture_Year + (Month|Site) + (1|Site:Animal)','Distribution','Normal','Link','log')
model.Rsquared

% Compute variances

% This is the "Variance explained by fixed-effect" and "Variance explained by random effect".
% MATLAB does not seem to compute this for the output.
% You can think of it this way: V_total = V_f + V_r + sigma^2  
%
% V_f: fixed-effect variance
% Intuitively, it captures how much of the variability in the fitted log-response is explained by the predictors Capture_Year, Month, and Sex.
% 
% V_r: random-effect variance
% This quantifies the contribution of the random intercepts (Animal-within-Site, or Site and Animal if you use nesting) to the overall variability in log-L.
% 
% sigma^2:   This is just model.Dispersion
% It accounts for the remaining unexplained variability after fitting both fixed and random effects.


yhat_fixed = predict(model, data, 'Conditional', false);
% Variance of the fixed‐only predictor
V_f = var(yhat_fixed);

% “Conditional” predictions (fixed + random effects)
yhat_full  = predict(model, data, 'Conditional', true);
% Variance due to random effects: difference between full and fixed‐only
V_r = var(yhat_full - yhat_fixed);

V_total = V_f + V_r + model.Dispersion;

% Marginal R2 is concerned with variance explained by fixed factors
R_2_marginal = V_f/V_total


% conditional R2 is concerned with variance explained by both fixed and random factors.
R_2_conditional  =  (V_f+V_r)/V_total
