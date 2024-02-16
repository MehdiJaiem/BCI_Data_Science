
clc;
clear;
close all;
%load('task_Day3to5.mat')
% defect electrode PO7
%dataEEG.label(7)=[];
load('task_Day3to5_Part2.mat'); %P7 deffect hence 23 channels
%Sampling frequency
fs = 250;
% Operating workingChannels
workingChannels = [7, 8, 9, 12,13,14, 15, 17, 18,22, 23]';

% Data extraction
BT_WS = burning_trials(workingChannels,:,:);
CT_WS = control_trials(workingChannels,:,:);
ET_WS = explosion_trials(workingChannels,:,:);

BTM_O2 = mean(BT_WS(3,:,:),3);
CTM_O2 = mean(control_trials(3,:,:),3);
ETM_O2 = mean(explosion_trials(3,:,:),3);

BT_std_O2 = std(BTM_O2);
CT_std_O2 = std(CTM_O2);
ET_std_O2 = std(ETM_O2);

%% Variance calculation for all channels
%Variance over all channels

%Burning trials variance across all rows(trials)
VAR_BT = var(BT_WS,'',2); 
VAR_BT=squeeze(VAR_BT);
VAR_BT_MEAN=mean(VAR_BT(:,:),2);
VAR_BT = VAR_BT';

%Control trials variance across all rows(trials)
VAR_CT = var(CT_WS,'',2);
VAR_CT=squeeze(VAR_CT);
VAR_CT_MEAN=mean(VAR_CT(:,:),2);
VAR_CT = VAR_CT';

%Explosion event variance across all rows(trials) 
VAR_ET = var(ET_WS,'',2);
VAR_ET= squeeze(VAR_ET);
VAR_ET_MEAN=mean(VAR_ET(:,:),2);
VAR_ET = VAR_ET';



%% PSD
    fig4=figure;
    
% [pwelchBT,fBT] = pwelch(BT_WS,[],[],[],fs); %squeezing does not make
% sense then

% pwelchBT=pwelchBT(1:30,:);
% PSD_BT_MEAN=mean(pwelchBT(:,:),2);
    
    [PSDB,fB] = pwelch(squeeze(BT_WS(:,:,:)),[],[],[],fs);
    [PSDC,fC]=pwelch(squeeze(CT_WS(:,:,:)),[],[],[],fs);
    [PSDE,fE]=pwelch(squeeze(ET_WS(:,:,:)),[],[],[],fs);



    PSDB=PSDB(1:30,:);
    fB=fB(1:30,:);
    PSDB_MEAN=mean(PSDB(:,:),2);
    PSDC=PSDC(1:30,:);
    fC=fC(1:30,:);
    PSDC_MEAN=mean(PSDC(:,:),2);
    PSDE=PSDE(1:30,:);
    fE=fE(1:30,:);
    PSDE_MEAN=mean(PSDE(:,:),2);
    

for j=1:11
    for i=1:35
PSDB_n(j,:,i) = pwelch(BT_WS(j,:,i),[],[],[],fs);
    end
    for i=1:158
PSDC_n(j,:,i) = pwelch(CT_WS(j,:,i),[],[],[],fs);
    end
    for i=1:38
PSDE_n(j,:,i) = pwelch(ET_WS(j,:,i),[],[],[],fs);
    end
end

PSDE_n=PSDE_n(:,1:40,:)
PSDB_n=PSDB_n(:,1:40,:);
PSDC_n=PSDC_n(:,1:40,:)
B_psd = zeros(440,35);
 C_psd = zeros(440,158);
 E_psd = zeros(440,38);
for i=1:11
    k = 40*(i-1);
for j=1:40
    B_psd(j+k,:) = PSDB_n(i,j,:); 
    C_psd(j+k,:) = PSDC_n(i,j,:);
    E_psd(j+k,:) = PSDE_n(i,j,:); 
end
end
    
     


%% DWT - Approximation coeffcient

for m=1:11
    %burning trials - Aproximation coefficient
    for i=1:35
    [Col1(m,:,i),Line1(m,:,i)] = wavedec(BT_WS(m,:,i),3,'db8');
    A1(m,:,i) = appcoef(Col1(m,:,i),Line1(m,:,i),'db8');
    end
    %control trials - Aproximation coefficient
    for i=1:158
    [Col2(m,:,i),Line2(m,:,i)] = wavedec(CT_WS(m,:,i),3,'db8');
    A2(m,:,i) = appcoef(Col2(m,:,i),Line2(m,:,i),'db8');
    end
    %explosion trials - Aproximation coefficient
    for i=1:38
    [Col3(m,:,i),Line3(m,:,i)] = wavedec(ET_WS(m,:,i),3,'db8');
    A3(m,:,i) = appcoef(Col3(m,:,i),Line3(m,:,i),'db8');
    end
end
 B_dwt = zeros(660,35);
 C_dwt = zeros(660,158);
 E_dwt = zeros(660,38);
for i=1:11
    k = 60*(i-1);
for j=1:60
    B_dwt(j+k,:) = A1(i,j,:); 
    C_dwt(j+k,:) = A2(i,j,:);
    E_dwt(j+k,:) = A3(i,j,:); 
end
end





%% Machine Learning SVM

E_psd=E_psd';
E_dwt=E_dwt';
B_psd=B_psd';
B_dwt=B_dwt';
C_psd=C_psd';
C_dwt=C_dwt';

%explode_feat = [VAR_ET, E_psd, E_dwt];
explode_feat = [VAR_ET, E_psd, E_dwt];
control_feat = [VAR_CT, C_psd, C_dwt];

%explode_feat = E_psd;
%burn_feat = B_psd;
%control_feat = C_psd;

%control_feat_mat = features(event, filtered_O2, 250, 'control');
%burning_feat = features(event, filtered_O2, 250, 'burning'); 

X = [explode_feat; control_feat];
X(:,2) = [];
Y = cell(size(X,1), 1);
Y(1:size(explode_feat, 1)) = {'explosion'};
Y(size(explode_feat, 1) + 1:end) = {'control'};


%% Shuffle Data

permutations = randperm(size(X, 1));

randomized_X = X(permutations, :);
randomized_Y = Y(permutations, :);

%% Cross Validation
indices = crossvalind('HoldOut', size(randomized_X, 1), 0.2);
onez = find(indices==1);
zeros = find(indices==0); 
X_train = randomized_X(onez, :);

Y_train = randomized_Y(onez, :);
X_test = randomized_X(zeros, :);
Y_test = randomized_Y(zeros, :);
%% ML

SVMModel = fitcsvm(X_train,Y_train,'Standardize',true, 'KernelFunction','linear','ClassNames',{'explosion','control'});

%sv = SVMModel.SupportVectors;
%
predictions = predict(SVMModel, X_test);

c = confusionchart(Y_test, predictions);


% tp = c.NormalizedValues(1,1)
% fp = c.NormalizedValues(1,2)
% tn = c.NormalizedValues(2,2)
% fn = c.NormalizedValues(2,1)
%         
% sens_svm= tp/(tp+fp);% accuracy
% acc_svm= (tp+tn)/(tp+fp+fn+tn);
% spec_svm= tn/(tn+fn);% accuracy

%% K FOLD

label_E = ones(1,38).*(-1); % Explosion-control label

label_C = ones(1,158); % Explosion-control label

label_B = ones(1,35).*(-1); % Burn-control label
%% transpose

VAR_BT=VAR_BT';
VAR_CT=VAR_CT';
VAR_ET=VAR_ET';

B_psd=B_psd';
C_psd=C_psd';
E_psd=E_psd';

B_dwt=B_dwt';
C_dwt=C_dwt';
E_dwt=E_dwt';

burn_feat = [VAR_BT; B_psd; B_dwt];
B_dvp=burn_feat;
explode_feat = [VAR_ET; E_psd; E_dwt];
E_dvp=explode_feat;
control_feat = [VAR_CT; C_psd; C_dwt];
C_dvp=control_feat;



%% %%%%%%%%% Burn VS Control Crossvalidation data partition %%%%%%%%%%%%



[BCvar_train, BCvar_traindata, BCvar_trainlabel, BCvar_test, BCvar_testdata, ...
    BCvar_testlabel]= datapartition(VAR_BT, label_B, VAR_CT, label_C);

[BCpsd_train, BCpsd_traindata, BCpsd_trainlabel, BCpsd_test, BCpsd_testdata, ...
    BCpsd_testlabel]= datapartition(B_psd, label_B, C_psd, label_C);

[BCdwt_train, BCdwt_traindata, BCdwt_trainlabel, BCdwt_test, BCdwt_testdata, ...
    BCdwt_testlabel]= datapartition(B_dwt, label_B, C_dwt, label_C);

[BCdvp_train, BCdvp_traindata, BCdvp_trainlabel, BCdvp_test, BCdvp_testdata, ...
    BCdvp_testlabel]= datapartition(B_dvp, label_B, C_dvp, label_C);


%% %%%%%%%%% Explosion VS Control Crossvalidation data partition %%%%%%%%%%%%

[ECvar_train, ECvar_traindata, ECvar_trainlabel, ECvar_test, ECvar_testdata, ...
    ECvar_testlabel]= datapartition(VAR_ET, label_E, VAR_CT, label_C);

[ECpsd_train, ECpsd_traindata, ECpsd_trainlabel, ECpsd_test, ECpsd_testdata, ...
    ECpsd_testlabel]= datapartition(E_psd, label_E, C_psd, label_C);

[ECdwt_train, ECdwt_traindata, ECdwt_trainlabel, ECdwt_test, ECdwt_testdata, ...
    ECdwt_testlabel]= datapartition(E_dwt, label_E, C_dwt, label_C);

[ECdvp_train, ECdvp_traindata, ECdvp_trainlabel, ECdvp_test, ECdvp_testdata, ...
    ECdvp_testlabel]= datapartition(E_dvp, label_E, C_dvp, label_C);


%% %%%%%%%%% Training %%%%%%%%%%
%%%%%%% Burning-Control variance%%%%%%%%%%%

figSize = [2 2 15 10];
fig9 = figure;
set(fig9, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);

[BCvar_model, BCvar_acc, BCvar_sens, BCvar_spec, BCvaravg_acc, BCvaravg_sens, BCvaravg_spec, ...
    BCvarstd_acc, BCvarstd_sens, BCvarstd_spec, BCvar_indx] = kFoldCrossVal(BCvar_traindata,BCvar_trainlabel);

BCvar_pred = predict(BCvar_model, BCvar_testdata);
BCvar_mat=confusionmat(BCvar_testlabel, BCvar_pred);

BCvar_pred = predict(BCvar_model, BCvar_testdata);
BCvar_mat=confusionmat(BCvar_testlabel, BCvar_pred);
BCvar_chart= confusionchart(BCvar_testlabel, BCvar_pred);

Accuracy = [BCvar_acc'; BCvaravg_acc; BCvarstd_acc];
Sensitivity = [BCvar_sens'; BCvaravg_sens; BCvarstd_sens];
Specificity = [BCvar_spec'; BCvaravg_spec; BCvarstd_spec];
k_fold = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; 'Mean'; 'Std'};

Table_BCvar = table(Accuracy, Sensitivity, Specificity, 'RowNames',k_fold)
savefig(fig9);
print(fig9, '-dpdf','burn_cont_var');



%% 
% %%%%%%% Burning-Control PSD%%%%%%%%%%%
% 
figSize = [2 2 15 10];
fig10 = figure;
set(fig10, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);

[BCpsd_model, BCpsd_acc, BCpsd_sens, BCpsd_spec, BCpsdavg_acc, BCpsdavg_sens, BCpsdavg_spec,...
    BCpsdstd_acc, BCpsdstd_sens, BCpsdstd_spec, BCpsd_indx] = kFoldCrossVal(BCpsd_traindata,BCpsd_trainlabel);

BCpsd_pred = predict(BCpsd_model, BCpsd_testdata);
BCpsd_mat=confusionmat(BCpsd_testlabel, BCpsd_pred);
BCpsd_chart= confusionchart(BCpsd_testlabel, BCpsd_pred);

Accuracy = [BCpsd_acc'; BCpsdavg_acc; BCpsdstd_acc];
Sensitivity = [BCpsd_sens'; BCpsdavg_sens; BCpsdstd_sens];
Specificity = [BCpsd_spec'; BCpsdavg_spec; BCpsdstd_spec];
k_fold = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; 'Mean'; 'Std'};

Table_BCpsd = table(Accuracy, Sensitivity, Specificity, 'RowNames',k_fold)
savefig(fig10);
print(fig10, '-dpdf','burn_cont_psd');
%%
% %%%%%%% Burning-Control DWT%%%%%%%%%%%
% 
figSize = [2 2 15 10];
fig11 = figure;
set(fig11, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
[BCdwt_model, BCdwt_acc, BCdwt_sens, BCdwt_spec, BCdwtavg_acc, BCdwtavg_sens, BCdwtavg_spec,...
    BCdwtstd_acc, BCdwtstd_sens, BCdwtstd_spec, BCdwt_indx] = kFoldCrossVal(BCdwt_traindata,BCdwt_trainlabel);

BCdwt_pred = predict(BCdwt_model, BCdwt_testdata);
BCdwt_mat=confusionmat(BCdwt_testlabel, BCdwt_pred);
BCdwt_chart= confusionchart(BCdwt_testlabel, BCdwt_pred);

Accuracy = [BCdwt_acc'; BCdwtavg_acc; BCdwtstd_acc];
Sensitivity = [BCdwt_sens'; BCdwtavg_sens; BCdwtstd_sens];
Specificity = [BCdwt_spec'; BCdwtavg_spec; BCdwtstd_spec];
k_fold = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; 'Mean'; 'Std'};

Table_BCdwt = table(Accuracy, Sensitivity, Specificity, 'RowNames',k_fold)
savefig(fig11);
print(fig11, '-dpdf','burn_cont_dwt');
%%
% %%%%%%% Burning-Control DWT%%%%%%%%%%%
% figSize = [2 2 15 10];
fig12 = figure;
set(fig12, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
[BCdvp_model, BCdvp_acc, BCdvp_sens, BCdvp_spec, BCdvpavg_acc, BCdvpavg_sens, BCdvpavg_spec,...
    BCdvpstd_acc, BCdvpstd_sens, BCdvpstd_spec, BCdvp_indx] = kFoldCrossVal(BCdvp_traindata,BCdvp_trainlabel);

BCdvp_pred = predict(BCdvp_model, BCdvp_testdata);
BCdvp_mat=confusionmat(BCdvp_testlabel, BCdvp_pred);
BCdvp_chart= confusionchart(BCdvp_testlabel, BCdvp_pred);

Accuracy = [BCdvp_acc'; BCdvpavg_acc; BCdvpstd_acc];
Sensitivity = [BCdvp_sens'; BCdvpavg_sens; BCdvpstd_sens];
Specificity = [BCdvp_spec'; BCdvpavg_spec; BCdvpstd_spec];
k_fold = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; 'Mean'; 'Std'};

Table_BCdvp = table(Accuracy, Sensitivity, Specificity, 'RowNames',k_fold)
savefig(fig12);
print(fig12, '-dpdf','explode_cont_dvp');

%% %%%%%%%%% Training %%%%%%%%%%
%%%%%%% Explosion-Control variance%%%%%%%%%%%

figSize = [2 2 15 10];
fig14 = figure;
set(fig14, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);

[ECvar_model, ECvar_acc, ECvar_sens, ECvar_spec, ECvaravg_acc, ECvaravg_sens, ECvaravg_spec, ...
    ECvarstd_acc, ECvarstd_sens, ECvarstd_spec, ECvar_indx] = kFoldCrossVal(ECvar_traindata,ECvar_trainlabel);

ECvar_pred = predict(ECvar_model, ECvar_testdata);
ECvar_mat=confusionmat(ECvar_testlabel, ECvar_pred);

ECvar_pred = predict(ECvar_model, ECvar_testdata);
ECvar_mat=confusionmat(ECvar_testlabel, ECvar_pred);
ECvar_chart= confusionchart(ECvar_testlabel, ECvar_pred);

Accuracy = [ECvar_acc'; ECvaravg_acc; ECvarstd_acc];
Sensitivity = [ECvar_sens'; ECvaravg_sens; ECvarstd_sens];
Specificity = [ECvar_spec'; ECvaravg_spec; ECvarstd_spec];
k_fold = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; 'Mean'; 'Std'};

Table_ECvar = table(Accuracy, Sensitivity, Specificity, 'RowNames',k_fold)
savefig(fig14);
print(fig14, '-dpdf','explode_cont_var');



%% 
% %%%%%%% Explosion-Control PSD%%%%%%%%%%%
% 
figSize = [2 2 15 10];
fig15 = figure;
set(fig15, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);

[ECpsd_model, ECpsd_acc, ECpsd_sens, ECpsd_spec, ECpsdavg_acc, ECpsdavg_sens, ECpsdavg_spec,...
    ECpsdstd_acc, ECpsdstd_sens, ECpsdstd_spec, ECpsd_indx] = kFoldCrossVal(ECpsd_traindata,ECpsd_trainlabel);

ECpsd_pred = predict(ECpsd_model, ECpsd_testdata);
ECpsd_mat=confusionmat(ECpsd_testlabel, ECpsd_pred);
ECpsd_chart= confusionchart(ECpsd_testlabel, ECpsd_pred);

Accuracy = [ECpsd_acc'; ECpsdavg_acc; ECpsdstd_acc];
Sensitivity = [ECpsd_sens'; ECpsdavg_sens; ECpsdstd_sens];
Specificity = [ECpsd_spec'; ECpsdavg_spec; ECpsdstd_spec];
k_fold = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; 'Mean'; 'Std'};

Table_ECpsd = table(Accuracy, Sensitivity, Specificity, 'RowNames',k_fold)
savefig(fig15);
print(fig15, '-dpdf','explode_cont_psd');
 %%
% %%%%%%% Explosion-Control DWT%%%%%%%%%%%
% 
figSize = [2 2 15 10];
fig16 = figure;
set(fig16, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
[ECdwt_model, ECdwt_acc, ECdwt_sens, ECdwt_spec, ECdwtavg_acc, ECdwtavg_sens, ECdwtavg_spec,...
    ECdwtstd_acc, ECdwtstd_sens, ECdwtstd_spec, ECdwt_indx] = kFoldCrossVal(ECdwt_traindata,ECdwt_trainlabel);

ECdwt_pred = predict(ECdwt_model, ECdwt_testdata);
ECdwt_mat=confusionmat(ECdwt_testlabel, ECdwt_pred);
ECdwt_chart= confusionchart(ECdwt_testlabel, ECdwt_pred);

Accuracy = [ECdwt_acc'; ECdwtavg_acc; ECdwtstd_acc];
Sensitivity = [ECdwt_sens'; ECdwtavg_sens; ECdwtstd_sens];
Specificity = [ECdwt_spec'; ECdwtavg_spec; ECdwtstd_spec];
k_fold = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; 'Mean'; 'Std'};

Table_ECdwt = table(Accuracy, Sensitivity, Specificity, 'RowNames',k_fold)
savefig(fig16);
print(fig16, '-dpdf','explode_cont_dwt');

% %%%%%%% Explosion-Control DVP%%%%%%%%%%%
% 
figSize = [2 2 15 10];
fig20= figure;
set(fig20, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
[ECdvp_model, ECdvp_acc, ECdvp_sens, ECdvp_spec, ECdvpavg_acc, ECdvpavg_sens, ECdvpavg_spec,...
    ECdvpstd_acc, ECdvpstd_sens, ECdvpstd_spec, ECdvp_indx] = kFoldCrossVal(ECdvp_traindata,ECdvp_trainlabel);

ECdvp_pred = predict(ECdvp_model, ECdvp_testdata);
ECdvp_mat=confusionmat(ECdvp_testlabel, ECdvp_pred);
ECdvp_chart= confusionchart(ECdvp_testlabel, ECdvp_pred);

Accuracy = [ECdvp_acc'; ECdvpavg_acc; ECdvpstd_acc]
Sensitivity = [ECdvp_sens'; ECdvpavg_sens; ECdvpstd_sens]
Specificity = [ECdvp_spec'; ECdvpavg_spec; ECdvpstd_spec]
k_fold = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; 'Mean'; 'Std'};

Table_ECdvp = table(Accuracy, Sensitivity, Specificity, 'RowNames',k_fold)
savefig(fig20);
print(fig20, '-dpdf','explode_cont_dvp');

%% Conclusion Table
feature = {'Burning (VAR)'; 'Burning (PSD)'; 'Burning (DWT)'; 'Burning (DVP)';...
    'Explosion (VAR)'; 'Explosion (PSD)'; 'Explosion (DWT)';...
    'Explosion (DVP)'};
Accuracy = [BCvaravg_acc BCpsdavg_acc BCdwtavg_acc BCdvpavg_acc ECvaravg_acc ECpsdavg_acc ECdwtavg_acc ECdvpavg_acc]';
Sensitivity= [BCvaravg_sens BCpsdavg_sens BCdwtavg_sens BCdvpavg_sens ECvaravg_sens ECpsdavg_sens ECdwtavg_sens ECdvpavg_sens ]';
Specificity= [ BCvaravg_spec BCpsdavg_spec BCdwtavg_spec BCdvpavg_spec ECvaravg_spec ECpsdavg_spec ECdwtavg_spec ECdvpavg_spec ]';
Table_compare = table(Accuracy, Sensitivity, Specificity, 'RowNames', feature)

