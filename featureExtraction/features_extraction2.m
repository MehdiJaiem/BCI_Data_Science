

fs = 250;

Channels = [5, 7, 8, 9, 14, 15, 17, 22, 23]';

burning_t = burning_trials(Channels,:,:);
control_t = control_trials(Channels,:,:);
explosion_t = explosion_trials(Channels,:,:);

%% Variance


burn_var = var(burning_t,'',2);
cont_var = var(control_t,'',2);
expl_var = var(explosion_t,'',2);


B_var = squeeze(burn_var);
C_var = squeeze(cont_var);
E_var = squeeze(expl_var);

%% Epoching

burn_mean = mean(burning_trials(5,:,:),3);
cont_mean = mean(control_trials(5,:,:),3);
expl_mean = mean(explosion_trials(5,:,:),3);

burn_std = std(burn_mean);
cont_std = std(cont_mean);
expl_std = std(expl_mean);

%% %%%%%%% std plotting %%%%%%%%%% 
%%%%%%% Explosion vs. Control
figSize = [2 2 17 12];
fig1 = figure;
set(fig1, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
x = 1:376;
x2 = [x, fliplr(x)];

curve1 = expl_mean + expl_std;
curve2 = expl_mean - expl_std;
inBetween = [curve1 fliplr(curve2)];
fill(x2, inBetween, 'r', 'facealpha', 0.3, 'edgecolor', 'none');
hold on
plot(expl_mean, 'r');
curve1 = cont_mean + cont_std;
curve2 = cont_mean - cont_std;
inBetween = [curve1 fliplr(curve2)];
fill(x2, inBetween, 'b', 'facealpha', 0.3, 'edgecolor', 'none');
plot(cont_mean, 'b');
axis ([0 376 -15 20]);
hold off
xlabel('Time (s)');
ylabel('Magintude (u.a)');
set(gca,'FontSize', 14);
title('Average on trials for explosion data on Channel O1');
legend('std of exploding box event','mean of exploding box','std of control',' mean of control');
savefig(fig1);
print(fig1, '-dpdf','epoching_expl');
%%%%%%% Burning vs. Control

fig2 = figure;
set(fig2, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
curve1 = burn_mean + burn_std;
curve2 = burn_mean - burn_std;
inBetween = [curve1 fliplr(curve2)];
fill(x2, inBetween, 'r', 'facealpha', 0.3, 'edgecolor', 'none');
hold on
plot(burn_mean, 'r');
curve1 = cont_mean + cont_std;
curve2 = cont_mean - cont_std;
inBetween = [curve1 fliplr(curve2)];
fill(x2, inBetween, 'b', 'facealpha', 0.3, 'edgecolor', 'none');
plot(cont_mean, 'b');
axis ([0 376 -15 20]);
hold off
xlabel('Time (s)');
ylabel('Magintude (u.a)');
set(gca,'FontSize', 14);
title('Average on trials for burning data on Channel O1');
legend('std of burning box event','mean of burning box','std of control',' mean of control');
savefig(fig2);
print(fig2, '-dpdf','epoching_burn');

%%
%%%%% Plot variance

figSize = [2 2 17 20];
fig3 = figure;

set(fig3, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);


subplot(2,1,1)
scatter((1:9),C_var,'b','filled');
xlabel('Channels');
ylabel('Variance (u.a)');
grid on;
hold on
title('Burning / control variance');
scatter((1:9),B_var,'r','filled');
legend('Control trials','Burning trials');
set(gca,'FontSize', 14);

subplot(2,1,2)
scatter((1:9),C_var,'b','filled');
xlabel('Channels');
ylabel('Variance (u.a)');
hold on
scatter((1:9),E_var,'r','filled');
grid on;
title('Explosion / control variance');
legend('Control trials','Explosion trials');
set(gca,'FontSize', 14);

savefig(fig3);
print(fig3, '-dpng','variance');
%% %%%%%%%%%%%%%% PSD %%%%%%%%%%%%%%%

for j=1:9
    for i=1:35
burn_psd_long(j,:,i) = pwelch(burning_t(j,:,i),[],[],[],fs);
    end
    for i=1:158
cont_psd_long(j,:,i) = pwelch(control_t(j,:,i),[],[],[],fs);
    end
    for i=1:38
expl_psd_long(j,:,i) = pwelch(explosion_t(j,:,i),[],[],[],fs);
    end
end
burn_psd = burn_psd_long(:,1:40,:);
cont_psd = cont_psd_long(:,1:40,:);
expl_psd = expl_psd_long(:,1:40,:);

%%

%%%%% plot PSD

figSize = [2 2 17 20];
fig4 = figure;
set(fig4, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);


subplot(2,1,1)
scatter((1:40),burn_psd(:,:,20),'r');
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
grid on;
hold on
scatter((1:40),cont_psd(:,:,20),'b');
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/Hz)');
set(gca,'FontSize', 14);
title('20^{th} trial PSD of burning and control');
legend('Burning trials','Control trials');

subplot(2,1,2)
scatter((1:40),expl_psd(:,:,20),'r');
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/Hz)');
grid on;
hold on
scatter((1:40),cont_psd(:,:,20),'b');
set(gca,'FontSize', 14);
title('20^{th} trial PSD of explosion and control');
legend('Explosion trials','Control trials');


savefig(fig4);
print(fig4, '-dpng','psd');
%% %%%%%%%%%% DWT %%%%%%%%% 

for j=1:9
    for i=1:35
    [C1(j,:,i),L1(j,:,i)] = wavedec(burning_trials(Channels(j),:,i),3,'db8');
    A_burn(j,:,i) = appcoef(C1(j,:,i),L1(j,:,i),'db8');
    end
  
    for i=1:158
    [C2(j,:,i),L2(j,:,i)] = wavedec(control_trials(Channels(j),:,i),3,'db8');
    A_cont(j,:,i) = appcoef(C2(j,:,i),L2(j,:,i),'db8');
    end
    for i=1:38
    [C3(j,:,i),L3(j,:,i)] = wavedec(explosion_trials(Channels(j),:,i),3,'db8');
    A_expl(j,:,i) = appcoef(C3(j,:,i),L3(j,:,i),'db8');
    end
end
%%
%%%%% plot DWT

figSize = [2 2 17 20];
fig5 = figure;
set(fig5, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);


subplot(2,1,1)
scatter(1:9,A_burn(:,:,20),'r','filled');
hold on
ylabel('Approximate coefficients (u.a)');
xlabel('Channels');
grid on;
scatter(1:9,A_cont(:,:,20),'b','filled');
title('20^{th} trial DWT approximation over channels');
legend('Burning trials', 'Control trials');
axis([1 9 -80 60]);
set(gca,'FontSize', 14);

subplot(2,1,2)
scatter(1:9,A_expl(:,:,20),'r','filled');
hold on
xlabel('Channels');
ylabel('Approximate coefficients (u.a)');
grid on;
scatter(1:9,A_cont(:,:,20),'b','filled');
hold off;
legend('Explosion trials','Control trials');
set(gca,'FontSize', 14);

savefig(fig5);
print(fig5, '-dpng','dwt');
%% %%%%%%%%%%% Concatenate Matrices %%%%%%%%%%%%%
 

 B_psd = zeros(360,35);
 C_psd = zeros(360,158);
 E_psd = zeros(360,38);
for i=1:9
    k = 40*(i-1);
for j=1:40
    B_psd(j+k,:) = burn_psd(i,j,:); 
    C_psd(j+k,:) = cont_psd(i,j,:);
    E_psd(j+k,:) = expl_psd(i,j,:); 
end
end

%%%%%%%%%%%%%%%%%%%%%
 B_dwt = zeros(540,35);
 C_dwt = zeros(540,158);
 E_dwt = zeros(540,38);
for i=1:9
    k = 60*(i-1);
for j=1:60
    B_dwt(j+k,:) = A_burn(i,j,:); 
    C_dwt(j+k,:) = A_cont(i,j,:);
    E_dwt(j+k,:) = A_expl(i,j,:); 
end
end

%% %%%%%%% Concatenation %%%%%%%%%%%%

Evar = permutation(E_var);
Epsd = permutation(E_psd);
Edwt = permutation(E_dwt);

Cvar = permutation(C_var);
Cpsd = permutation(C_psd);
Cdwt = permutation(C_dwt);

Bvar = permutation(B_var);
Bpsd = permutation(B_psd);
Bdwt = permutation(B_dwt);

%% %%%%%%% Label vector %%%%%%%%%
label_E = ones(1,35).*(-1); % Explosion-control label

label_C = ones(1,35); % Explosion-control label

label_B = ones(1,35).*(-1); % Burn-control label

%% %%%%%%%%% Crossvalidation data partition %%%%%%%%%%%%


[ECvar_train, ECvar_traindata, ECvar_trainlabel, ECvar_test, ECvar_testdata, ...
    ECvar_testlabel]= datapartition(Evar, label_E, Cvar, label_C);

[ECpsd_train, ECpsd_traindata, ECpsd_trainlabel, ECpsd_test, ECpsd_testdata, ...
    ECpsd_testlabel]= datapartition(Epsd, label_E, Cpsd, label_C);

[ECdwt_train, ECdwt_traindata, ECdwt_trainlabel, ECdwt_test, ECdwt_testdata, ...
    ECdwt_testlabel]= datapartition(Edwt, label_E, Cdwt, label_C);

[BCvar_train, BCvar_traindata, BCvar_trainlabel, BCvar_test, BCvar_testdata, ...
    BCvar_testlabel]= datapartition(Bvar, label_B, Cvar, label_C);

[BCpsd_train, BCpsd_traindata, BCpsd_trainlabel, BCpsd_test, BCpsd_testdata, ...
    BCpsd_testlabel]= datapartition(Bpsd, label_B, Cpsd, label_C);

[BCdwt_train, BCdwt_traindata, BCdwt_trainlabel, BCdwt_test, BCdwt_testdata, ...
    BCdwt_testlabel]= datapartition(Bdwt, label_B, Cdwt, label_C);


%% %%%%%%%%% Training %%%%%%%%%%
%%%%%%%Explosion-Control variance%%%%%%%%%%%

figSize = [2 2 15 10];
fig6 = figure;
set(fig6, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);

[ECvar_model, ECvar_acc, ECvar_sens, ECvar_spec, ECvaravg_acc, ECvaravg_sens, ECvaravg_spec, ...
    ECvarstd_acc, ECvarstd_sens, ECvarstd_spec, ECvar_indx] = kFoldCrossVal(ECvar_traindata,ECvar_trainlabel);

ECvar_pred = predict(ECvar_model, ECvar_testdata);
ECvar_mat=confusionmat(ECvar_testlabel, ECvar_pred);

%%

figSize = [2 2 15 10];
fig6 = figure;
set(fig6, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
set(gca,'FontSize', 18);
ECvar_chart= confusionchart(ECvar_testlabel, ECvar_pred);

Accuracy = [ECvar_acc'; ECvaravg_acc; ECvarstd_acc];
Sensitivity = [ECvar_sens'; ECvaravg_sens; ECvarstd_sens];
Specificity = [ECvar_spec'; ECvaravg_spec; ECvarstd_spec];
k_fold = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; 'Mean'; 'Std'};

Table_ECvar = table(Accuracy, Sensitivity, Specificity, 'RowNames',k_fold)

savefig(fig6);
print(fig6, '-dpdf','expl_cont_var');

%%
%%%%%%%Explosion-Control PSD%%%%%%%%%%%

figSize = [2 2 15 10];
fig7 = figure;
set(fig7, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
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
savefig(fig7);
print(fig7, '-dpdf','expl_cont_psd');

%%
%%%%%%%Explosion-Control DWT%%%%%%%%%%%

figSize = [2 2 15 10];
fig8 = figure;
set(fig8, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
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
savefig(fig8);
print(fig8, '-dpdf','expl_cont_dwt');
%%
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
BCvar_chart= confusionchart(BCvar_testlabel, BCvar_pred);

Accuracy = [BCvar_acc'; BCvaravg_acc; BCvarstd_acc];
Sensitivity = [BCvar_sens'; BCvaravg_sens; BCvarstd_sens];
Specificity = [BCvar_spec'; BCvaravg_spec; BCvarstd_spec];
k_fold = {'1'; '2'; '3'; '4'; '5'; '6'; '7'; '8'; '9'; '10'; 'Mean'; 'Std'};

Table_BCvar = table(Accuracy, Sensitivity, Specificity, 'RowNames',k_fold)
savefig(fig9);
print(fig9, '-dpdf','burn_cont_var');

%%
%%%%%%% Burning-Control PSD%%%%%%%%%%%

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
%%%%%%% Burning-Control DWT%%%%%%%%%%%

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
feature = {'Explosion-Control variance'; 'Explosion-Control PSD'; 'Explosion-Control DWT';...
    'Burning-Control variance'; 'Burning-Control PSD'; 'Burning-Control DWT'};
Mean_acc = [ECvaravg_acc ECpsdavg_acc ECdwtavg_acc BCvaravg_acc BCpsdavg_acc BCdwtavg_acc]';
Mean_sens = [ECvaravg_sens ECpsdavg_sens ECdwtavg_sens BCvaravg_sens BCpsdavg_sens BCdwtavg_sens]';
Mean_spec = [ECvaravg_spec ECpsdavg_spec ECdwtavg_spec BCvaravg_spec BCpsdavg_spec BCdwtavg_spec]';
Table_compare = table(Mean_acc, Mean_sens, Mean_spec, 'RowNames', feature)