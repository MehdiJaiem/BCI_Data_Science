clc;
clear;
close all;
load('task_Day3to5.mat')
% defect electrode PO7
dataEEG.label(7)=[];
%% Original Signal - No Filtering or Signal Processing
%sampling frequency
fs = 250;
figSize = [2 2 17 10];
%extracting and experimenting with 9th trial. 
Channel_O2_trial = dataEEG.trial{1,1}(9,:);
%time vector extraction
t= dataEEG.time{1,1}-dataEEG.time{1,1}(1); 
%Trial time signal plot
fig1 = figure;
set(fig1, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
plot(t,Channel_O2_trial)
xlim([0, 2055]);
grid on;
title('Not filtered signal');
xlabel('time (s)');
ylabel('Magnitude (a.u.)');
%PSD Plot
fig2 = figure;
set(fig2, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
%Welch's power spectral density estimate (PSD)
pwelch(Channel_O2_trial,[],[],[],fs);
title('PSD of O2 Electrode Signal');
savefig(fig1);
%T6_TS: trial 6 time signal
print(fig1, '-dpdf','T6_TS');
%trial 6 power signal density
savefig(fig2);
print(fig2, '-dpdf','T6_PSD');
%% Filtering
%Highpass Chebyshev using 9th trial data 
low = 1;
Rp = 3;
Rs = 60;
Wp = 1;
Ws = 1/10;
Nfreq = fs/2;
[n1,Ws1] = cheb2ord(Wp/Nfreq,Ws/Nfreq,Rp,Rs);
[b1,a1] = cheby2(n1,Rs,Ws1,'high');
filtered_O2=filtfilt(b1,a1,double(Channel_O2_trial));  

fig3 = figure;
set(fig3, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
subplot(2,1,1);
plot(t,filtered_O2);
xlim([0, 2050]);
grid on;
title('Electrode Signal post HP filtering using Chebyshev type 2');
xlabel('time (s)'); 
ylabel('Magnitude (a.u.)');

%Lowpass Chebyshev
high = 40;
Wp = high;
Ws = high+round(high/10);
[n2,Ws2] = cheb2ord(Wp/Nfreq,Ws/Nfreq,Rp,Rs);
[b2,a2] =cheby2(n2,Rs,Ws2,'low');
filtered_O2=filtfilt(b2,a2,filtered_O2);

subplot(2,1,2);
plot(t,filtered_O2);
xlim([0, 2050]);
grid on;
title('Electrode HP filtered Signal post LP filtering using Chebyshev type 2');
xlabel('time (s)');
ylabel('Magnitude (a.u.)');
savefig(fig3);
print(fig3, '-dpdf','chebyshev');

%% Filtering all data trials with cheby2 filter
%create a null matrix with the size of trial cell where we would implement
%the newly filtered signal
allData=dataEEG.trial{1,1};
allDataFiltered_Fz = zeros(size(allData,1),size(allData,2));
%each row represents a trial and therefore a for loop over 
%size(dataEEG.trial{1,1},1) is neessary
for k=1:size(allData,1)
    %High Pass Filtering
    allDataFiltered_Fz(k,:)=filtfilt(b1,a1,double(allData(k,:))); 
    %Low Pass Filtering of pre-High Pass filtered signal
    allDataFiltered_Fz(k,:)=filtfilt(b2,a2,allDataFiltered_Fz(k,:));
end

%% PSD 
fig4 = figure;
set(fig4, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
subplot(2,1,1);
pwelch(Channel_O2_trial,[],[],[],fs);
title('PSD-Original Signal');
subplot(2,1,2);
pwelch(filtered_O2,[],[],[],fs);
title('PSD-Chebyshev Filtered O2 Signal');
savefig(fig4);
print(fig4, '-dpdf','PSD original VS chebytchev type2');

%% CAR
avg = zeros(1,length(t));
%size of the original signal
CAR = zeros(23,length(t));
% compute the average of the signal at all EEG electrodes
for i=1:length(t)
    avg(i) = mean(dataEEG.trial{1,1}(:,i));
end
for i=1:length(t)
    %j is the number of electrodes/trials17
    for j=1:23
        %subtract the average from the EEG signal at every electrode 
        % for every time point
        CAR(j,i) = dataEEG.trial{1,1}(j,i)-avg(i);
    end
end

%% Filtering of CAR data for O2 and Fz
%Channel 17: Fz data channel without CAR
Fz_DataSignal = dataEEG.trial{1,1}(17,:);
Cz_DataSignal = dataEEG.trial{1,1}(16,:);
%O2 channel data post CAR
Avg_O2 = CAR(9,:);
%Fz data channel post CAR
avg_Fz = CAR(17,:);
avg_Cz = CAR(16,:);

%high pass Chebytchev on normal Fz (without CAR)
DataFiltered_Fz=filtfilt(b1,a1,double(Fz_DataSignal)); 
DataFiltered_Cz=filtfilt(b1,a1,double(Cz_DataSignal)); 
%low pass Chebytchev on previously high passfiltered Fz 
DataFiltered_Fz=filtfilt(b2,a2,DataFiltered_Fz);
DataFiltered_Cz=filtfilt(b2,a2,DataFiltered_Cz);

%high pass Chebytchev on post CAR O2 
filtered_Avg_O2=filtfilt(b1,a1,double(Avg_O2)); 
%low pass Chebytchev on previously high passfiltered post CAR O2  
filtered_Avg_O2=filtfilt(b2,a2,filtered_Avg_O2);

%high pass Chebytchev on post CAR Fz
filtered_Avg_Fz=filtfilt(b1,a1,double(avg_Fz));
%low pass Chebytchev on previously high passfiltered post CAR Fz 
filtered_Avg_Fz=filtfilt(b2,a2,filtered_Avg_Fz); 

%high pass Chebytchev on post CAR Czz
filtered_Avg_Cz=filtfilt(b1,a1,double(avg_Cz));
%low pass Chebytchev on previously high passfiltered post CAR Czz 
filtered_Avg_Cz=filtfilt(b2,a2,filtered_Avg_Cz); 
%% CAR plotting
fig5 = figure;
set(fig5, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
% Plot of the O2 signal without Filtering/ without CAR+Filtering
subplot(2,3,1)
plot(t,Channel_O2_trial); 
xlim([0, 2050]);
grid on;
xlabel('Time (s)');
ylabel('Magnitude (a.u.)');
subtitle('Original O2 signal');

% Plot of the O2 signal without CAR+Filtering
subplot(2,3,2)
plot(t,filtered_O2); 
xlim([0, 2050]);
grid on;
xlabel('Time (s)');
ylabel('Magnitude (a.u.)');
subtitle('Chebyshev filtered O2 Signal');
subplot(2,3,3)

% Plot of the O2 signal WITH CAR+Filtering
plot(t,filtered_Avg_O2); xlim([0, 2050]);
grid on;
xlabel('Time (s)');
ylabel('Magnitude (a.u.)');
subtitle('Chebyschev on avg O2 Signal');

% Plot of the Fz signal without Filtering/ without CAR+Filtering
subplot(2,3,4)
plot(t,Fz_DataSignal); 
xlim([0, 2050]);
grid on;
xlabel('Time (s)');
ylabel('Magnitude (a.u.)');
subtitle('Original Fz signal');

% Plot of the Fz signal without CAR+Filtering
subplot(2,3,5)
plot(t,DataFiltered_Fz); 
xlim([0, 2050]);
xlabel('Time (s)');
ylabel('Magnitude (a.u.)');
subtitle('Chebyshev on Fz signal');
grid on;
% Plot of the Fz signal WITH CAR+Filtering
subplot(2,3,6)
plot(t,filtered_Avg_Fz); 
xlim([0, 2050]);
grid on;
xlabel('Time (s)');
ylabel('Magnitude (a.u.)');
subtitle('Chebyshev on avg Fz Signal');

savefig(fig5);
print(fig5, '-dpdf','PP1');
%%
fig6 = figure;
set(fig6, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
% Plot of the O2 signal without Filtering/ without CAR+Filtering
subplot(2,2,1)
plot(t,Channel_O2_trial)
xlim([0, 2055]);
grid on;
title('Unfiltered O2 Electrode Signal');
xlabel('time (s)');
ylabel('Magnitude (a.u.)');
yl1=yline(-4906,'--','Bottom Amp');
yl1.LabelHorizontalAlignment = 'left';
subplot(2,2,2)
pwelch(Channel_O2_trial,[],[],[],fs);
hold on
plot(62.5,54.9644,'r*')
hold off
title('Unfiltered PSD of O2 Electrode Signal');
subplot(2,2,3)
plot(t,filtered_O2);
xlim([0, 2050]);
grid on;
title('O2 Electrode Signal filtering using Chebyshev type 2');
xlabel('time (s)');
ylabel('Magnitude (a.u.)');
subplot(2,2,4)
pwelch(filtered_O2,[],[],[],fs);
title('O2 PSD filtering using Chebyshev type 2');
xl1=xline(1,'--','1Hz');
xl1.LabelVerticalAlignment = 'bottom';
xline(40,'--','40Hz');

savefig(fig6);
print(fig6, '-dpdf','PP2');



fig7 = figure;
figSize = [2 2 17 10];
set(fig7, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);

subplot(1,3,1)
plot(t,filtered_Avg_O2); xlim([0, 2050]);
grid on;
xlabel('Time (s)');
ylabel('Magnitude (a.u.)');
subtitle('Chebyschev on avg O2 Signal');

subplot(1,3,2)
plot(t,filtered_Avg_Fz); 
xlim([0, 2050]);
grid on;
xlabel('Time (s)');
ylabel('Magnitude (a.u.)');
subtitle('Chebyshev on avg Fz Signal');

subplot(1,3,3)
plot(t,filtered_Avg_Cz); 
xlim([0, 2050]);
grid on;
xlabel('Time (s)');
ylabel('Magnitude (a.u.)');
subtitle('Chebyshev and CAR on avg Cz Signal');

savefig(fig7);
print(fig7, '-dpdf','PP3');


envenrrrrr=0;

%% Epoching
% epoching(-0.3,0.8,250,event,filtered_O2, 1);


% %% Feature Extraction
% exploding_feat_mat = features(event, filtered_O2, 250, 'exploding');
% control_feat_mat = features(event, filtered_O2, 250, 'control');
% burning_feat_mat = features(event, filtered_O2, 250, 'burning'); 
% 
% X = [burning_feat_mat; control_feat_mat];
% X(:,2) = [];
% Y = cell(size(X,1), 1);
% Y(1:size(burning_feat_mat, 1)) = {'burning'};
% Y(size(burning_feat_mat, 1) + 1:end) = {'control'};
% 
% %% Shuffle Data
% 
% permutations = randperm(size(X, 1));
% 
% randomized_X = X(permutations, :)
% randomized_Y = Y(permutations, :)
% 
% %% Cross Validation
% indices = crossvalind('HoldOut', size(randomized_X, 1), 0.2);
% kfold_set = find(indices==1);
% holdout_set = find(indices == 0);
% 
% inideces_kfold = crossvalind('KFold', size(kfold_set, 1), 10);
% 
% Models = {};
% accuracies = [];
% charts = {};
% measurements=["F1","precision","recall","accuracy"];
% 
%     F1 = [];
%     recalls= [];
%     precisions= [];
% 
% 
% 
% for s = 1:10
%     within_fold = find(inideces_kfold==s);
%     outside_fold = find(inideces_kfold~=s);
% 
% 
%     
%     SVMModel = perform_svm(randomized_X, randomized_Y, within_fold, outside_fold);
%     
%     holdout_predicts = predict(SVMModel, randomized_X(holdout_set, :));
%     [f1, precision, recall, accuracy] = printClassMetrics(randomized_Y(holdout_set, :), holdout_predicts)
%     
%     chart = confusionchart(randomized_Y(holdout_set, :), holdout_predicts)
%     % loss = crossentropy(randomized_Y(holdout_set, :), holdout_predicts)
%     accuracies(end+1) = accuracy;
%     F1(end+1) = f1;
%     recalls(end+1) = recall;
%     precisions(end+1) = precision;
% 
%    Models{end+1} = SVMModel;
% 
% 
%     
% 
% pause(1)
% end
% 
% % plot best confusion
% 
% 
% [M,I] = max(accuracies);
% 
% 
% 
% 
% ones = find(indices==1);
% zeros = find(indices==0); 
% %% ML
% 
% function [SVMModel, predictions, c] = perform_svm(randomized_X, randomized_Y, ones, zeros)
% 
% X_train = randomized_X(ones, :);
% Y_train = randomized_Y(ones, :);
% X_test = randomized_X(zeros, :);
% Y_test = randomized_Y(zeros, :);
% 
% 
% 
% SVMModel = fitcsvm(X_train,Y_train,'Standardize',true, 'KernelFunction','linear','ClassNames',{'burning','control'});
% 
% %sv = SVMModel.SupportVectors;
% %
% predictions = predict(SVMModel, X_test);
% 
% c = confusionchart(Y_test, predictions);
% 
% 
% 
% end
% 
% 
% 
% %figure
% %gscatter(X(:,:),X(:,:),Y)
% %hold on
% %plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
% %legend("exploding", "control", "Support Vector")

%%


%% Functions
function [feat_mat] = features(event, filtered_O2, fs, condition)
    % select event
    if strcmp(condition, 'exploding')
        box_epoch = epoching(0,0.66,250,event,filtered_O2, 0);
    elseif strcmp(condition, 'burning')
        [~,box_epoch,~] = epoching(0,0.66,250,event,filtered_O2, 0);
    else
        [~,~,box_epoch] = epoching(0,0.66,250,event,filtered_O2, 0);
    end
    % Variance with windowing
    feat_var = NaN(size(box_epoch,1), 4);
    feat_var(:,1) = var(box_epoch,0,2);
    windowed_box_epoch = epoching(0,0.22,250,event,filtered_O2, 0);
    windowed_box_epoch = windowed_box_epoch(1:size(box_epoch,1)); % take no more samples than matrix size
    feat_var(:,2) = var(windowed_box_epoch,0,2);
    windowed_box_epoch = epoching(0.22,0.44,250,event,filtered_O2, 0);
    windowed_box_epoch = windowed_box_epoch(1:size(box_epoch,1));
    feat_var(:,3) = var(windowed_box_epoch,0,2);
    windowed_box_epoch = epoching(0.44,0.66,250,event,filtered_O2, 0);
    windowed_box_epoch = windowed_box_epoch(1:size(box_epoch,1));
    feat_var(:,4) = var(windowed_box_epoch,0,2);
    
    % PSD
    [feat_psd, f_vec] = pwelch(box_epoch',[],[],[],fs);
    feat_psd = feat_psd';
    [~,idx_down] = min(abs(f_vec-1));
    [~,idx_up] = min(abs(f_vec-30));
    feat_psd = feat_psd(:,idx_down:idx_up);
    
    % DWT
    feat_dwt = [];
    for i = 1:size(box_epoch,1)
        feat_dwt = [feat_dwt; wavedec(box_epoch(i,:),3,'db8')];
    end
    
    % Feature vector
    feat_mat = [feat_var, feat_psd, feat_dwt];
    feat_mat = normalize(feat_mat);
end

function [exploding_box_epoch,burning_box_epoch,control_condition_epoch] = epoching(epoch_start,epoch_end,fs,event, filtered_O2, do_plot)
    exploding_box_epoch = [];
    control_condition_epoch = []; 
    burning_box_epoch = [];
    
    epoch_start_samples = round(epoch_start*fs);
    epoch_end_samples = round(epoch_end*fs); 
    
    for idx = 1:length(event)
        switch event(idx).values
            case 'S12' % Exploding Box Marker
                disp('S12')
                exploding_box_epoch = [exploding_box_epoch; filtered_O2(event(idx).samples+epoch_start_samples:event(idx).samples+epoch_end_samples)];
            case 'S13' % Control Condition Marker 
                %disp('S13')
                control_condition_epoch = [control_condition_epoch; filtered_O2(event(idx).samples+epoch_start_samples:event(idx).samples+epoch_end_samples)];
            case 'S14' % Burning Box Marker
                %disp('S14') 
                burning_box_epoch = [burning_box_epoch; filtered_O2(event(idx).samples+epoch_start_samples:event(idx).samples+epoch_end_samples)];
        end
    end
    
    if do_plot
    figure(10)
    for idx = 1:size(exploding_box_epoch,1)
        plot(epoch_start:1/fs:epoch_end, exploding_box_epoch(idx,:))
        hold on 
    end
    hold off
    title('Explosion Epoched trials')
    xlabel('Time (s)')
    ylabel('Amplitude (a.u)')
    saveas(gcf, 'Explosion_Epoched_trials.jpg')
    
    
    figure(11)
    for idx = 1:size(control_condition_epoch,1)
        plot(epoch_start:1/fs:epoch_end, control_condition_epoch(idx,:))
        hold on 
    end
    hold off
    title('NO Explosion Epoched trials')
    xlabel('Time (s)')
    ylabel('Amplitude (a.u)')
    saveas(gcf, 'NO_Explosion_Epoched_trials.jpg')
    end
end
