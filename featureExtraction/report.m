
clc;
clear;
close all;
load('task_Day3to5.mat')
% defect electrode PO7
dataEEG.label(7)=[];
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
%% std plotting 
% Explosion vs. Control

fig1=figure;

x = 1:376;
x2 = [x, fliplr(x)];
%average/mean for explosion
curve1 = ETM_O2 + ET_std_O2;
curve2 = ETM_O2 - ET_std_O2;

inBetween = [curve1 fliplr(curve2)];

fill(x2, inBetween, 'r', 'facealpha', 0.3, 'edgecolor', 'none');
hold on
plot(ETM_O2, 'r');



%Average for control
curve1 = CTM_O2 + CT_std_O2;
curve2 = CTM_O2 - CT_std_O2;
inBetween = [curve1 fliplr(curve2)];
fill(x2, inBetween, 'b', 'facealpha', 0.3, 'edgecolor', 'none');
plot(CTM_O2, 'b');
axis ([0 376 -20 20]);
hold off
xlabel('Time (s)');
ylabel('Magintude (u.a)');
set(gca,'FontSize', 10);
title('Average visualization for explosion trials VS control trials on Channel O2');
legend('std of exploding box event','mean of exploding box','std of control',' mean of control');
set(gcf, 'PaperPositionMode','auto','Units','Centimeters','Position',[2 2 20 15],'PaperSize', [20 15]);
savefig(fig1);
print(fig1, '-dpdf','MEAN Explosion trials VS control trials');
%%%%%%% Burning vs. Control

fig2 = figure;
%Average for burning trials
curve1 = BTM_O2 + BT_std_O2;
curve2 = BTM_O2 - BT_std_O2;
inBetween = [curve1 fliplr(curve2)];
fill(x2, inBetween, 'r', 'facealpha', 0.3, 'edgecolor', 'none');
hold on
plot(BTM_O2, 'r');
curve1 = CTM_O2 + CT_std_O2;
curve2 = CTM_O2 - CT_std_O2;
inBetween = [curve1 fliplr(curve2)];
fill(x2, inBetween, 'b', 'facealpha', 0.3, 'edgecolor', 'none');
plot(CTM_O2, 'b');
axis ([0 376 -15 20]);
hold off
xlabel('Time (s)');
ylabel('Magintude (u.a)');
set(gca,'FontSize', 10);
title('Average visualization for burning trials VS control trials on Channel O2');
legend('std of burning box event','mean of burning box','std of control',' mean of control');
set(gcf, 'PaperPositionMode','auto','Units','Centimeters','Position',[2 2 20 15],'PaperSize', [20 15]);
savefig(fig2);
print(fig2, '-dpdf','MEAN burning trials VS control trials');

%% Variance calculation for all channels
%Variance over all channels

%Burning trials variance across all rows(trials)
VAR_BT = var(BT_WS,'',2); 
VAR_BT=squeeze(VAR_BT);
VAR_BT_MEAN=mean(VAR_BT(:,:),2);
VAR_BT_FEAT=(mean(VAR_BT(:,:),1))';

%Control trials variance across all rows(trials)
VAR_CT = var(CT_WS,'',2);
VAR_CT=squeeze(VAR_CT);
VAR_CT_MEAN=mean(VAR_CT(:,:),2);
VAR_CT_FEAT = (mean(VAR_CT(:,:),1))';

%Explosion event variance across all rows(trials) 
VAR_expl = var(ET_WS,'',2);
VAR_expl= squeeze(VAR_expl);
VAR_ET_MEAN=mean(VAR_expl(:,:),2);
VAR_ET_FEAT=(mean(VAR_expl(:,:),1))';

%% Variance O2 Channel


%Variance O2 Channel
fig3=figure;
subplot(2,2,1)
a=scatter((1:11),VAR_CT,'g','filled');
xticks(0:1:11);
xlabel('Channels');
ylabel('Variance');
grid on;
hold on
title('Control Trial VS Burning Trial ');
b=scatter((1:11),VAR_BT,'black','filled');
legend([a(1),b(1)],'Control trials','Burning trials');
hold off
grid off

subplot(2,2,2)
c=scatter((1:11),VAR_CT_MEAN,'r','filled');
xticks(0:1:11);
xlabel('Channels');
ylabel('Variance');
grid on;
hold on
title('MEAN: Control Trial VS Burning Trial ');
d=scatter((1:11),VAR_BT_MEAN,'b','filled');
legend([c(1),d(1)],'Control trials Mean','Burning trials Mean');
hold off
grid off



subplot(2,2,3)
a=scatter((1:11),VAR_CT,'g','filled');
xticks(0:1:11);
xlabel('Channels');
ylabel('Variance');
grid on;
hold on
title('Control Trial VS explosion Trial ')
b=scatter((1:11),VAR_expl,'black','filled');
legend([a(1),b(1)],'Control trials','Explosion trials');
hold off
grid off

subplot(2,2,4)
c=scatter((1:11),VAR_CT_MEAN,'r','filled');
xticks(0:1:11);
xlabel('Channels');
ylabel('Variance');
grid on;
hold on
title('MEAN: Control Trial VS Burning Trial ');
d=scatter((1:11),VAR_ET_MEAN,'b','filled');
legend([c(1),d(1)],'Control trials Mean','Burning trials Mean');
hold off
grid off

set(gcf, 'PaperPositionMode','auto','Units','Centimeters','Position',[2 2 20 15],'PaperSize', [20 15]);
savefig(fig3);
print(fig3, '-dpdf','Variance Interpretation over all channels');




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
    subplot(2,1,1)
    p1=scatter(fB,pow2db(PSDB_MEAN),'g','filled');
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    grid on
    hold on
    title('PSD-Burning trials VS PSD-Control trials')
    p2=scatter(fC,pow2db(PSDC_MEAN),'black','filled');
    legend([p1(1),p2(1)],'Burning trials','Control trials');

    subplot(2,1,2)
    p2=scatter(fC,pow2db(PSDC_MEAN),'black','filled');
    xlabel('Frequency (Hz)')
    ylabel('Power (dB)')
    grid on
    hold on
    title('PSD-Explosion trials VS PSD-Control trials')
    p3=scatter(fE,pow2db(PSDE_MEAN),'g','filled');
    legend([p2(1),p3(1)],'Control trials','Explosion trials');
    set(gcf, 'PaperPositionMode','auto','Units','Centimeters','Position',[2 2 20 15],'PaperSize', [20 15]);
    savefig(fig4);
    print(fig4, '-dpdf','PSD Interpretation');

    
     


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

fig5 = figure;
% DWT - Plot of the Approximation coefficient during trial 1 for channel O2
% for the 3 different types of trials 

subplot(2,1,1)
plot(1:60,A1(3,:,1));
hold on
ylabel('Approximate coefficients (u.a)');
xlabel('Coefficient Index');
grid on;
plot(1:60,A2(3,:,1));
title('Trial 1: DWT  over O2 channels');
legend('Burning trials', 'Control trials');
axis([1 60 -80 80]);
set(gca,'FontSize', 10);

subplot(2,1,2)
plot(1:60,A3(3,:,1));
hold on
xlabel('Coefficient Index');
ylabel('Approximate coefficients (u.a)');
grid on;
plot(1:60,A2(3,:,1));
hold off;
legend('Explosion trials','Control trials');
axis([1 60 -80 80]);
set(gca,'FontSize', 10);

set(gcf, 'PaperPositionMode','auto','Units','Centimeters','Position',[2 2 20 15],'PaperSize', [20 15]);
savefig(fig5);
print(fig5, '-dpdf','DWT Interpretation');



%% Machine Learning SVM


explode_feat = [VAR_expl, PSDE, DWTE]
burn_feat
control_feat

control_feat_mat = features(event, filtered_O2, 250, 'control');
burning_feat_mat = features(event, filtered_O2, 250, 'burning'); 

X = [burning_feat_mat; control_feat_mat];
X(:,2) = [];
Y = cell(size(X,1), 1);
Y(1:size(burning_feat_mat, 1)) = {'burning'};
Y(size(burning_feat_mat, 1) + 1:end) = {'control'};


%% Shuffle Data

permutations = randperm(size(X, 1));

randomized_X = X(permutations, :);
randomized_Y = Y(permutations, :);

%% Cross Validation
indices = crossvalind('HoldOut', size(randomized_X, 1), 0.2);
ones = find(indices==1);
zeros = find(indices==0); 
X_train = randomized_X(ones, :);

Y_train = randomized_Y(ones, :);
X_test = randomized_X(zeros, :);
Y_test = randomized_Y(zeros, :);
%% ML

SVMModel = fitcsvm(X_train,Y_train,'Standardize',true, 'KernelFunction','linear','ClassNames',{'exploding','control'});

%sv = SVMModel.SupportVectors;
%
predictions = predict(SVMModel, X_test);

c = confusionchart(Y_test, predictions)

%figure
%gscatter(X(:,:),X(:,:),Y)
%hold on
%plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
%legend("exploding", "control", "Support Vector")


