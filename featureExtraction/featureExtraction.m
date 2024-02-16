
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
figSize = [2 2 17 10];
set(fig1, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
subplot(1,2,1)
x = ((1:376)/250)-0.5;
x2 = [x, fliplr(x)];
%average/mean for explosion
curve1 = ETM_O2 + ET_std_O2;
curve2 = ETM_O2 - ET_std_O2;

inBetween = [curve1 fliplr(curve2)];

fill(x2, inBetween, 'r', 'facealpha', 0.3, 'edgecolor', 'none');
hold on

plot(x,ETM_O2, 'r');



%Average for control
curve1 = CTM_O2 + CT_std_O2;
curve2 = CTM_O2 - CT_std_O2;
inBetween = [curve1 fliplr(curve2)];
fill(x2, inBetween, 'b', 'facealpha', 0.3, 'edgecolor', 'none');

plot(x,CTM_O2, 'b');
axis ([-0.4960 1 -20 20]);
text(-0.196,-1.81388,'\leftarrow N1(-0.196,-2.90665)',fontsize=8)
text(0.18,8.99116,'\leftarrow P1(0.18,8.99116)',fontsize=8)
text(0.312,-7.38247,'\leftarrow N2(0.312,-7.38247)',fontsize=8)
hold off
xlabel('Time (s)');
ylabel('Magintude (u.a)');
% xL = xlim;
% yL = ylim;
% line([0 0], yL);  %x-axis
% line(xL, [0 0]);  %y-axis
set(gca,'FontSize', 10);
title('Average on Explosion VS Control Trials-Channel O2',fontsize=8);
legend('std of exploding box event','mean of exploding box','std of control',' mean of control');
set(gcf, 'PaperPositionMode','auto','Units','Centimeters','Position',[2 2 20 15],'PaperSize', [20 15]);

%%%%%%% Burning vs. Control


subplot(1,2,2)
%Average for burning trials
curve1 = BTM_O2 + BT_std_O2;
curve2 = BTM_O2 - BT_std_O2;
inBetween = [curve1 fliplr(curve2)];
fill(x2, inBetween, 'r', 'facealpha', 0.3, 'edgecolor', 'none');
hold on
plot(x,BTM_O2, 'r');
curve1 = CTM_O2 + CT_std_O2;
curve2 = CTM_O2 - CT_std_O2;
inBetween = [curve1 fliplr(curve2)];
fill(x2, inBetween, 'b', 'facealpha', 0.3, 'edgecolor', 'none');
plot(x,CTM_O2, 'b');
axis ([-0.4960 1 -20 20]);
text(-0.024,5.21834,'\leftarrow P1(-0.024,5.21834)',fontsize=8)
text(0.308,-10.2401,'\leftarrow N1(0.308,-10.2401)',fontsize=8)
text(0.512,3.45559,'\leftarrow P2(0.512,3.45559)',fontsize=8)
hold off
xlabel('Time (s)');
ylabel('Magintude (u.a)');

% xL = xlim;
% yL = ylim;
% line([0 0], yL);  %x-axis
% line(xL, [0 0]);  %y-axis
set(gca,'FontSize', 10);
title('Average on Burning VS Control Trials-Channel O2',fontsize=8);
legend('std of burning box event','mean of burning box','std of control',' mean of control');
set(gcf, 'PaperPositionMode','auto','Units','Centimeters','Position',[2 2 25 10],'PaperSize', [25 10]);
savefig(fig1);
print(fig1, '-dpdf','MEAN');

%% Variance calculation for all channels
%Variance over all channels

%Burning trials variance across all rows(trials)
VAR_BT = var(BT_WS,'',2); 
VAR_BT=squeeze(VAR_BT);
VAR_BT_MEAN=mean(VAR_BT(:,:),2);
%Control trials variance across all rows(trials)
VAR_CT = var(CT_WS,'',2);
VAR_CT=squeeze(VAR_CT);
VAR_CT_MEAN=mean(VAR_CT(:,:),2);
%Explosion event variance across all rows(trials) 
VAR_expl = var(ET_WS,'',2);
VAR_expl= squeeze(VAR_expl);
VAR_ET_MEAN=mean(VAR_expl(:,:),2);

%% Variance O2 Channel


%Variance O2 Channel
figSize = [2 2 18 20];
fig3 = figure;
set(fig3, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);
subplot(2,2,1)
a=scatter((1:11),VAR_CT,'g','filled');
xticks(0:1:11);
xlabel('Channels');
ylabel('Variance');
grid on;
hold on
title('VAR: Control Trial VS Burning Trial ',fontsize=8);
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
title('VAR-MEAN: Control Trial VS Burning Trial ',fontsize=8);
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
title('VAR: Control Trial VS explosion Trial',fontsize=8)
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
title('VAR-MEAN: Control Trial VS Explosion Trial',fontsize=8);
d=scatter((1:11),VAR_ET_MEAN,'b','filled');
legend([c(1),d(1)],'Control trials Mean','Explosion trials Mean');
hold off
grid off

set(gcf, 'PaperPositionMode','auto','Units','Centimeters','Position',[2 2 20 15],'PaperSize', [20 15]);
savefig(fig3);
print(fig3, '-dpdf','VAR');




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
A12=mean(A1,2);
A12=squeeze(A12);
A13=mean(A12,2);

A21=mean(A2,2);
A21=squeeze(A21);
A23=mean(A21,2);

A31=mean(A3,2);
A31=squeeze(A31);
A33=mean(A31,2);

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

%% ff
fig6 = figure;
% DWT - Plot of the Approximation coefficient during trial 1 for channel O2
% for the 3 different types of trials 


plot(1:60,A1(3,:,1));
hold on
ylabel('Approximate coefficients (u.a)');
xlabel('Coefficient Index');
grid on;
plot(1:60,A2(3,:,1));
plot(1:60,A3(3,:,1));
title('Trial 1: DWT  over O2 channels',fontsize=8);
legend('Burning', 'Control', 'Explosion',fontsize=5);
axis([1 60 -80 80]);
set(gcf, 'PaperPositionMode','auto','Units','Centimeters','Position',[2 2 20 15],'PaperSize', [20 15]);
savefig(fig6);
print(fig6, '-dpdf','DWT');

%% gg
figSize = [2 2 20 15];
fig7 = figure;
set(fig7, 'Units', 'Centimeters', 'PaperPositionMode', 'Auto',...
    'PaperUnits', 'Centimeters','Position', figSize,...
    'PaperSize', [figSize(3) figSize(4)]);

subplot(2,1,1)

scatter(1:11,A13);
hold on
scatter((1:11),A23);
scatter((1:11),A33);
legend('Burning','Control','Explosion')
xlabel('Working Channels');
ylabel('Approximate coefficients (u.a)');
title('Average Approximate coefficients over all trials over all Channels',FontSize=10)

xticks(0:1:11);
grid on;
subplot(2,1,2)
plot(1:60,A1(3,:,1));
hold on
ylabel('Approximate coefficients (u.a)');
xlabel('Coefficient Index');
grid on;
plot(1:60,A2(3,:,1));
plot(1:60,A3(3,:,1));
title('Trial 1: App Coeffcient  over O2 channels',fontsize=10);
legend('Burning', 'Control', 'Explosion',fontsize=5);
axis([1 60 -80 80]);
set(gcf, 'PaperPositionMode','auto','Units','Centimeters','Position',[2 2 20 15],'PaperSize', [20 15]);




savefig(fig7);
print(fig7, '-dpdf','MEAN DWT');

