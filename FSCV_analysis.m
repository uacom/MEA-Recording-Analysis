function [IT1,spectrogram1,IT2,spectrogram2] = FSCV_analysis(filepath)

% 95% CI with shaded area for dopamine recording
% FOR DEBUGGING filepath = 'NAcc data.xlsx';
%%% filepath = 'NAcc data.xlsx';  sheet1 = five IT traces, sheet2 =
%%% colormap data; 
%%% filepath = 'NAcc data1.xlsx';
%%% [IT1,spectrogram1,IT2,spectrogram2] = FSCV_analysis(filepath);
IT1 = xlsread(filepath, 'Sheet1'); 
spectrogram1 = xlsread(filepath, 'Sheet2');
IT2 = xlsread(filepath, 'Sheet3');
spectrogram2 = xlsread(filepath, 'Sheet4');

%% take excel data and store in meaningful format

for i=1:150
    meanIT1(i) = mean(IT1(i,:)); %#ok<*AGROW>
    meanIT2(i) = mean(IT2(i,:));
end

for i=1:length(IT1(1,:))
    sortedColormap1{i} = spectrogram1((1000*i - 999):(1000*i),:);
    sortedColormap2{i} = spectrogram2((1000*i - 999):(1000*i),:);
end

%% average the colormaps together
if length(IT1(1,:)) ~= 1
    meanColormap1 = zeros(1000,150);
    meanColormap2 = zeros(1000,150);
    for i=2:length(IT1(1,:))
        meanColormap1 = meanColormap1 + sortedColormap1{i-1} + sortedColormap1{i};
        meanColormap2 = meanColormap2 + sortedColormap2{i-1} + sortedColormap1{i};
    end

    meanColormap1 = meanColormap1 ./ length(IT1(1,:));
    meanColormap2 = meanColormap2 ./ length(IT2(1,:));
else
    meanColormap1 = sortedColormap1{1};
    meanColormap2 = sortedColormap2{1};
end

%% calculate 95% CI's
for i=1:150
    SEM1(i) = std(IT1(i,:)/sqrt(length(IT1(i,:))));
    ts = tinv([0.025 0.975],length(IT1(i,:))-1);
    CIpos1(i) = mean(IT1(i,:)) + ts(2)*SEM1(i); %#ok<*NASGU>
    CIneg1(i) = mean(IT1(i,:)) + ts(1)*SEM1(i);
end
for i=1:150
    SEM2(i) = std(IT1(i,:)/sqrt(length(IT2(i,:))));
    ts = tinv([0.025 0.975],length(IT2(i,:))-1);
    CIpos2(i) = mean(IT2(i,:)) + ts(2)*SEM2(i);
    CIneg2(i) = mean(IT2(i,:)) + ts(1)*SEM2(i);
end

t = 0.1:0.1:15;

load('demonColormap')
[cmap] = demonColormap; % Create colormap


%% plot data from first condition

figure(1)
subplot (2,1,1)
t2 = [t,fliplr(t)];
inBetweenPos = [meanIT1, fliplr(meanIT1+SEM1)];
fill(t2,inBetweenPos,[1 0.7 0.7],'EdgeColor','none');
hold on;
inBetweenNeg = [meanIT1-SEM1, fliplr(meanIT1)];
fill(t2,inBetweenNeg,[1 0.7 0.7],'EdgeColor','none');
hold on

plot(t,meanIT1+SEM1,'color',[1 0.7 0.7])
hold on
plot(t,meanIT1-SEM1,'color',[1 0.7 0.7])
hold on

plot(t,meanIT1,'r','LineWidth',2)
hold on
% plot(t,CIpos1,'color',[1 0.5 0.5])
% hold on
% plot(t,CIneg1,'color',[1 0.5 0.5])
hold off
xlabel('Time (s)');
ylabel('Current (nA)');
title('Dopamine Redox - First Condition');

subplot(2,1,2)
imagesc(meanColormap1);
set(gca,'Ydir','Normal');
xlabel('Time (s)')
ylabel('Point')
title('Colormap')
%[cmap] = buildcmap('bkymcg');
colormap(cmap)
xticks([0 50 100 150])
xticklabels({'0','5','10','15'})
colorbar

%% plot data from second condition

figure(2)
subplot (2,1,1)
t2 = [t,fliplr(t)];
inBetweenPos = [meanIT2, fliplr(meanIT2+SEM2)];
fill(t2,inBetweenPos,[0.7 0.7 1],'EdgeColor','none');
hold on;
inBetweenNeg = [meanIT2-SEM2, fliplr(meanIT2)];
fill(t2,inBetweenNeg,[0.7 0.7 1],'EdgeColor','none');
hold on

plot(t,meanIT2+SEM2,'color',[0.7 0.7 1])
hold on
plot(t,meanIT2-SEM2,'color',[0.7 0.7 1])
hold on

plot(t,meanIT2,'b','LineWidth',2)
hold on
% plot(t,CIpos2,'color',[0.5 0.5 1])
% hold on
% plot(t,CIneg2,'color',[0.5 0.5 1])
hold off
xlabel('Time (s)');
ylabel('Current (nA)');
title('Dopamine Redox - Second Condition');

subplot(2,1,2)
imagesc(meanColormap2);
set(gca,'Ydir','Normal');
xlabel('Time (s)')
ylabel('Point')
title('Colormap')
%[cmap] = buildcmap('bkymcg');
colormap(cmap)
xticks([0 50 100 150])
xticklabels({'0','5','10','15'})
colorbar

%% merge the two dopamine current plots and perform t-test, determine which intervals are significant

figure(3)

% plot first condition
t2 = [t,fliplr(t)];
inBetweenPos = [meanIT1, fliplr(meanIT1+SEM1)];
fill(t2,inBetweenPos,[1 0.7 0.7],'EdgeColor','none');
hold on;
inBetweenNeg = [meanIT1-SEM1, fliplr(meanIT1)];
fill(t2,inBetweenNeg,[1 0.7 0.7],'EdgeColor','none');
hold on

plot(t,meanIT1+SEM1,'color',[1 0.7 0.7])
hold on
plot(t,meanIT1-SEM1,'color',[1 0.7 0.7])
hold on

plot(t,meanIT1,'r','LineWidth',2)
% hold on
% plot(t,CIpos1,'color',[1 0.5 0.5])
% hold on
% plot(t,CIneg1,'color',[1 0.5 0.5])
hold off
xlabel('Time (s)');
ylabel('Current (nA)');
title('Dopmaine Redox - both conditions');
hold on

% plot second condition on top of first condition in different color
t2 = [t,fliplr(t)];
inBetweenPos = [meanIT2, fliplr(meanIT2+SEM2)];
fill(t2,inBetweenPos,[0.7 0.7 1],'EdgeColor','none');
hold on;
inBetweenNeg = [meanIT2-SEM2, fliplr(meanIT2)];
fill(t2,inBetweenNeg,[0.7 0.7 1],'EdgeColor','none');
hold on

plot(t,meanIT2+SEM2,'color',[0.7 0.7 1])
hold on
plot(t,meanIT2-SEM2,'color',[0.7 0.7 1])
hold on

plot(t,meanIT2,'b','LineWidth',2)
% hold on
% plot(t,CIpos2,'color',[0.5 0.5 1])
% hold on
% plot(t,CIneg2,'color',[0.5 0.5 1])
% hold off

% t test
H = zeros(150,1);
p = zeros(150,1);
for i=1:150
    [H(i), p(i)] = ttest2(IT1(i,:), IT2(i,:));
end

%% script to add bar over significant data

sigVector = find(H==1);
ii = 1;
iii = 1;
j = 1;
i = 2;
ingSig = 0;
selVector{:} = 0;
while 1
    if i > length(sigVector)
        break;
    else
    if sigVector(i)-sigVector(i-1) == 1
        while 1
            if sigVector(i)-sigVector(i-1) == 1
                selVector{ii}(iii+1) = sigVector(i)/10;
                selVector{ii}(iii) = sigVector(i-1)/10;
            else
                ii = ii+1;
                iii = 1;
                break;
            end
            if i+1 < length(sigVector)
                i = i+1;
            else
                break;
            end
            iii = iii+1;
        end
        i = i+1;
    else
        indSig(j) = sigVector(i);
        j = j+1;
        i = i+1;
    end
    end
end

yloc = 1.1*max(max((meanIT1+SEM1+1)));

hold on

if isempty(selVector) == 0
    for i = 1:length(selVector)
        plot([selVector{i}(1) (selVector{i}(length(selVector{i})))], [yloc yloc], 'k', 'LineWidth', 2)
        text(((selVector{i}(1) + selVector{i}(length(selVector{i})))/2),yloc+0.1,'*','FontName','Times New Roman','FontSize',13,'FontWeight','bold')
    end
else
    if isempty(sigVector)
        text(10,5,'NS','FontName','Times New Roman','FontSize',13,'FontWeight','bold')
    end
    
    if length(sigVector) == 1
        ylim([0 yloc+1])
        text(sigVector/10,yloc+0.1,'*','FontName','Times New Roman','FontSize',13,'FontWeight','bold')
    end
end





end