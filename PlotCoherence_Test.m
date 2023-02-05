function LINE = PlotCoherence_Test(COH1,COH2,Cohf,frang,Flag)
% 
% COH1 = squeeze(CauValue_S1(:,2,1,:));
% COH2 = squeeze(CauValue_S1(:,1,2,:));
% Cohf = Cauf_S1;
% frang = f_rang_cau;
% Flag = 0';

%
% plot the shadeerrorbar figure for Two Coherence Values and test the
% significant difference between them in all and specific frequency range
%
% input:
%      COH1, the first coherence/causality values with format: variables*fre
%      COH2, the second coherence/causality values with format: variables*fre
%      Cohf, the diftribution of frequency band
%      frang, specific frequency range
%      Flag, a flag for plotting COH or GC cpmparision,
%            Note: 1 for COH and 0 for GC
%

% check the format
if size(COH1,1) > size(COH1,2)
    COH1 = COH1';
end
if size(COH2,1) > size(COH2,2)
    COH2 = COH2';
end

COH_mean(1,:) = mean(COH1,1);
COH_mean(2,:) = mean(COH2,1);
COH_sem(1,:) = std(COH1,0,1)/sqrt(size(COH1,1));
COH_sem(2,:) = std(COH2,0,1)/sqrt(size(COH2,1));

% t test
for i = 1 : 1 : length(Cohf)
    [H(i),P(i)] = ttest2(COH1(:,i),COH2(:,i),0.05);
end

if Flag == 1
    LINE = ShadedErrorBar_Qiu(COH_mean,COH_sem,Cohf,'Coherence');
else
    LINE = ShadedErrorBar_Qiu(COH_mean,COH_sem,Cohf,'Causality');
end
xlim([0, length(COH_mean)/2]) % SQ add

sigVector = find(H==1);
ii = 1;
iii = 1;
j = 1;
i = 2;
ingSig = 0;
selVector{:} = 0;
while 1
    if i >= length(sigVector)
        break;
    else
    if sigVector(i)-sigVector(i-1) == 1
        while 1
            if sigVector(i)-sigVector(i-1) == 1
                selVector{ii}(iii+1) = Cohf(sigVector(i));
                selVector{ii}(iii) = Cohf(sigVector(i-1));
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

%yloc = 0.9*min(min([COH_mean-COH_sem]));
yloc = 1.1*max(max([COH_mean+COH_sem]));

hold on

for i = 1:length(selVector)
    plot([selVector{i}(1) (selVector{i}(length(selVector{i})))], [yloc yloc], 'k', 'LineWidth', 2)
    text(((selVector{i}(1) + selVector{i}(length(selVector{i})))/2),yloc+0.01,'*','FontName','Times New Roman','FontSize',13,'FontWeight','bold')
end

% for i = 1 : 1 : length(Cohf)
%     if H(i)==1
%         plot(Cohf(i),yloc,'b*')
%         
%     end
% end

floc = find(Cohf>=frang(1) & Cohf<=frang(2));
COH1_S = mean(COH1(:,floc),2);
COH2_S = mean(COH2(:,floc),2);
[H_f,~] = ttest2(COH1_S,COH2_S,0.05);

COH_S_mean(:,1) = mean(COH1_S);
COH_S_mean(:,2) = mean(COH2_S);
COH_S_sem(1) = std(COH1_S)/sqrt(length(COH1_S));
COH_S_sem(2) = std(COH2_S)/sqrt(length(COH2_S));

figure
set(gcf,'position',[200 200 400 400],'color','w')
handles.bar = bar(COH_S_mean,'edgecolor','w', 'linewidth',0.5);
handles.bar.FaceColor = [1 0 0];

hold on
% plot errorbar
for col = 1 : 1: length(COH_S_mean)
    if col ==1
        line([1 1],[COH_S_mean(1),COH_S_mean(1)+COH_S_sem(1)],'color',[1 0 0],'linewidth',2);
        line([0.8 1.2],[COH_S_mean(1)+COH_S_sem(1),COH_S_mean(1)+COH_S_sem(1)],'color',[1 0 0],'linewidth',2);
    else
        line([2 2],[COH_S_mean(2),COH_S_mean(2)+COH_S_sem(2)],'color',[1 0 0],'linewidth',2);
        line([1.8 2.2],[COH_S_mean(2)+COH_S_sem(2),COH_S_mean(2)+COH_S_sem(2)],'color',[1 0 0],'linewidth',2);
    end
end
hold on
if Flag == 1
    ylabel('Coherence','FontName','Times New Roman','FontSize',15,'FontWeight','bold')
else
    ylabel('Causality','FontName','Times New Roman','FontSize',15,'FontWeight','bold')
end
xlim([0.5 2+0.5]);
set(gca,'xticklabel',[])
set(gca,'xticklabel',{'State1';'State2';},'FontName','Times New Roman','FontSize',15,'FontWeight','bold')

% PLOT sig
line([1 1],[COH_S_mean(1)+COH_S_sem(1),1.1*max([COH_S_mean(1)+COH_S_sem(1) COH_S_mean(2)+COH_S_sem(2)])],'color','k','linewidth',3)
line([2 2],[COH_S_mean(2)+COH_S_sem(2),1.1*max([COH_S_mean(1)+COH_S_sem(1) COH_S_mean(2)+COH_S_sem(2)])],'color','k','linewidth',3)
line([1 2],[1.1*max([COH_S_mean(1)+COH_S_sem(1) COH_S_mean(2)+COH_S_sem(2)]),1.1*max([COH_S_mean(1)+COH_S_sem(1) COH_S_mean(2)+COH_S_sem(2)])],'color','k','linewidth',3)
if H_f==1
    text(1.45,1.125*max([COH_S_mean(1)+COH_S_sem(1) COH_S_mean(2)+COH_S_sem(2)]),'*','FontName','Times New Roman','FontSize',13,'FontWeight','bold')
else
    text(1.45,1.15*max([COH_S_mean(1)+COH_S_sem(1) COH_S_mean(2)+COH_S_sem(2)]),'NS','FontName','Times New Roman','FontSize',13,'FontWeight','bold')
end
