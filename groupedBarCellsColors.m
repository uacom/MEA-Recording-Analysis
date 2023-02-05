function [] = groupedBarCellsColors(data,varargin)
%% Same as groupedBar function, but with data in cells because of differing lengths
% INPUT: data = cell array, where each cell represents one column to graph

%%% load('C:\Users\Shenfeng Qiu\Documents\MATLAB\mycode\Graphing Scripts\graphing.mat')
%%% data = grouped2;
%%% groupedBarCellsColors(grouped2);
if nargin == 2
    if varargin{1} == 1
        noDots = 1;
    else
        noDots = 0;
    end
else
    noDots = 0;
end

% Creates bar graph with errorbars of data
numCol = length(data);
avg = zeros(1,numCol);
sem = zeros(1,numCol);
for i=1:numCol
    avg(i) = mean(data{i});
    sem(i) = buzsem(data{i});
end

corder = get(gca,'colororder');
for i=1:numCol
    switcher = mod(i,7);
    switch switcher
        case 1
            b = bar(i,avg(i),'EdgeColor',corder(switcher,:), 'FaceColor', 'none', 'LineWidth', 1.5);
            b.BarWidth = 0.6;

            hold on
            errorbar(i,avg(i),sem(i),'.','Color',corder(switcher,:),'LineWidth',1.5)
        case 2
            b = bar(i,avg(i),'EdgeColor',corder(switcher,:), 'FaceColor', 'none', 'LineWidth', 1.5);
            b.BarWidth = 0.6;

            hold on
            errorbar(i,avg(i),sem(i),'.','Color',corder(switcher,:),'LineWidth',1.5)
        case 3
            b = bar(i,avg(i),'EdgeColor',corder(switcher,:), 'FaceColor', 'none', 'LineWidth', 1.5);
            b.BarWidth = 0.6;

            hold on
            errorbar(i,avg(i),sem(i),'.','Color',corder(switcher,:),'LineWidth',1.5)
        case 4
            b = bar(i,avg(i),'EdgeColor',corder(switcher,:), 'FaceColor', 'none', 'LineWidth', 1.5);
            b.BarWidth = 0.6;

            hold on
            errorbar(i,avg(i),sem(i),'.','Color',corder(switcher,:),'LineWidth',1.5)
        case 5
            b = bar(i,avg(i),'EdgeColor',corder(switcher,:), 'FaceColor', 'none', 'LineWidth', 1.5);
            b.BarWidth = 0.6;

            hold on
            errorbar(i,avg(i),sem(i),'.','Color',corder(switcher,:),'LineWidth',1.5)
        case 6
            b = bar(i,avg(i),'EdgeColor',corder(switcher,:), 'FaceColor', 'none', 'LineWidth', 1.5);
            b.BarWidth = 0.6;

            hold on
            errorbar(i,avg(i),sem(i),'.','Color',corder(switcher,:),'LineWidth',1.5)
        case 0
            b = bar(i,avg(i),'EdgeColor',corder(7,:), 'FaceColor', 'none', 'LineWidth', 1.5);
            b.BarWidth = 0.6;

            hold on
            errorbar(i,avg(i),sem(i),'.','Color',corder(7,:),'LineWidth',1.5)
    end

end

numDots{numCol} = [];
for i=1:numCol
    numDots{i} = length(data{i});
end

if noDots == 0
    spacer = b.BarWidth/4;
    xScat{length(numDots)} = [];
    for i=1:length(numDots)
        for ii=1:numDots{i}
            xScat{i}(ii) = 1 + spacer;
            spacer = -spacer;
        end
    end
    for i=1:numCol
        switcher = mod(i,7);
        switch switcher
            case 1
                scatter(xScat{i}+(i-1),data{i},'o','filled','MarkerFaceColor',corder(switcher,:));
            case 2
                scatter(xScat{i}+(i-1),data{i},'o','filled','MarkerFaceColor',corder(switcher,:));
            case 3
                scatter(xScat{i}+(i-1),data{i},'o','filled','MarkerFaceColor',corder(switcher,:));
            case 4
                scatter(xScat{i}+(i-1),data{i},'o','filled','MarkerFaceColor',corder(switcher,:));
            case 5
                scatter(xScat{i}+(i-1),data{i},'o','filled','MarkerFaceColor',corder(switcher,:));
            case 6
                scatter(xScat{i}+(i-1),data{i},'o','filled','MarkerFaceColor',corder(switcher,:));
            case 0
                scatter(xScat{i}+(i-1),data{i},'o','filled','MarkerFaceColor',corder(7,:));
        end
    end
end
% for i=1:numCol
%     switch rem(i,7)
%         case 1
%             scatter(xScat{i}+(i-1),data{i},'o','filled','MarkerFaceColor','k');
%         case 2
%             scatter(xScat{i}+(i-1),data{i},'v','filled','MarkerFaceColor','k');
%         case 3
%             scatter(xScat{i}+(i-1),data{i},'s','filled','MarkerFaceColor','k');
%         case 4
%             scatter(xScat{i}+(i-1),data{i},'d','filled','MarkerFaceColor','k');
%         case 5
%             scatter(xScat{i}+(i-1),data{i},'^','filled','MarkerFaceColor','k');
%         case 6
%             scatter(xScat{i}+(i-1),data{i},'p','filled','MarkerFaceColor','k');
%         case 0
%             scatter(xScat{i}+(i-1),data{i},'h','filled','MarkerFaceColor','k');
%     end
% end
end