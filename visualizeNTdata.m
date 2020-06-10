%%
%Further Visualizations of data from '01 ANALYZE NEURAL TRACES'
%change i to number of wells/2 minus 1
%so for 8 wells, 2 replicates 
% num=2 and n=3
%n = number of treatments
%num = number of replicates
% 9 wells
num = 9;
%Eric June 1 Imaging 54 wells
%6 treatments
n = 6;

hmap2 = hmpeaks([1:num*n],:)';
hmap2std = hmpeakssd([1:num*n],:)';
hmap2var = hmpeaksvar([1:num*n],:)';

for i = 0:n-1
    hmap2ave(:,i+1) = nanmean(hmap2(:,i*num+1:(i+1)*num),2); 
    %hmap2std2(:,i+1) = nanstd(hmap2(:,i*num+1:(i+1)*num),0,2);

    %hmap2avestd(:,i+1) = nanmean(hmap2std(:,i*num+1:(i+1)*num),2);
    hmap2stdvar(:,i+1) = sqrt(nanmean(hmap2var(:,i*num+1:(i+1)*num),2)); 

    meanan(:,i+1) = nanmean(anDistm(:,i*num+1:(i+1)*num),2);
end

usethe = hmap2stdvar;

hmap2avenorm = hmap2ave ./ repmat(hmap2ave(1,:),length(hmap2ave),1);

meanan2 = repmat(meanan,length(usethe),1);
hmap2sem = usethe ./ sqrt(meanan2);

figure(1);
imagesc(hmap2);
box off

figure(2);
imagesc(flipud(hmap2ave'));
box off

figure(3);
p = [1 2 3 4 5 6 7 8 9 10 11 12 24 36 48];
colormap winter
cmap = colormap;
c = 4;

errorneg = (hmap2ave-hmap2sem)';
errorpos = (hmap2ave+hmap2sem)';
% errorbar(hmap2ave,hmap2std)
% jbfill(1:length(hmap2ave(:,1)),(hmap2ave-hmap2std)',(hmap2ave+hmap2std)',[189 189 189]/255,[189 189 189]/255,0,.25);
% getget = [1 2];
for i = 1:n %1 to n or % then change all getget(i) to just i
%for i = 1:length(getget)
    jbfill(1:length(hmap2ave(:,i)),errorneg(i,:),errorpos(i,:),coloring(i,:),cmap(i*c,:),0,.25); %fig5
    hold on;
    plot(hmap2ave(:,i),'color',cmap(i*(c),:),'LineWidth',1.5);
    xlim([0,length(hmap2ave(:,1))+1])
    hold on;
    box off
end

hold on;
for i = 1:length(p)
    area([(p(i)-0.5) (p(i)+0.5)], [1.5 1.5], 1.025, 'FaceColor',[227 26 28]/255,'EdgeColor','none');
    alpha(0.1);
    hold on;
end

hold on;
for i = 1:length(p)
    area([(p(i)-0.5) (p(i)+0.5)], [1.025 1.025], 1, 'FaceColor',cmap(i*c,:),'EdgeColor','none');
    hold on;
end

ylim([1,1.2])
%%
figure(4);
for j = 1:n*num 
    subplotv(num,n,j);
    area([25 75], [230 230], 0.9, 'FaceColor',[227 26 28]/255,'EdgeColor','none');
    alpha(0.1);
    hold on;
    
    for i = 1:length(p)%max(unique(max(AllCycler)))
        plot(nanmean(AllSqIntNorm(:,find(AllCycler(1,:) == p(i) & AllWellr(1,:) == j)),2),'color',cmap(i*c,:))
        ylim([1,1.5])
        xlim([0,150])
        hold on;
        box off
        coloring(i,:) = cmap(i*c,:);
        axis off
        %xticklabels({'0','5','10','15'})
    end
end

set(gcf, 'Position', [100, 500, 1750, 600])
%%
figure(5);
replicate = num;

for j = 0:n-1
    subplot(1,n,j+1);
    area([25 75], [230 230], 0.9, 'FaceColor',[227 26 28]/255,'EdgeColor','none');
    alpha(0.1);
    hold on;
    
    for i = 1:length(p)%max(unique(max(AllCycler)))
        one(:,1) = nanmean(AllSqIntNorm(:,find(AllCycler(1,:) == p(i) & AllWellr(1,:) == (j*replicate)+1)),2);
        two(:,1) = nanmean(AllSqIntNorm(:,find(AllCycler(1,:) == p(i) & AllWellr(1,:) == (j*replicate)+2)),2);
        three(:,1) = nanmean(AllSqIntNorm(:,find(AllCycler(1,:) == p(i) & AllWellr(1,:) == (j*replicate)+3)),2);
        four(:,1) = nanmean([one two three],2);
        plot(four,'color',cmap(i*c,:))
        ylim([1,1.4])
        xlim([0,150])
        hold on;
        box off
        coloring(i,:) = cmap(i*c,:);
        %axis off
    end
end

set(gcf, 'Position', [100, 100, 1750, 125])

figure(6);
rgbplot(coloring)
colormap(coloring)
colorbar('Ticks',[])

%%
figure(8);
usetimepoint = 48;
%normhm2ave = hmap2ave(usetimepoint,:) - hmap2ave(usetimepoint,1);

bar(hmap2ave(usetimepoint,:)); %change : to x:n for specific wells
hold on;
errorbar(hmap2ave(usetimepoint,:),hmap2sem(usetimepoint,:),'.','color','k') %change : to x:n for specific wells
ylim([1,1.15])

%%
%for 3 replicates (num) (may need to add or subtract well# if more or less)
for i = 0:n-1
    well1 = length(AllWellr(find(AllCycler(1,:) == 1 & AllWellr(1,:) == i*num+1)));
    well2 = length(AllWellr(find(AllCycler(1,:) == 1 & AllWellr(1,:) == i*num+2)));
    well3 = length(AllWellr(find(AllCycler(1,:) == 1 & AllWellr(1,:) == i*num+3)));
    getthem(i+1,:) = well1+well2+well3;
end

%for 11 groups (n) (may need to add or subtract repelem if more or less)
gWells = horzcat(repelem(1,getthem(1)),...
    repelem(2,getthem(2)),...
    repelem(3,getthem(3)),...
    repelem(4,getthem(4)),...
    repelem(5,getthem(5)),...
    repelem(6,getthem(6)),...
    repelem(7,getthem(7)),...
    repelem(8,getthem(8)));

%stats for usetimepoint cycle
[p1,tbl1,stats1] = anova1(peaks(find(AllCycler(1,:) == usetimepoint)),gWells);
checkStats1 = multcompare(stats1,'CType','bonferroni','Alpha',0.005);

% [p1,tbl1,stats1] = anova1(peaks(find(AllCycler(1,:) == usetimepoint & AllWellr(1,:) <= 12)),gWells);
% checkStats1 = multcompare(stats1,'CType','bonferroni','Alpha',0.005);


