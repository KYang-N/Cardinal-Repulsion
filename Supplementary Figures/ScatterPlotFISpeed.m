%% Scatter plot of minimum FI against maximum absolute drift speed

clear
DataDir = '';
load([DataDir,'TwoLayerMinFI.mat']);
TwoLayerFI = reshape(FIMatrix,1,numel(FIMatrix));
load([DataDir,'OneLayerMinFI.mat']);
OneLayerFI = reshape(FIMatrix,1,numel(FIMatrix));
load([DataDir,'TwoLayerDriftSpeed.mat'])
TwoLayerSpeed = reshape(DriftSpeed,1,numel(DriftSpeed));
load([DataDir,'OneLayerDriftSpeed.mat'])
OneLayerSpeed = reshape(DriftSpeed,1,numel(DriftSpeed));

f1 = figure;
figure(f1)
plot(OneLayerSpeed,OneLayerFI,'.','MarkerSize',8,'MarkerEdgeColor',[0.3,.53,.37]);
box off
set(gca,'FontSize',10,'TickLength',[0.025,0.01]);
xlim([0 max(OneLayerSpeed)*5/4]);
xValid =  [ones(length(OneLayerSpeed),1),OneLayerSpeed'];
b = xValid\OneLayerFI';
hold on
yFit = b(2)*[0;OneLayerSpeed']+b(1);
plot([0;OneLayerSpeed'],yFit,'LineWidth',1.1,'Color',[0.2,0.353,0.67]);
ylim([0 40]);
hold off
[r,p] = corr(OneLayerSpeed',OneLayerFI');
title(['r = ',num2str(r,'%.2f'),', p = ',num2str(p,'%0.1e')],'FontSize',6);
set(gcf,'Units','Centimeters','Position',[2,2,6,6]);
set(gca,'FontSize',10,'LooseInset',[0 0 0 0],'TickDir','out','TickLength',[0.025,0.01],'LineWidth',0.8);
xlabel('Speed ($^\circ$/s)','Interpreter','latex');
ylabel('FI (deg^{-2})');
xlim([0 4])

f2 = figure;
figure(f2)
plot(TwoLayerSpeed,TwoLayerFI,'.','MarkerSize',8,'MarkerEdgeColor',[0.3,.53,.37]);
box off
set(gca,'FontSize',10,'TickLength',[0.025,0.01]);
xlim([0 max(TwoLayerSpeed)*5/4]);
xValid =  [ones(length(TwoLayerSpeed),1),TwoLayerSpeed'];
b = xValid\TwoLayerFI';
hold on
yData = TwoLayerFI;
yFit = b(2)*[0;TwoLayerSpeed']+b(1);
plot([0;TwoLayerSpeed'],yFit,'LineWidth',1.1,'Color',[0.2,0.353,0.67]);
ylim([0 50]);
hold off
[r,p] = corr(TwoLayerSpeed',TwoLayerFI');
title(['r = ',num2str(r,'%.2f'),', p = ',num2str(p,'%0.1e')],'FontSize',6);
set(gcf,'Units','Centimeters','Position',[2,2,6,6]);
set(gca,'FontSize',10,'LooseInset',[0 0 0 0],'TickDir','out','TickLength',[0.025,0.01],'LineWidth',0.8);
xlabel('Speed ($^\circ$/s)','Interpreter','latex');
ylabel('FI (deg^{-2})');
xlim([0 max(TwoLayerSpeed)+0.05])