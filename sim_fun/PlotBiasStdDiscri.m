function f1 = PlotBiasStdDiscri(DecodedOrientation,DecodeTime,StimTime,SampleInput,LegendFlag,Color)
% DecodedOrientation is a N*M*T matrix, where N is the number of
% realizations, M is the number of input samples, T is the number of sample
% time points.
%% Bias
if nargin == 5
Color = [0.4706    0.7765    0.4745;
    0.2549    0.6706    0.3647;
    0.1373    0.5176    0.2627];
end
NInputSample = size(DecodedOrientation,2);
RepTime = size(DecodedOrientation,1);

Bias = DecodedOrientation-SampleInput;

Bias(Bias>pi) = Bias(Bias>pi) - 2*pi;
Bias(Bias<-pi) = Bias(Bias<-pi) + 2*pi;
Bias = Bias/pi*180/2;

BiasMean = zeros(length(DecodeTime),NInputSample);
StdDv = BiasMean;
SEM = StdDv;

for Time = 1:length(DecodeTime)
    BiasMean(Time,:) = mean(Bias(:,:,Time),1);
    StdDv(Time,:) = std(Bias(:,:,Time),0,1);
    SEM(Time,:) = StdDv(Time,:)/sqrt(RepTime);
end

f1 = figure;
figure(f1)
subplot(1,2,1)
for Time = 1:length(DecodeTime)
    Upper = BiasMean(Time,:) + SEM(Time,:);
    Lower = BiasMean(Time,:) - SEM(Time,:);
    fill([SampleInput, fliplr(SampleInput)],[Upper, fliplr(Lower)],Color(Time,:),...
        'LineStyle','none','FaceAlpha',0.3);
    hold on
    plot(SampleInput,BiasMean(Time,:),'LineWidth',1.5,'Color',Color(Time,:));
end
hold off
yline(0,'LineStyle','--','Color','k','LineWidth',1.0);
xlabel('$\theta\; (^\circ)$','Interpreter','latex');
ylabel('Bias $(^\circ)$','Interpreter','latex');
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
xlim([0 2*pi]);
ylim([min(BiasMean(end,:)-1),max(BiasMean(end,:)+1)])
set(gca,'FontSize',10,'LineWidth',0.8,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0, ...
    'LooseInset',[0 0 0 0]);
box off
if LegendFlag
legend(gca,'',[num2str((DecodeTime(1)-StimTime)),'s'],'',[num2str((DecodeTime(2)-StimTime)),'s'],'',...
    [num2str((DecodeTime(3)-StimTime)),'s']);
legend('boxoff');
end

%% Std
Std = zeros(length(DecodeTime),NInputSample);
for Time = 1:length(DecodeTime)
    Std(Time,:) = std(Bias(:,:,Time),0,1);
    MeanSquare = mean(Bias(:,:,Time),1).^2;
    StdSample = sqrt(Bias(:,:,Time).^2-MeanSquare);
    StdDv(Time,:) = std(StdSample,0,1);
    SEM(Time,:) = StdDv(Time,:)/sqrt(RepTime);
end

figure(f1)
subplot(1,2,2)
for Time = 1:length(DecodeTime)
    Upper = Std(Time,:) + SEM(Time,:);
    Lower = Std(Time,:) - SEM(Time,:);
    fill([SampleInput, fliplr(SampleInput)],[Upper, fliplr(Lower)],Color(Time,:),...
        'LineStyle','none','FaceAlpha',0.3);
    hold on
    plot(SampleInput,Std(Time,:),'LineWidth',1.5,'Color',Color(Time,:));
end
hold off
set(gca,'FontSize',10,'LineWidth',0.8,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0, ...
    'LooseInset',[0 0 0 0]);
% legend('',[num2str((DecodeTime(1)-StimTime)/1000),'s'],'',[num2str((DecodeTime(2)-StimTime)/1000),'s'],'',...
%     [num2str((DecodeTime(3)-StimTime)/1000),'s']);
% legend('boxoff')
ylim([0 max(Std(end,:))+0.2]);
xlim([0 2*pi]);
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
box off
xlabel('$\theta\;(^\circ)$','Interpreter','latex');
ylabel('SD $(^\circ)$','Interpreter','latex');


% %% Discrimination Threshold
% subplot(1,3,3)
% for Time = 1:length(DecodeTime)
%     Discrimination = ComputeDT(BiasMean(Time,:),sqrt(Std(Time,:)),SampleInput);
%     plot(SampleInput,Discrimination,'Color',Color(Time,:),'LineWidth',1.5);
%     hold on
% end
% hold off
% set(gca,'FontSize',10,'LineWidth',.8,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0, ...
%     'LooseInset',[0 0 0 0]);
% xlim([0 2*pi]);
% xticks(0:pi/2:2*pi);
% xticklabels({'0','','90','','180'});
% yticks([]);
% box off
% xlabel('$\theta\; (^\circ)$','Interpreter','latex');
% ylabel('Threshold (a.u.)','Interpreter','latex');
set(f1,'Units','centimeters','Position',[4,2,9,4]);
