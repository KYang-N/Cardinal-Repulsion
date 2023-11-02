%% Simulation of the one-layer network with Poisson noise
%% Fisher information

clear
close
LoadFlag = 1;
SaveFlag = 0;
%% Define parameters

if LoadFlag == 0
    seed = 123;
    rng(seed)

    MemoryNet.alpha = 5e-4;
    MemoryNet.beta = 2.4e-3;
    MemoryNet.Mode = 'Both';

    MemoryNet.N = 300;

    % Simulation parameters
    DynParams.AddNoise = 1;
    DynParams.dt = 1e-3;
    DynParams.StimTime = 0.5;
    DynParams.RepTime = 3e3;
    DynParams.PFTime = 5;
    DynParams.DecodeTime = [1,2.5,4] + DynParams.StimTime;
%% Dynamics of the full SDE

    p = parpool(DynParams.Parallel);
    [FI,MemoryNet,DynParams] = OneLayerNetworkFI(MemoryNet,DynParams);
    delete(p);
    if SaveFlag
        disp('Saving the data.')
        save(['OneLayerFIalpha',strrep(num2str(MemoryNet.alpha),'.','p'), ...
            'beta',strrep(num2str(MemoryNet.beta),'.','p'),'Mode',num2str(MemoryNet.Mode),'.mat']);
    end
end
%%
if LoadFlag
    clear %#ok<*UNRCH>
    DataDir = '';
    load([DataDir,'OneLayerFIalpha0p0005beta0p0024ModeBoth.mat']);
end
close
Color = [0.6667 0.6667 0.6667;0.3333 0.3333 0.3333;0 0 0];

SampleInput = 0:DynParams.dSample:2*pi;
FigOutDir = '';
f = figure;
figure(f)
for Time = 1:3
    %     plot(SampleInput,FI(Time,:)/max(FI(1,:)),'Color',Color(Time,:),'LineWidth',1.5);
    plot(SampleInput,FI(Time,:),'Color',Color(Time,:),'LineWidth',1.5);
    hold on
end
hold off
xlabel('$\theta$ ($^\circ$)','Interpreter','latex');
ylabel('FI (deg^{-2})');
set(gca,'FontSize',10,'LineWidth',0.8,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0, ...
    'LooseInset',[0 0 0 0]);
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
box off
% yticks([0,1]);
% ylim([0,1.05]);
set(gcf,'Unit','Centimeters','Position',[2,2,3.7,3.7]);
