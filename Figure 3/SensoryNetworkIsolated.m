%% Implement the sensory network in islolation

clear
close
%%
N = 300;
Tau_syn = .01;
dtheta = 2*pi/N; % We map 180 degs to 2*pi
theta = 0:dtheta:2*pi-dtheta;
%% Generate connectivity

JI = .35;
JE = .60;
alpha = -0.07;
beta = 0;

lambda = 0.36*pi;
lambdaI = 1.1*pi;
Conn = zeros(N,N);

for i = 1:N
    thetaCurrent = theta(i);
    JEModulated = JE*ConnModulation(thetaCurrent,alpha);
    JIModulated = JI*ConnModulation(thetaCurrent,beta);
    J = (-JI+JEModulated*exp(-((theta-pi)/lambda).^2))/(2*pi);
    Conn(i,:) = dtheta*circshift(J,[0 round(-pi/dtheta-1+i)]);
end
%% Dynamics

dt = 1e-3;
tmax = 0.5; % s
PFTime = tmax;
NInputSample = 50;
dSample = 2*pi/NInputSample;
SampleInput = 0:dSample:2*pi;
NInputSample = NInputSample + 1;

step = round(tmax/dt);
epsilon = 0.2;

C = 4;
mu = 0.3*pi;
PlotTime = 0.5;

% Transfer function
NE = 2; thE = 0.1; sigE = 6; maxfE = 100;
q = @(x) maxfE*(x-thE).^NE./(sigE^NE+(x-thE).^NE).*(x>thE);

TC = zeros(N,NInputSample);
mSensory_old = zeros(N,1);
sSensory_old = zeros(N,1);

tic
for jj = 1:NInputSample
    I0 = C*(1-2*epsilon+2*epsilon*exp(-((theta-pi)/mu).^2));
    I0 = circshift(I0,round(((-pi+SampleInput(jj))/dtheta)));
    m = zeros(N,tmax/dt);
    for ii = 1:step
        sSensory_new = sSensory_old + dt*1/Tau_syn*(-sSensory_old+mSensory_old);
        SensoryInputI = Conn*sSensory_old + I0';
        mSensory_new = q(SensoryInputI);
        mSensory_old = mSensory_new;
        sSensory_old = sSensory_new;
        if ii == round(PlotTime/dt)
            TC(:,jj) = mSensory_new;
        end
    end
end
toc
%% Plot result

f1 = figure;
figure(f1)
plot(theta,Conn(1:20:end,:)/dtheta,'LineWidth',1,'Color','#448983');
xlim([0,2*pi]);
yticks([-0.05 0 0.05])
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
set(gca,'FontSize',10,'TickLength',[0.025,0.01],'TickDir','out','XTickLabelRotation',0, ...
    'LooseInset',[0 0 0 0],'LineWidth',0.8);
box off
xlabel('$\psi$ ($^\circ$)','Interpreter','latex');
set(gcf,'Unit','Centimeters','Position',[2,2,5,3]);
ylabel('Conn.')

%%
f2 = figure;
figure(f2)
NeuronIdx = 1:20:N;
plot(SampleInput,TC(NeuronIdx,:),'LineWidth',1);
set(gca,'FontSize',10);
xlim([0 2*pi]);
box off
ylim([0,60]);
xlabel('$\theta$ ($^\circ$)','Interpreter','latex');
set(gca,'FontSize',10,'TickDir','out','TickLength',[0.025,.01],'LineWidth',.8, ...
    'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
ylabel('FR (Hz)')
set(gcf,'Unit','Centimeters','Position',[2,2,5,3]);

f3 = figure;
figure(f3)
PFSensory = ComputePF(TC,SampleInput);
histogram(PFSensory,30,'EdgeAlpha',0.5,'FaceColor',[0.5,0.5,0.5]);
box off
xlim([0 2*pi]);
ylabel('Counts','FontSize',10);
set(gca,'FontSize',10,'TickDir','out','TickLength',[0.025,.01],'LineWidth',.8, ...
        'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
xlabel('$\tilde{\theta}$','FontSize',10,'Interpreter','latex');
set(gcf,'Unit','Centimeters','Position',[2,2,5,2.5]);


WidthColor = '#448983';
f4 = figure;
figure(f4)
Widths = zeros(1,N);
for k = 1:N
    Widths(k) = ComputeTuningWidth(SampleInput,TC(k,:));
end
plot(theta,Widths,'LineWidth',2,'Color',WidthColor);
set(gca,'FontSize',10,'TickDir','out','TickLength',[0.025,.01],'LineWidth',.8, ...
        'XTickLabelRotation',0,'LooseInset',[0 0 0 0]);
ylim([0 55]);
xlim([0 2*pi]);
box off
xlabel('\psi')
xticks(0:pi/2:2*pi);
xticklabels({'0','','90','','180'});
ylabel('Width')
set(gcf,'Unit','Centimeters','Position',[2,2,5,2.5]);

%% Response manifold

RSteadyState = TC - mean(TC,2);
[U,S,V] = svd(RSteadyState');
fSVD2 = figure;

figure(fSVD2)
plot(V(:,1)'*RSteadyState,V(:,2)'*RSteadyState,'-','LineWidth',2,'Color', ...
    [0.7,0.7,0.7],'Marker','.','MarkerSize',14,'MarkerEdgeColor','#448983');
ylim([min(V(:,1)'*RSteadyState) max(V(:,1)'*RSteadyState)]);
xlim([min(V(:,1)'*RSteadyState) max(V(:,1)'*RSteadyState)]);
axis square;
axis off
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);

%%
fSVD3 = figure;
figure(fSVD3)

figure(fSVD3)
plot3(V(:,1)'*RSteadyState,V(:,2)'*RSteadyState,V(:,3)'*RSteadyState,'-', ...
    'LineWidth',3,'Color','#448983');

axis square;
grid on
xlim([-260 260])
ylim([-260 260])
zlim([-260 260])
xticklabels({'','',''})
yticklabels({'','',''})
zticklabels({'','',''})
set(gcf,'Unit','Centimeters','Position',[2,2,5,5]);
view([-19 43])
