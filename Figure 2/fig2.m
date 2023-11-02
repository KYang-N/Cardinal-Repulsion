%% Initialize
clear;
clc;

% Save/load conditional distribution table, not required
save = false;
load = false;
outpath = '/';
fname = 'fig2_con';

%% Constants
len = 10000;
x = linspace(0, 2*pi, len+1);
prior = my_pdf(x);

%% Calculate mapping F
m = linspace(0, 2*pi, len+1);
fmap = zeros(size(x));

for i = 1:len+1
    fmap(i) = cumulative(i, len, prior) * (2*pi);
end

%% Create the conditional probability matrix
if ~load
    kappa = 1/250; % Set sensory noise (adjust as needed)
    sens_noise = 1/kappa;
    con = get_con_with_given_noise_level(sens_noise, fmap, fmap, x, m);
    
    if save
        data = struct('x', x, 'm', m, 'kappa', kappa, 'con', con); %#ok<*UNRCH> 
        savejson('', data, fullfile(outpath, [fname '.json']));
    end
else
    % Load con and overwrite x, m, length
    data = loadjson(fullfile(outpath, [fname '.json']));
    con = data.con;
    x = data.x;
    m = data.m;
    kappa = data.kappa;
    n_samples = length(x);
end

%% Visualize likelihood functions (cols of con matrix) in stimulus space
figure;
subplot(2, 1, 1);
plot(x/(2*pi)*360, prior);
title('prior');
xticks([]);
subplot(2, 1, 2);
likidx = round(linspace(0, 1, 13) * len + 1);
plot(x/(2*pi)*360, con(:, likidx));
xlabel('x = 2\theta');
xticks([0, 90, 180, 270, 360]);
title('likelihood p(m|x) in the stimulus space');

%% Visualize the con matrix
figure;
imagesc(con);
colorbar('Location', 'eastoutside');
title('\boldmath$con$(i,j) = $p(\mathbf{m}_j|\mathbf{x}_i)$', 'Interpreter', 'latex');

%% Visualize conditional probabilities (rows of con matrix)
figure;
subplot(2, 1, 1);
plot(x/(2*pi)*360, con(likidx, :));
xticks([0, 90, 180, 270, 360]);
xlabel('m');
title('distributions p(m|x) in the stimulus space');

%% Compute distributions of inferred x_h given x
n_post = 13;
idx_post = round(linspace(0, 1, n_post) * len + 1);
post_all = zeros(n_post, len+1);

%% Compute b_all and p_b_all
b_all = zeros(size(x));
for j = 1:(len+1)
    posterior = con(:, j)' .* prior;
    b_all(j) = bayes_mean(x, posterior);
    b_all(j) = b_all(j) - floor(b_all(j) / (2*pi)) * 2*pi;
end

b_all = round(((b_all - x) / (2*pi)) * (len+1));

p_b_all = zeros(len+1);
for j = 1:(len+1)
    p_b_all(j, b_all(j) + round((len+1)/2)) = 1;
end

for i = 1:n_post
    offset = mod(idx_post(i) - round((len+1)/2), len+1);
    con_centered = [con(idx_post(i), offset+1:end), con(idx_post(i), 1:offset)];
    convo = conv(con_centered, p_b_all(idx_post(i), :), 'same');
    post_all(i, :) = [convo(end-offset+1:end), convo(1:end-offset)];
    post_all(i, :) = post_all(i, :) / sum(post_all(i, :));
end

subplot(2, 1, 2);
plot(x/(2*pi)*360, post_all');
xticks([0, 90, 180, 270, 360]);
xlabel('$$\hat{x}$$','Interpreter','Latex');
title('distributions $$p(\hat{x}|x)$$ in the stimulus space','Interpreter','Latex');
%% Run simulation

% Run simulation without memory noise first
n_samples = 48;
points = linspace(0, 1, n_samples+1) * 2*pi;
n_trials = 10000;
n_epochs = 1;
out_epochs = [1]; % Indices of epochs to store results for
x_out_all_evo_nomemnoise = zeros(length(points), length(out_epochs), n_trials);

for i = 1:length(points)
    tic
    x_out_nomemnoise = run_inference_trials_evo(points(i), ...
                                                con, prior, ...
                                                out_epochs, n_trials, ...
                                                n_epochs, x, m, 0);
    x_out_all_evo_nomemnoise(i, :, :) = x_out_nomemnoise;
    elapsed_time = toc;
    fprintf('%d out of %d: %.2f sec\n', i, n_samples, elapsed_time);
end

[x_est_all_evo_nomemnoise, ...
 bias_all_evo_nomemnoise, ...
 var_all_evo_nomemnoise, ...
 std_all_evo_nomemnoise] = compute_bias_var_evo(x_out_all_evo_nomemnoise, points);

%% Add memory noise and simulate again
memsd = min(std_all_evo_nomemnoise(1, :));
fprintf('Memory standard deviation: %.4f radians or %.2f degrees\n', ...
         memsd, memsd * (180 / (2*pi)));

n_epochs = 5;
out_epochs = [1, 3, 5];
x_out_all_evo = zeros(length(points), length(out_epochs), n_trials);
x_sens_all_evo = zeros(length(points), length(out_epochs), n_trials);

for i = 1:length(points)
    tic
    [x_out, x_sens] = run_inference_trials_evo(points(i), ...
                                               con, prior, ...
                                               out_epochs, n_trials, ...
                                               n_epochs, x, m, memsd);
    x_out_all_evo(i, :, :) = x_out;
    x_sens_all_evo(i, :, :) = x_sens;
    elapsed_time = toc;
    fprintf('%d out of %d: %.2f sec\n', i, n_samples, elapsed_time);
end

[x_est_all_evo, ...
 bias_all_evo, ...
 var_all_evo, ...
 std_all_evo] = compute_bias_var_evo(x_out_all_evo, points);

%% Compute memory output distributions analytically
mem_noise = normpdf(x - max(x) / 2, 0, memsd);
conv_all = zeros(size(post_all));

for i = 1:n_post
    conv_all(i, :) = conv(mem_noise, post_all(i, :), 'same');
    conv_all(i, :) = conv_all(i, :) / sum(conv_all(i, :));
end
%% Fig 2 plots

% Load required tick label formatter
x_formatter = @(x) sprintf('%d', x);

% Define color formatting
COLOR = [0, 0, 0];
color_conv = @(color_range) COLOR + (color_range + 0.1)^0.5 - 0.2;

% Create figure and axes for the first set of subplots
figure;
subplot(2, 1, 1);
for c = 1:3
    plot(points, bias_all_evo(c, :) / (2*pi) * 180, ...
        'Color', color_conv((c-1)/2.2));
    hold on;
end
xticks([]);
yticks([-1, 0, 1]);
ylim([-1, 1]);
title('Bias');
ylabel('Bias ($^\circ$)', 'Interpreter',"latex");

subplot(2, 1, 2);
for c = 1:3
    circ_sd_c = std_all_evo(c, :) * (1/(2*pi) * 180);
    plot(points / (2*pi) * 180, circ_sd_c, ['' ...
        'Color'], color_conv((c-1)/2.2));
    hold on;
end
xticks([0, 45, 90, 135, 180]);
title('Standard deviation');
xlabel('Orientation ($$^\circ$$)','Interpreter',"latex");
ylabel('S.d. ($$^\circ$$)', 'Interpreter',"latex");
yticks([0, 3, 6]);
ylim([0, 8]);
set(gcf, 'Position', [100, 100, 400, 800]);


% Convert x to theta (x = 2*theta) and into degrees
x_to_theta = (x / (2*pi)) * 180;

% Calculate x axis limits xlow and xhigh
xl = floor(3 * (length(x) - 1) / 8) + 1;
xh = floor(7 * (length(x) - 1) / 8) + 1;

% Normalize the prior
prior_norm = prior / (3 * pi); % adjust as needed if you change the prior

% Create figure and axes for the second set of subplots
figure;
subplot(4, 1, 1);
xl = 3 * (length(x) - 1) / 8;
xh = 7 * (length(x) - 1) / 8 + 1;
prior_norm = prior / (3 * pi);
plot(x_to_theta(xl:xh), prior_norm(xl:xh), 'Color', 'black');
xticks([90, 135]);
yticks([0, 0.5]);
xlim([x_to_theta(xl), x_to_theta(xh)]);
ylim([0, 0.5]);
ylabel('Density');
title('Prior');

subplot(4, 1, 2);
likidx = round(linspace(0, 1, 13) * length(x));
likidx = likidx(6:11);
con_test_plot = con(likidx, :).';
c = [1    .55   .00
      .58 .39   .83
     1.00 .55   .00
      .39 .75  1.00
      .18 .55   .34
      .39 .75  1.00];
% set(gca, 'ColorOrder', c);
for i = 1:6
    plot(x_to_theta, con_test_plot(:,i) ./ max(con_test_plot(:,i)), 'Color', c(i,:));
    hold on;
end
xticks([90, 135]);
xlim([x_to_theta(xl), x_to_theta(xh)]);
yticks([0, 1]);
ylim([0, 1.01]);
ylabel('Likelihood');
title('Likelihood functions');

subplot(4, 1, 3);
yyaxis right
counts = zeros(1, 6);
for i = 6:11
    edges = [75 + (i - 6) * 15 - 15, 75 + (i - 6) * 15 + 15];
    h = histogram(x_sens_all_evo((i-1) * 4 + 1, 1, :) / (2*pi) * 180, ...
        'BinLimits', edges, 'NumBins', 60, 'DisplayStyle','stairs', 'EdgeColor',[.5,.5,.5]);
    counts(i - 5) = max(h.Values);
    hold on;
end
histmax = ceil(max(counts)/100)*100;
ylim([0, histmax]);
yticks([]);
xlim([x_to_theta(xl), x_to_theta(xh)]);
yyaxis left
set(gca, 'ColorOrder', c);
post_test_plot = post_all(6:11, :).';
for i = 1:6
    plot(x_to_theta, post_test_plot(:,i) * max(counts) / max(post_test_plot,[],'all'), ...
        'Color', c(i,:));
    hold on;
end
ylim([0, histmax]);
yticks([0, histmax]);
xticks([90, 135]);
xlim([x_to_theta(xl), x_to_theta(xh)]);
ylabel('Frequency / Scaled Density');
title('Distributions of $\hat{\theta}$ given $\theta$', 'Interpreter',"latex");

subplot(4, 1, 4);
yyaxis right
for i = 6:11
    edges = [75 + (i - 6) * 15 - 15, 75 + (i - 6) * 15 + 15];
    h = histogram(x_out_all_evo((i-1) * 4 + 1, 1, :) / (2*pi) * 180, ...
        'BinLimits', edges, 'NumBins', 60, 'DisplayStyle','stairs', 'EdgeColor',[.5,.5,.5]);
    hold on;
end
yticks([]);
ylim([0, histmax]);
xlim([x_to_theta(xl), x_to_theta(xh)]);
yyaxis left
set(gca, 'ColorOrder', c);
conv_test_plot = conv_all(6:11, :).';
for i = 1:6
    plot(x_to_theta, conv_test_plot(:,i) * max(counts) / max(post_test_plot,[],'all'), ...
        'Color', c(i,:));
    hold on;
end
yticks([0, histmax]);
ylim([0, histmax]);
xticks([90, 135]);
xlim([x_to_theta(xl), x_to_theta(xh)]);
xlabel('Orientation ($^\circ$)', 'Interpreter',"latex");
ylabel('Frequency / Scaled Density');
title('Distributions of $\theta_{i+1}$ given $\theta_i$', 'Interpreter',"latex");

set(gcf, 'Position', [100, 100, 400, 1000]);