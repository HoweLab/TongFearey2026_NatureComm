%% figure_closedloop_standard_task_velocity_structuring.m
% Publication-oriented figure generator for CLOSED-LOOP standard task
%
% Generates:
% 1. Representative-simulation colormap of position-binned velocity across trials
% 2. Representative-simulation colormap of position-binned dSPN-iSPN imbalance across trials
% 3. Early vs late position-binned velocity, with 95% CI across simulations
% 4. Early vs late position-binned acceleration, with 95% CI across simulations
%
% Core model:
% - Continuous locomotion until reward
% - Local dSPN/iSPN imbalance modulates velocity/vigor
% - Learning signal = smoothed, lagged signed acceleration
% - No reward-related dopamine signal

clear; clc; close all;

%% =========================
% USER SETTINGS
% =========================
nSims            = 50;     % number of simulations for CI plots
baseSeed         = 3;
repSimIndex      = 1;      % representative simulation for colormaps

nTrials          = 150;
nEarly           = 10;
nLate            = 10;

saveFigures      = false;
saveDir          = pwd;

%% =========================
% GLOBAL PARAMETERS
% =========================
% Time / track
dt               = 0.05;      % s
trackLength      = 1.5;       % m
dxBin            = 0.02;      % m
posEdges         = 0:dxBin:trackLength;
posCenters       = posEdges(1:end-1) + dxBin/2;
nBins            = numel(posCenters);

maxTrialTime     = 150;       % includes post-reward dwell
maxSteps         = ceil(maxTrialTime / dt);

% Network
nNeurons         = 1000;
nD               = 500;
nI               = nNeurons - nD;
dIdx             = 1:nD;
iIdx             = (nD+1):nNeurons;
cellType         = [ones(nD,1); -ones(nI,1)];   % +1 dSPN, -1 iSPN

% Position fields
sigmaPos         = 0.12;

% Fixed velocity and positive-acceleration input weights
vNormMax         = 0.35;
baselineMean     = 0.02;
baselineSD       = 0.01;
noiseSD          = 0.02;

% Learning
etaPos           = 0.004;
tauEligibility   = 1.0;
lambdaE          = exp(-dt / tauEligibility);

daLagSec         = 0.20;
daLagSteps       = round(daLagSec / dt);

tauDA            = 0.80;
lambdaDA         = exp(-dt / tauDA);

% Reward / dwell
rewardPos         = trackLength;
postRewardDwell   = 3.0;

% Continuous locomotion controller
vCruiseBase       = 0.16;
tauV              = 0.35;
sigmaRunNoise     = 0.015;
vMaxCap           = 0.35;
vMinPreReward     = 0.03;

% Reward-zone slowing before reward
rewardZoneStart   = 1.40;
rewardSlowGain    = 0.08;

% After reward, allow stopping
tauStopPostReward = 0.20;

% Feedback from local imbalance to velocity
kImbVel           = 0.12;

%% =========================
% STORAGE ACROSS SIMULATIONS
% =========================
earlyVelAll   = nan(nSims, nBins);
lateVelAll    = nan(nSims, nBins);
earlyAccAll   = nan(nSims, nBins);
lateAccAll    = nan(nSims, nBins);

rep = struct();

%% =========================
% HELPERS
% =========================
clip01  = @(x) min(max(x,0),1);
gauss1d = @(x,mu,sig) exp(-((x-mu).^2) ./ (2*sig.^2));

%% =========================
% MAIN SIMULATION LOOP
% =========================
for sim = 1:nSims
    rng(baseSeed + sim - 1);

    % -------------------------
    % INITIALIZE NETWORK
    % -------------------------
    muPos    = rand(nNeurons,1) * trackLength;

    wVel     = rand(nNeurons,1);

    wAcc     = zeros(nNeurons,1);
    % accStrongMask    = rand(nNeurons,1) < 0.30;
    % wAcc(accStrongMask)  = 0.5 + 0.5*rand(sum(accStrongMask),1);
    % wAcc(~accStrongMask) = 0.0 + 0.1*rand(sum(~accStrongMask),1);

    wPos     = rand(nNeurons,1);

    baseline = baselineMean + baselineSD * randn(nNeurons,1);

    % -------------------------
    % PER-TRIAL STORAGE
    % -------------------------
    taskVelByTrial = nan(nTrials, nBins);
    taskAccByTrial = nan(nTrials, nBins);
    taskImbByTrial = nan(nTrials, nBins);

    % -------------------------
    % TRIAL LOOP
    % -------------------------
    for tr = 1:nTrials

        x = nan(maxSteps,1);
        v = nan(maxSteps,1);
        a = nan(maxSteps,1);
        R = nan(nNeurons, maxSteps);
        imbTime = nan(1, maxSteps);
        accSmooth = nan(1, maxSteps);

        ePos = zeros(nNeurons,1);
        daTrace = 0;

        rewardDelivered = false;
        rewardTime = nan;

        x(1) = 0;
        v(1) = 0.05;
        a(1) = 0;

        for t = 1:maxSteps
            xCurr = x(t);
            tSec  = (t-1)*dt;

            P = gauss1d(xCurr, muPos, sigmaPos);

            if t == 1
                vCurr = v(1);
                aCurr = 0;
            else
                vCurr = v(t);
                aCurr = a(t);
            end

            vNorm = min(vCurr / vNormMax, 1);
            aNorm = aCurr / 1.5;
            aPlus = max(aNorm, 0);

            R(:,t) = max(0, ...
                wPos .* P + ...
                wVel .* vNorm + ...
                wAcc .* aPlus + ...
                baseline + noiseSD * randn(nNeurons,1));

            % Weighted local imbalance readout
            Pd = P(dIdx);
            Pi = P(iIdx);

            if sum(Pd) > eps
                dLocal = sum(R(dIdx,t) .* Pd) / sum(Pd);
            else
                dLocal = mean(R(dIdx,t));
            end

            if sum(Pi) > eps
                iLocal = sum(R(iIdx,t) .* Pi) / sum(Pi);
            else
                iLocal = mean(R(iIdx,t));
            end

            imbLocal = dLocal - iLocal;
            imbTime(t) = imbLocal;

            % Reward detection
            if ~rewardDelivered && xCurr >= rewardPos
                rewardDelivered = true;
                rewardTime = tSec;
            end

            % ----------------------------------------
            % CLOSED-LOOP VELOCITY CONTROLLER
            % ----------------------------------------
            if ~rewardDelivered
                vTarget = vCruiseBase;

                if xCurr >= rewardZoneStart
                    frac = (xCurr - rewardZoneStart) / max(trackLength - rewardZoneStart, eps);
                    vTarget = vTarget - rewardSlowGain * frac;
                end

                % d/i balance modulates velocity/vigor
                vTarget = vTarget + kImbVel * imbLocal;
                vTarget = min(max(vTarget, vMinPreReward), vMaxCap);

                if t < maxSteps
                    vNext = vCurr + (dt/tauV) * (vTarget - vCurr) + sigmaRunNoise * randn;
                    vNext = min(max(vNext, vMinPreReward), vMaxCap);
                end

            else
                if t < maxSteps
                    if tSec <= rewardTime + postRewardDwell
                        vNext = max(0, vCurr - (dt/tauStopPostReward) * vCurr + sigmaRunNoise * randn);
                        vNext = min(max(vNext, 0), vMaxCap);
                    else
                        vNext = 0;
                    end
                end
            end

            if t < maxSteps
                v(t+1) = vNext;
                x(t+1) = x(t) + v(t+1) * dt;
                a(t+1) = (v(t+1) - v(t)) / dt;
            end

            % ----------------------------------------
            % LEARNING SIGNAL
            % ----------------------------------------
            daTrace = lambdaDA * daTrace + (1 - lambdaDA) * aNorm;
            accSmooth(t) = daTrace;

            if t > daLagSteps
                daSignal = accSmooth(t - daLagSteps);
                daSignal = max(min(daSignal, 1.5), -1.5);
            else
                daSignal = 0;
            end

            % Plasticity
            ePos = lambdaE * ePos + (1 - lambdaE) * P;
            dWPos = etaPos * cellType .* daSignal .* ePos;
            wPos = clip01(wPos + dWPos);

            % End trial after post-reward dwell
            if rewardDelivered && tSec >= rewardTime + postRewardDwell
                x = x(1:t);
                v = v(1:t);
                a = a(1:t);
                R = R(:,1:t);
                imbTime = imbTime(1:t);
                break;
            end
        end

        % Robust truncation
        validT = find(~isnan(x), 1, 'last');
        x = x(1:validT);
        v = v(1:validT);
        a = a(1:validT);
        R = R(:,1:validT);
        imbTime = imbTime(1:validT);

        if max(x) < rewardPos
            warning('Simulation %d, trial %d did not reach reward. max(x)=%.3f', sim, tr, max(x));
        end

        % Position-binned summaries
        velBin = nan(1,nBins);
        accBin = nan(1,nBins);
        imbBin = nan(1,nBins);

        for b = 1:nBins
            if b < nBins
                inBin = x >= posEdges(b) & x < posEdges(b+1);
            else
                inBin = x >= posEdges(b) & x <= posEdges(b+1);
            end

            if any(inBin)
                velBin(b) = mean(v(inBin), 'omitnan');
                accBin(b) = mean(a(inBin), 'omitnan');

                dBin = mean(R(dIdx, inBin), 'all');
                iBin = mean(R(iIdx, inBin), 'all');
                imbBin(b) = dBin - iBin;
            end
        end

        taskVelByTrial(tr,:) = velBin;
        taskAccByTrial(tr,:) = accBin;
        taskImbByTrial(tr,:) = imbBin;
    end

    % -------------------------
    % EARLY/LATE SUMMARY FOR THIS SIM
    % -------------------------
    earlyVelAll(sim,:) = mean(taskVelByTrial(1:nEarly,:), 1, 'omitnan');
    lateVelAll(sim,:)  = mean(taskVelByTrial(end-nLate+1:end,:), 1, 'omitnan');

    earlyAccAll(sim,:) = mean(taskAccByTrial(1:nEarly,:), 1, 'omitnan');
    lateAccAll(sim,:)  = mean(taskAccByTrial(end-nLate+1:end,:), 1, 'omitnan');

    % -------------------------
    % REPRESENTATIVE SIMULATION
    % -------------------------
    if sim == repSimIndex
        rep.taskVelByTrial = taskVelByTrial;
        rep.taskImbByTrial = taskImbByTrial;
    end

    fprintf('Completed simulation %d / %d\n', sim, nSims);
end

%% =========================
% SUMMARY STATISTICS
% =========================
earlyVel_mean = mean(earlyVelAll, 1, 'omitnan');
lateVel_mean  = mean(lateVelAll, 1, 'omitnan');
earlyVel_ci   = prctile(earlyVelAll, [2.5 97.5], 1);
lateVel_ci    = prctile(lateVelAll,  [2.5 97.5], 1);

earlyAcc_mean = mean(earlyAccAll, 1, 'omitnan');
lateAcc_mean  = mean(lateAccAll, 1, 'omitnan');
earlyAcc_ci   = prctile(earlyAccAll, [2.5 97.5], 1);
lateAcc_ci    = prctile(lateAccAll,  [2.5 97.5], 1);

%% =========================
% FIGURE PANELS
% =========================

% Panel 1: representative velocity colormap
fig1 = figure('Color','w','Name','Panel1_VelocityColormap','Position',[100 100 540 520]);
imagesc(posCenters, 1:nTrials, rep.taskVelByTrial);
set(gca, 'YDir', 'normal');
xlabel('Position (m)');
ylabel('Trial');
title('Velocity across learning');
colorbar;

% Panel 2: representative imbalance colormap
fig2 = figure('Color','w','Name','Panel2_ImbalanceColormap','Position',[100 100 540 520]);
imagesc(posCenters, 1:nTrials, rep.taskImbByTrial);
set(gca, 'YDir', 'normal');
xlabel('Position (m)');
ylabel('Trial');
title('dSPN - iSPN across learning');
colorbar;

% Panel 3: early vs late velocity
fig3 = figure('Color','w','Name','Panel3_EarlyLateVelocity','Position',[100 100 540 420]);
hold on;
plot_ci_shaded(posCenters, earlyVel_mean, earlyVel_ci, [0.5 0.5 0.5], '-');
plot_ci_shaded(posCenters, lateVel_mean,  lateVel_ci,  [0 0 0], '-');
xlabel('Position (m)');
ylabel('Velocity (m/s)');
title(sprintf('Velocity: early (first %d) vs late (last %d)', nEarly, nLate));
legend({'Early','Late'}, 'Location', 'best');
xlim([0 trackLength]);

% Panel 4: early vs late acceleration
fig4 = figure('Color','w','Name','Panel4_EarlyLateAcceleration','Position',[100 100 540 420]);
hold on;
plot_ci_shaded(posCenters, earlyAcc_mean, earlyAcc_ci, [0.5 0.5 0.5], '-');
plot_ci_shaded(posCenters, lateAcc_mean,  lateAcc_ci,  [0 0 0], '-');
yline(0,'k--');
xlabel('Position (m)');
ylabel('Acceleration (m/s^2)');
title(sprintf('Acceleration: early (first %d) vs late (last %d)', nEarly, nLate));
legend({'Early','Late'}, 'Location', 'best');
xlim([0 trackLength]);

%% =========================
% OPTIONAL SAVE
% =========================
if saveFigures
    exportgraphics(fig1, fullfile(saveDir, 'panel1_velocity_colormap.pdf'), 'ContentType','vector');
    exportgraphics(fig2, fullfile(saveDir, 'panel2_imbalance_colormap.pdf'), 'ContentType','vector');
    exportgraphics(fig3, fullfile(saveDir, 'panel3_early_late_velocity.pdf'), 'ContentType','vector');
    exportgraphics(fig4, fullfile(saveDir, 'panel4_early_late_acceleration.pdf'), 'ContentType','vector');
end

%% =========================
% LOCAL FUNCTIONS
% =========================
function plot_ci_shaded(x, yMean, yCI, lineColor, lineStyle)
    fill([x fliplr(x)], [yCI(1,:) fliplr(yCI(2,:))], lineColor, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none');
    plot(x, yMean, 'Color', lineColor, 'LineWidth', 2.2, 'LineStyle', lineStyle);
end