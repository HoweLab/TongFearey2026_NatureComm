%% figure_openloop_standard_task_corrected_with_kinematics_ctx1only.m
% Publication-oriented figure generator for the OPEN-LOOP standard 2-context task
%
% Context-1-only revision:
% - Colormap uses context 1 only
% - Early/late imbalance summaries use context 1 only
% - Post-learning activity curves use late context 1 trials only
% - Velocity and acceleration panels use late context 1 trials only
% - Position-tuned vs non-position-tuned groups are defined from late context-1 activity

clear; clc; close all;

%% =========================
% USER SETTINGS
% =========================
nSims            = 50;
baseSeed         = 100;
repSimIndex      = 1;

nTrials          = 150;
nEarly           = 10;
nLate            = 10;

groupFrac        = 0.25;
saveFigures      = false;
saveDir          = pwd;

%% =========================
% GLOBAL MODEL PARAMETERS
% =========================
dt               = 0.05;
trackLength      = 1.5;
dxBin            = 0.02;
posEdges         = 0:dxBin:trackLength;
posCenters       = posEdges(1:end-1) + dxBin/2;
nBins            = numel(posCenters);

maxTrialTime     = 150;
maxSteps         = ceil(maxTrialTime / dt);

% Network
nNeurons         = 1000;
nD               = 500;
nI               = 500;
dIdx             = 1:nD;
iIdx             = (nD+1):nNeurons;
cellType         = [ones(nD,1); -ones(nI,1)];

% Inputs
sigmaPos         = 0.12;
sigmaDist        = 0.12;
ctxGroupP        = [0.35, 0.35, 0.30];

vNormMax         = 0.35;

baselineMean     = 0.02;
baselineSD       = 0.01;
noiseSD          = 0.02;

% Learning
%0.004 default
etaPos           = 0.008;
etaDist          = 0.004;
tauEligibility   = 1.0;
lambdaE          = exp(-dt / tauEligibility);

daLagSec         = 0.20;
daLagSteps       = round(daLagSec / dt);
tauDA            = 0.80;
lambdaDA         = exp(-dt / tauDA);

% Reward DA
includeRewardDA   = false;
rewardDAPeakDelay = 0.50;
rewardDAWidth     = 0.35;
rewardDAGain      = 0.60;
postRewardDwell   = 3.0;

% Standard task gains
posGain          = 1.0;
distGain         = 0.20;
posLearnGain     = 1.0;
distLearnGain    = 0.20;

% Kinematics
vMaxMean         = 0.20;
vMaxSD           = 0.02;
accelLenMean     = 0.20;
accelLenSD       = 0.03;
decelLenMean     = 0.20;
decelLenSD       = 0.03;

%% =========================
% STORAGE ACROSS SIMULATIONS
% =========================
earlyImbAll         = nan(nSims, nBins);
lateImbAll          = nan(nSims, nBins);

lateImbPosTunedAll  = nan(nSims, nBins);
lateImbNonPosAll    = nan(nSims, nBins);

preDAll             = nan(nSims, nBins);
preIAll             = nan(nSims, nBins);
postDAll            = nan(nSims, nBins);
postIAll            = nan(nSims, nBins);

lateVelAll          = nan(nSims, nBins);
lateAccAll          = nan(nSims, nBins);

rep = struct();

%% =========================
% MAIN SIMULATION LOOP
% =========================
for sim = 1:nSims
    rng(baseSeed + sim - 1);

    % -------------------------
    % INITIALIZE NETWORK
    % -------------------------
    muPosCtx1 = rand(nNeurons,1) * trackLength;
    muPosCtx2 = rand(nNeurons,1) * trackLength;
    muDist    = rand(nNeurons,1) * trackLength;

    ctxGroup  = randsample(1:3, nNeurons, true, ctxGroupP);

    gPosCtx1 = nan(nNeurons,1);
    gPosCtx2 = nan(nNeurons,1);
    for ii = 1:nNeurons
        switch ctxGroup(ii)
            case 1
                gPosCtx1(ii) = 0.8 + 0.2*rand;
                gPosCtx2(ii) = 0.02 + 0.08*rand;
            case 2
                gPosCtx1(ii) = 0.02 + 0.08*rand;
                gPosCtx2(ii) = 0.8 + 0.2*rand;
            case 3
                gPosCtx1(ii) = 0.45 + 0.20*rand;
                gPosCtx2(ii) = 0.45 + 0.20*rand;
        end
    end

    wVel = rand(nNeurons,1);

    wAcc = zeros(nNeurons,1);
    % accStrongMask = rand(nNeurons,1) < 0.30;
    % wAcc(accStrongMask)  = 0.5 + 0.5*rand(sum(accStrongMask),1);
    % wAcc(~accStrongMask) = 0.0 + 0.1*rand(sum(~accStrongMask),1);

    baseline = baselineMean + baselineSD * randn(nNeurons,1);

    wPosCtx1 = rand(nNeurons,1);
    wPosCtx2 = rand(nNeurons,1);
    wDist    = rand(nNeurons,1);

    wPosCtx1_0 = wPosCtx1;
    wDist_0    = wDist;

    % -------------------------
    % PER-TRIAL STORAGE
    % -------------------------
    dPosByTrial   = nan(nTrials, nBins);
    iPosByTrial   = nan(nTrials, nBins);
    imbPosByTrial = nan(nTrials, nBins);
    ctxByTrial    = nan(nTrials,1);

    velPosByTrial = nan(nTrials, nBins);
    accPosByTrial = nan(nTrials, nBins);

    lateNeuronPosSum_ctx1   = zeros(nNeurons, nBins);
    lateNeuronPosCount_ctx1 = zeros(1, nBins);

    % -------------------------
    % PRE-LEARNING ACTIVITY CURVES
    % -------------------------
    [preD, preI] = compute_population_activity_curve_smooth( ...
        wPosCtx1_0, wDist_0, wVel, wAcc, baseline, ...
        gPosCtx1, muPosCtx1, muDist, sigmaPos, sigmaDist, ...
        posGain, distGain, posCenters, ...
        dIdx, iIdx, vNormMax);

    preDAll(sim,:) = preD;
    preIAll(sim,:) = preI;

    % -------------------------
    % MAIN TRIAL LOOP
    % -------------------------
    for tr = 1:nTrials
        ctx = mod(tr-1,2) + 1;
        ctxByTrial(tr) = ctx;

        vMax   = max(0.10, vMaxMean + vMaxSD * randn);
        accLen = min(0.35, max(0.08, accelLenMean + accelLenSD * randn));
        decLen = min(0.35, max(0.08, decelLenMean + decelLenSD * randn));

        xAbs = nan(maxSteps,1);
        v    = nan(maxSteps,1);
        a    = nan(maxSteps,1);
        distFromStart = nan(maxSteps,1);
        R    = nan(nNeurons, maxSteps);
        daTime = nan(1, maxSteps);
        rewardDATime = nan(1, maxSteps);
        accSmooth = nan(1, maxSteps);

        ePosCtx1 = zeros(nNeurons,1);
        ePosCtx2 = zeros(nNeurons,1);
        eDist    = zeros(nNeurons,1);
        daTrace  = 0;

        rewardDelivered = false;
        rewardTime = nan;

        startPos = 0;
        rewardTriggerDist = trackLength;
        decelStartDist = max(0, rewardTriggerDist - decLen);

        xAbs(1) = startPos;
        v(1) = 0;
        a(1) = 0;
        distFromStart(1) = 0;

        for t = 1:maxSteps
            dCurr = distFromStart(t);
            xCurr = xAbs(t);
            tSec  = (t-1) * dt;

            if ~rewardDelivered
                if dCurr <= accLen
                    s = min(max(dCurr / accLen, 0), 1);
                    vTarget = 0.5 * (1 - cos(pi*s)) * vMax;
                elseif dCurr < decelStartDist
                    vTarget = vMax;
                elseif dCurr < rewardTriggerDist
                    s = min(max((dCurr - decelStartDist) / max(decLen, eps), 0), 1);
                    vTarget = 0.5 * (1 + cos(pi*s)) * vMax;
                else
                    vTarget = 0;
                end
            else
                vTarget = 0;
            end

            v(t) = max(0, vTarget + 0.005 * randn);

            if rewardDelivered
                v(t) = 0;
            end

            if t == 1
                a(t) = 0;
            else
                a(t) = (v(t) - v(t-1)) / dt;
            end

            if ~rewardDelivered && dCurr >= rewardTriggerDist
                rewardDelivered = true;
                rewardTime = tSec;
                v(t) = 0;
                if t > 1
                    a(t) = (v(t) - v(t-1)) / dt;
                else
                    a(t) = 0;
                end
            end

            vNorm = min(v(t) / vNormMax, 1);
            aNorm = a(t) / 1.5;
            aPlus = max(aNorm, 0);

            if ctx == 1
                P = gPosCtx1 .* exp(-((xCurr - muPosCtx1).^2) ./ (2*sigmaPos^2));
                wPosActive = wPosCtx1;
            else
                P = gPosCtx2 .* exp(-((xCurr - muPosCtx2).^2) ./ (2*sigmaPos^2));
                wPosActive = wPosCtx2;
            end

            D = exp(-((dCurr - muDist).^2) ./ (2*sigmaDist^2));

            R(:,t) = max(0, ...
                posGain  * (wPosActive .* P) + ...
                distGain * (wDist .* D) + ...
                wVel .* vNorm + ...
                wAcc .* aPlus + ...
                baseline + noiseSD * randn(nNeurons,1));

            daTrace = lambdaDA * daTrace + (1 - lambdaDA) * aNorm;
            accSmooth(t) = daTrace;

            if t > daLagSteps
                kinDASignal = accSmooth(t - daLagSteps);
                kinDASignal = max(min(kinDASignal, 1.5), -1.5);
            else
                kinDASignal = 0;
            end

            rewardDASignal = 0;
            if includeRewardDA && rewardDelivered
                dtReward = tSec - rewardTime;
                rewardDASignal = rewardDAGain * exp(-((dtReward - rewardDAPeakDelay)^2) / (2*rewardDAWidth^2));
            end
            rewardDATime(t) = rewardDASignal;

            daSignal = kinDASignal + rewardDASignal;
            daTime(t) = daSignal;

            if ctx == 1
                ePosCtx1 = lambdaE * ePosCtx1 + (1 - lambdaE) * P;
                ePosCtx2 = lambdaE * ePosCtx2;
            else
                ePosCtx2 = lambdaE * ePosCtx2 + (1 - lambdaE) * P;
                ePosCtx1 = lambdaE * ePosCtx1;
            end
            eDist = lambdaE * eDist + (1 - lambdaE) * D;

            if ctx == 1
                dWPosCtx1 = etaPos * posLearnGain * cellType .* daSignal .* ePosCtx1;
                wPosCtx1 = min(max(wPosCtx1 + dWPosCtx1, 0), 1);
            else
                dWPosCtx2 = etaPos * posLearnGain * cellType .* daSignal .* ePosCtx2;
                wPosCtx2 = min(max(wPosCtx2 + dWPosCtx2, 0), 1);
            end

            dWDist = etaDist * distLearnGain * cellType .* daSignal .* eDist;
            wDist = min(max(wDist + dWDist, 0), 1);

            if t < maxSteps
                if rewardDelivered
                    xAbs(t+1) = xAbs(t);
                    distFromStart(t+1) = distFromStart(t);
                else
                    xAbs(t+1) = mod(xAbs(t) + v(t)*dt, trackLength);
                    distFromStart(t+1) = distFromStart(t) + v(t)*dt;
                end
            end

            if rewardDelivered && tSec >= rewardTime + postRewardDwell
                xAbs = xAbs(1:t);
                v = v(1:t);
                a = a(1:t);
                distFromStart = distFromStart(1:t);
                R = R(:,1:t);
                daTime = daTime(1:t);
                break;
            end
        end

        validT = find(~isnan(xAbs), 1, 'last');
        xAbs = xAbs(1:validT);
        v    = v(1:validT);
        a    = a(1:validT);
        R    = R(:,1:validT);

        dPos   = nan(1,nBins);
        iPos   = nan(1,nBins);
        velPos = nan(1,nBins);
        accPos = nan(1,nBins);

        for b = 1:nBins
            if b < nBins
                inBin = xAbs >= posEdges(b) & xAbs < posEdges(b+1);
            else
                inBin = xAbs >= posEdges(b) & xAbs <= posEdges(b+1);
            end

            if any(inBin)
                dPos(b)   = mean(R(dIdx, inBin), 'all');
                iPos(b)   = mean(R(iIdx, inBin), 'all');
                velPos(b) = mean(v(inBin), 'omitnan');
                accPos(b) = mean(a(inBin), 'omitnan');
            end
        end

        dPosByTrial(tr,:)   = dPos;
        iPosByTrial(tr,:)   = iPos;
        imbPosByTrial(tr,:) = dPos - iPos;

        velPosByTrial(tr,:) = velPos;
        accPosByTrial(tr,:) = accPos;

        isLateTrial = tr > (nTrials - nLate);
        if isLateTrial && ctx == 1
            for b = 1:nBins
                if b < nBins
                    inBin = xAbs >= posEdges(b) & xAbs < posEdges(b+1);
                else
                    inBin = xAbs >= posEdges(b) & xAbs <= posEdges(b+1);
                end

                if any(inBin)
                    lateNeuronPosSum_ctx1(:,b) = lateNeuronPosSum_ctx1(:,b) + mean(R(:, inBin), 2);
                    lateNeuronPosCount_ctx1(b) = lateNeuronPosCount_ctx1(b) + 1;
                end
            end
        end
    end

    % Context-1-only trial indices
    ctx1TrialIdx = find(ctxByTrial == 1);
    earlyCtx1Idx = ctx1TrialIdx(1:min(nEarly, numel(ctx1TrialIdx)));
    lateCtx1Idx  = ctx1TrialIdx(max(1, numel(ctx1TrialIdx)-nLate+1):end);

    % Post-learning activity curves: context 1 only
    postDAll(sim,:) = mean(dPosByTrial(lateCtx1Idx,:), 1, 'omitnan');
    postIAll(sim,:) = mean(iPosByTrial(lateCtx1Idx,:), 1, 'omitnan');

    % Kinematic curves: context 1 only
    lateVelAll(sim,:) = mean(velPosByTrial(lateCtx1Idx,:), 1, 'omitnan');
    lateAccAll(sim,:) = mean(accPosByTrial(lateCtx1Idx,:), 1, 'omitnan');

    % Early/late imbalance: context 1 only
    earlyImbAll(sim,:) = mean(imbPosByTrial(earlyCtx1Idx,:), 1, 'omitnan');
    lateImbAll(sim,:)  = mean(imbPosByTrial(lateCtx1Idx,:), 1, 'omitnan');

    % Define subgroups from post-learning position tuning
    lateNeuronPosMean_ctx1 = lateNeuronPosSum_ctx1 ./ max(lateNeuronPosCount_ctx1, 1);

    centeredCurves = lateNeuronPosMean_ctx1 - mean(lateNeuronPosMean_ctx1, 2, 'omitnan');
    posTuningStrength = std(centeredCurves, 0, 2, 'omitnan');

    [~, sortPos] = sort(posTuningStrength, 'descend');
    nGroup = round(groupFrac * nNeurons);

    posTunedGroup = sortPos(1:nGroup);
    nonPosGroup   = sortPos(end-nGroup+1:end);

    dPosTuned = intersect(posTunedGroup, dIdx);
    iPosTuned = intersect(posTunedGroup, iIdx);
    dNonPos   = intersect(nonPosGroup, dIdx);
    iNonPos   = intersect(nonPosGroup, iIdx);

    latePosTunedCurve = nan(1,nBins);
    lateNonPosCurve   = nan(1,nBins);

    for b = 1:nBins
        if ~isempty(dPosTuned) && ~isempty(iPosTuned)
            latePosTunedCurve(b) = ...
                mean(lateNeuronPosMean_ctx1(dPosTuned, b), 'omitnan') - ...
                mean(lateNeuronPosMean_ctx1(iPosTuned, b), 'omitnan');
        end

        if ~isempty(dNonPos) && ~isempty(iNonPos)
            lateNonPosCurve(b) = ...
                mean(lateNeuronPosMean_ctx1(dNonPos, b), 'omitnan') - ...
                mean(lateNeuronPosMean_ctx1(iNonPos, b), 'omitnan');
        end
    end

    lateImbPosTunedAll(sim,:) = latePosTunedCurve;
    lateImbNonPosAll(sim,:)   = lateNonPosCurve;

    if sim == repSimIndex
        rep.imbPosByTrial = imbPosByTrial;
        rep.ctxByTrial = ctxByTrial;
        rep.ctx1TrialIdx = find(ctxByTrial == 1);
        rep.ctx1ImbPosByTrial = imbPosByTrial(rep.ctx1TrialIdx, :);
    end

    fprintf('Completed simulation %d / %d\n', sim, nSims);
end

%% =========================
% SUMMARY STATISTICS
% =========================
preD_mean  = mean(preDAll, 1, 'omitnan');
preI_mean  = mean(preIAll, 1, 'omitnan');
postD_mean = mean(postDAll, 1, 'omitnan');
postI_mean = mean(postIAll, 1, 'omitnan');

preD_ci    = prctile(preDAll,  [2.5 97.5], 1);
preI_ci    = prctile(preIAll,  [2.5 97.5], 1);
postD_ci   = prctile(postDAll, [2.5 97.5], 1);
postI_ci   = prctile(postIAll, [2.5 97.5], 1);

earlyImb_mean = mean(earlyImbAll, 1, 'omitnan');
lateImb_mean  = mean(lateImbAll, 1, 'omitnan');
earlyImb_ci   = prctile(earlyImbAll, [2.5 97.5], 1);
lateImb_ci    = prctile(lateImbAll,  [2.5 97.5], 1);

latePosTuned_mean = mean(lateImbPosTunedAll, 1, 'omitnan');
lateNonPos_mean   = mean(lateImbNonPosAll, 1, 'omitnan');
latePosTuned_ci   = prctile(lateImbPosTunedAll, [2.5 97.5], 1);
lateNonPos_ci     = prctile(lateImbNonPosAll,  [2.5 97.5], 1);

lateVel_mean = mean(lateVelAll, 1, 'omitnan');
lateVel_ci   = prctile(lateVelAll, [2.5 97.5], 1);

lateAcc_mean = mean(lateAccAll, 1, 'omitnan');
lateAcc_ci   = prctile(lateAccAll, [2.5 97.5], 1);

allSubgroupVals = [latePosTuned_ci(:); lateNonPos_ci(:)];
ylSub = [min(allSubgroupVals), max(allSubgroupVals)];

%% =========================
% FIGURE PANELS
% =========================

fig1a = figure('Color','w','Name','Panel1a_PreLearningActivity','Position',[100 100 520 420]);
hold on;
plot_ci_shaded(posCenters, preD_mean, preD_ci, [0.1 0.55 0.1], '-');
plot_ci_shaded(posCenters, preI_mean, preI_ci, [0.75 0.1 0.75], '-');
xlabel('Position (m)');
ylabel('Population activity (a.u.)');
title('Pre-learning');
legend({'dSPN','iSPN'}, 'Location', 'best');
xlim([0 trackLength]);

fig1b = figure('Color','w','Name','Panel1b_PostLearningActivity','Position',[100 100 520 420]);
hold on;
plot_ci_shaded(posCenters, postD_mean, postD_ci, [0.1 0.55 0.1], '-');
plot_ci_shaded(posCenters, postI_mean, postI_ci, [0.75 0.1 0.75], '-');
xlabel('Position (m)');
ylabel('Population activity (a.u.)');
title(sprintf('Post-learning context 1 (last %d context-1 trials)', nLate));
legend({'dSPN','iSPN'}, 'Location', 'best');
xlim([0 trackLength]);

fig2 = figure('Color','w','Name','Panel2_ImbalanceColormap','Position',[100 100 520 520]);
imagesc(posCenters, 1:size(rep.ctx1ImbPosByTrial,1), rep.ctx1ImbPosByTrial);
set(gca, 'YDir', 'normal');
xlabel('Position (m)');
ylabel('Context 1 trial #');
title('dSPN - iSPN across learning');
colorbar;

fig3 = figure('Color','w','Name','Panel3_EarlyLateImbalance','Position',[100 100 520 420]);
hold on;
plot_ci_shaded(posCenters, earlyImb_mean, earlyImb_ci, [0.5 0.5 0.5], '-');
plot_ci_shaded(posCenters, lateImb_mean,  lateImb_ci,  [0 0 0], '-');
yline(0,'k--');
xlabel('Position (m)');
ylabel('dSPN - iSPN');
title('Early vs late imbalance (context 1)');
legend({'Early ctx1','Late ctx1'}, 'Location', 'best');
xlim([0 trackLength]);

fig4 = figure('Color','w','Name','Panel4_TunedVsNonTuned_PostOnly','Position',[100 100 620 420]);
hold on;
plot_ci_shaded(posCenters, latePosTuned_mean, latePosTuned_ci, [0 0 0], '-');
plot_ci_shaded(posCenters, lateNonPos_mean,   lateNonPos_ci,   [0.2 0.45 0.85], '-');
yline(0,'k--');
xlabel('Position (m)');
ylabel('dSPN - iSPN');
title('Post-learning context 1: position-tuned vs non-position-tuned');
legend({'Position-tuned','Non-position-tuned'}, 'Location', 'best');
xlim([0 trackLength]);
ylim(ylSub);

fig5 = figure('Color','w','Name','Panel5_PostLearningVelocity','Position',[100 100 520 420]);
hold on;
plot_ci_shaded(posCenters, lateVel_mean, lateVel_ci, [0.1 0.1 0.1], '-');
xlabel('Position (m)');
ylabel('Velocity (m/s)');
title(sprintf('Late context-1 velocity profile (last %d context-1 trials)', nLate));
xlim([0 trackLength]);

fig6 = figure('Color','w','Name','Panel6_PostLearningAcceleration','Position',[100 100 520 420]);
hold on;
plot_ci_shaded(posCenters, lateAcc_mean, lateAcc_ci, [0.1 0.1 0.1], '-');
yline(0,'k--');
xlabel('Position (m)');
ylabel('Acceleration (m/s^2)');
title(sprintf('Late context-1 acceleration profile (last %d context-1 trials)', nLate));
xlim([0 trackLength]);

%% =========================
% OPTIONAL SAVE
% =========================
if saveFigures
    exportgraphics(fig1a, fullfile(saveDir, 'panel1a_pre_learning_activity.pdf'), 'ContentType','vector');
    exportgraphics(fig1b, fullfile(saveDir, 'panel1b_post_learning_activity_ctx1.pdf'), 'ContentType','vector');
    exportgraphics(fig2,  fullfile(saveDir, 'panel2_imbalance_colormap.pdf'), 'ContentType','vector');
    exportgraphics(fig3,  fullfile(saveDir, 'panel3_early_late_imbalance_ctx1.pdf'), 'ContentType','vector');
    exportgraphics(fig4,  fullfile(saveDir, 'panel4_tuned_vs_nontuned_postonly_ctx1.pdf'), 'ContentType','vector');
    exportgraphics(fig5,  fullfile(saveDir, 'panel5_post_learning_velocity_ctx1.pdf'), 'ContentType','vector');
    exportgraphics(fig6,  fullfile(saveDir, 'panel6_post_learning_acceleration_ctx1.pdf'), 'ContentType','vector');
end

%% =========================
% LOCAL FUNCTIONS
% =========================
function [dCurve, iCurve] = compute_population_activity_curve_smooth( ...
    wPosActive, wDist, wVel, wAcc, baseline, ...
    gPos, muPos, muDist, sigmaPos, sigmaDist, ...
    posGain, distGain, posCenters, ...
    dIdx, iIdx, vNormMax)

    nBins = numel(posCenters);
    dCurve = nan(1,nBins);
    iCurve = nan(1,nBins);

    n1 = round(nBins * 0.20);
    n2 = round(nBins * 0.60);
    n3 = nBins - n1 - n2;

    rise = linspace(0,1,max(n1,2));
    plateau = ones(1,max(n2,1));
    fall = linspace(1,0.05,max(n3,2));

    vProfile = [0.20*0.5*(1-cos(pi*rise)), 0.20*plateau, 0.20*fall];
    vProfile = vProfile(1:nBins);

    aProfile = [0 diff(vProfile)] / 0.05;

    for b = 1:nBins
        x = posCenters(b);
        d = x;

        P = gPos .* exp(-((x - muPos).^2) ./ (2*sigmaPos^2));
        D = exp(-((d - muDist).^2) ./ (2*sigmaDist^2));

        vNorm = min(vProfile(b) / vNormMax, 1);
        aPlus = max(aProfile(b) / 1.5, 0);

        R = max(0, ...
            posGain  * (wPosActive .* P) + ...
            distGain * (wDist .* D) + ...
            wVel .* vNorm + ...
            wAcc .* aPlus + ...
            baseline);

        dCurve(b) = mean(R(dIdx));
        iCurve(b) = mean(R(iIdx));
    end
end

function plot_ci_shaded(x, yMean, yCI, lineColor, lineStyle)
    fill([x fliplr(x)], [yCI(1,:) fliplr(yCI(2,:))], lineColor, ...
        'FaceAlpha', 0.18, 'EdgeColor', 'none');
    plot(x, yMean, 'Color', lineColor, 'LineWidth', 2.2, 'LineStyle', lineStyle);
end