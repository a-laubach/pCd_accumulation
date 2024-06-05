%%---------------Particulate Cd Accumulation Calculator--------------------
% Searches for points of cadmium accumulation by assuming subsurface
% particulate cadmium profiles are controlled by a combination of
% remineralization and an accumulation process.
% 
% First fits profiles with a 1-D remineralization equation that allows two
% pools of remineralizing cadmium plus a background inert pool:
%
%  Cd = Cd_inert + (Cd_lab0 * e^[k_lab * z]) +  (Cd_ref0 * e^[k_ref * z])
% 
% where z is depth below the euphotic zone and k values are decay constants.
%
% Accumulation is determined through rolling point-by-point removal of data
% and new remineralization fits to remaining data. Points for which removal
% creates better remineralization fits to remaining data are considered
% Cd-specific accumulation if in addition 1) original remineralization
% fits do not have R2 values above a threshold value 2) removed values are
% higher than remineralization fit values 3) the fraction of excess P is 
% below a threshold value.
%
%--------------------------------------------------------------------------
% This implementation requires the MATLAB Curve Fitting Toolbox for full 
% functionality. Developed in version 2023a. 
%
% A Laubach
% April 2024
% alaubach@ucsc.edu
%--------------------------------------------------------------------------


%% load data - adjust datapath as necessary
% pump small size fraction (SSF) pCd and pP from GP15 Station 33
load('GP15_Stn33_example.mat')


%% calculation options
opts.max_fit_depth = 2000; % [m]
opts.max_point_removal = 5; % max # of points to remove
opts.remin_R2_threshold = 0.9; % above this we do not consider accumulation
opts.max_qual_flag = 2; % follows SeaDataNet scheme
opts.do_P_removal = 1; % [0 or 1] option to do pP point removal/accumulation
opts.P_acc_threshold = 0.2; % frac P that disqualifies Cd-specific accumulation


%% plot output options
% [0 or 1]
plt.remin_only_prof = 1; % original remineralization profiles
plt.removal_prof = 1; % final removal profiles


%% Cd remineralization fit

% locate data - apply quality flag and depth constraints
remi = find(data.Cd_QF<=opts.max_qual_flag & data.Depth>=data.z_euph ...
           & data.Depth<=opts.max_fit_depth & data.P_QF<=opts.max_qual_flag);

Cd_meas = data.Cd_pM(remi);
Cd_max = max(data.Cd_pM(data.Cd_QF<=opts.max_qual_flag));

% define inert pool as minimum measured concentration
Cd_inert = min(data.Cd_pM(data.Depth>=data.z_euph & data.Cd_QF<=opts.max_qual_flag));

% calculate depth below euphotic zone for each data point
z_fit = data.Depth(remi) - data.z_euph(remi);

% remineralization fit to data
[Cd_remin_fit, Cd_remin_gof, Cd_remin_out] = fit(z_fit,Cd_meas-Cd_inert,...
    'exp2','Lower',[0,-0.5,0,-0.5],...
    'Upper',[Cd_max,0,Cd_max,0]);



%% P remineralization fit

% locate data - apply quality flag and depth constraints
P_meas = data.P_nM(remi);
P_max = max(data.P_nM(data.P_QF<=opts.max_qual_flag));

% define inert pool as minimum measured concentration
P_inert = min(data.P_nM(data.Depth>=data.z_euph & data.P_QF<=opts.max_qual_flag));

% calculate depth below euphotic zone for each data point
z_fit = data.Depth(remi) - data.z_euph(remi);

% remineralization fit to data
[P_remin_fit, P_remin_gof, P_remin_out] = fit(z_fit,P_meas-P_inert,...
    'exp2','Lower',[0,-0.5,0,-0.5],...
    'Upper',[P_max,0,P_max,0]);



%% calculate remineralization Cd:P R2

Cd_rem_vals = expVals(Cd_remin_fit,z_fit,Cd_inert);
P_rem_vals = expVals(P_remin_fit,z_fit,P_inert);

CdP_rem_vals = Cd_rem_vals./P_rem_vals;

r = corrcoef(CdP_rem_vals,data.CdP_mmolmol(remi));

CdP_rem_R2 = r(2)^2;


%% remineralization profile plots
if plt.remin_only_prof == 1

    z_pred = 0:opts.max_fit_depth-data.z_euph(1); % depth below euph zone
    z_real = z_pred+data.z_euph(1); % actual depth
    Cd_pred = expVals(Cd_remin_fit,z_pred,Cd_inert);
    P_pred = expVals(P_remin_fit,z_pred,P_inert);
    CdP_pred = Cd_pred./P_pred;

    figure
    sgtitle('Remineralization Fit - No Points Removed','FontSize',24)
    subplot(131)
    plot(data.Cd_pM,data.Depth,'ok','MarkerSize',8,'LineWidth',2)
    set(gca,'YDir','Reverse','XAxisLocation','top','FontSize',16)
    ylim([0,opts.max_fit_depth])
    xlabel('Cd [pM]')
    ylabel('Depth [m]')
    hold on
    yline(data.z_euph(1),':k','LineWidth',1)
    plot(Cd_pred,z_real,'--k','LineWidth',2)

    subplot(132)
    plot(data.P_nM,data.Depth,'ok','MarkerSize',8,'LineWidth',2)
    set(gca,'YDir','Reverse','XAxisLocation','top','FontSize',16)
    ylim([0,opts.max_fit_depth])
    xlabel('P [nM]')
    hold on
    yline(data.z_euph(1),':k','LineWidth',1)
    plot(P_pred,z_real,'--k','LineWidth',2)

    subplot(133)
    plot(data.CdP_mmolmol,data.Depth,'ok','MarkerSize',8,'LineWidth',2)
    set(gca,'YDir','Reverse','XAxisLocation','top','FontSize',16)
    ylim([0,opts.max_fit_depth])
    hold on
    xlabel('Cd:P [mmol:mol]')
    yline(data.z_euph(1),':k','LineWidth',1)
    plot(CdP_pred,z_real,'--k','LineWidth',2)

    clear z_pred z_real Cd_pred P_pred CdP_pred

end


%% Cd accumulation

% only calc for stations below the remin R2 threshold
if CdP_rem_R2 < opts.remin_R2_threshold
    [z_fit, sorti] = sort(z_fit); % confirm depths in order
    Cd_meas = Cd_meas(sorti);
    
    npts = length(Cd_meas);

    % count each removal iteration
    m = 0;

    % roll through each point except first point below euphotic zone 
    for pt = 2:npts

        % remove selected point and subsequent points up to max #
        for nrem = 1:min([npts,opts.max_point_removal])
            
            % must have at least 4 points to constrain fit
            if (pt-1) + nrem > npts - 4
                continue
            end

            m = m+1;

            % depths to remove
            z_rem = z_fit(pt:pt-1+nrem);

            % depths to keep
            [z_keep, zi] = setdiff(z_fit,z_rem);

            % recalc remineralization fit to remaining points
            [loc_fit{m}, loc_gof{m}, loc_out{m}] = fit(z_keep,...
                Cd_meas(zi)-Cd_inert,'exp2','Lower',[0,-0.5,0,-0.5],...
                'Upper',[Cd_max,0,Cd_max,0]);
            loc_z_acc{m} = z_rem;
            loc_R2(m) = loc_gof{m}.rsquare;
            loc_rmse(m) = loc_gof{m}.rmse;

            % check if removal was of high Cd vals or low Cd vals
            conc_remin = expVals(Cd_remin_fit,z_rem,Cd_inert);
            conc_acc = expVals(loc_fit{m},z_rem,Cd_inert);
            loc_conc{m} = Cd_meas(pt:pt-1+nrem) - conc_acc;


            % disqualify if removal was of low vals 
            if sum(loc_conc{m}>0) < nrem 
                loc_rmse(m) = Inf;
            end

        end
        
    end

    % choose best refit option
    [rmse_acc, acci] = min(loc_rmse);
    
    % keep original fit if it has lowest rmse
    if rmse_acc>Cd_remin_gof.rmse
        Cd_acc.fit = Cd_remin_fit;
        Cd_acc.gof = Cd_remin_gof;
        Cd_acc.out = Cd_remin_out;
        Cd_acc.z = NaN;
        Cd_acc.acc_conc = 0;

    % otherwise use best refit    
    else
        Cd_acc.fit = loc_fit{acci};
        Cd_acc.gof = loc_gof{acci};
        Cd_acc.out = loc_out{acci};
        Cd_acc.z = loc_z_acc{acci}+data.z_euph(1);
        Cd_acc.acc_conc = loc_conc{acci};

        % relative improvement over original no removal fit
        Cd_acc.relimp_R2 = (loc_gof{acci}.rsquare - Cd_remin_gof.rsquare)...
                            ./Cd_remin_gof.rsquare;
        Cd_acc.relimp_rmse = (loc_gof{acci}.rmse - Cd_remin_gof.rmse)...
                            ./Cd_remin_gof.rmse;


    end

    % clean up local variables
    clear loc_fit loc_gof loc_out loc_z_acc loc_R2 loc_rmse loc_conf conf_conc
end


%% P accumulation

% only calc for stations below the remin R2 threshold
if CdP_rem_R2 < opts.remin_R2_threshold
    [z_fit, sorti] = sort(z_fit); % confirm depths in order
    P_meas = P_meas(sorti);
    
    npts = length(P_meas);

    m = 0;

    % roll through each point except first point below euphotic zone 
    for pt = 2:npts

        % remove selected point and subsequent points up to max #
        for nrem = 1:min([npts,opts.max_point_removal])
            
            % must have at least 4 points to constrain fit
            if (pt-1) + nrem > npts - 4
                continue
            end

            m = m+1;

            % depths to remove
            z_rem = z_fit(pt:pt-1+nrem);

            % depths to keep
            [z_keep, zi] = setdiff(z_fit,z_rem);

            % recalc remineralization fit to remaining points
            [loc_fit{m}, loc_gof{m}, loc_out{m}] = fit(z_keep,...
                P_meas(zi)-P_inert,'exp2','Lower',[0,-0.5,0,-0.5],...
                'Upper',[P_max,0,P_max,0]);
            loc_z_acc{m} = z_rem;
            loc_R2(m) = loc_gof{m}.rsquare;
            loc_rmse(m) = loc_gof{m}.rmse;

            % check if removal was of high P vals or low P vals
            conc_remin = expVals(P_remin_fit,z_rem,P_inert);
            conc_acc = expVals(loc_fit{m},z_rem,P_inert);
            loc_conc{m} = P_meas(pt:pt-1+nrem) - conc_acc;


            % disqualify if removal was of low values
            if sum(loc_conc{m}>0) < nrem 
                loc_rmse(m) = Inf;
            end

        end
        
    end

    % choose best refit option
    [rmse_acc, acci] = min(loc_rmse);
    
    % keep original fit if it has lowest rmse
    if rmse_acc>P_remin_gof.rmse
        P_acc.fit = P_remin_fit;
        P_acc.gof = P_remin_gof;
        P_acc.out = P_remin_out;
        P_acc.z = NaN;
        P_acc.acc_conc = 0;

    % otherwise use best refit    
    else
        P_acc.fit = loc_fit{acci};
        P_acc.gof = loc_gof{acci};
        P_acc.out = loc_out{acci};
        P_acc.z = loc_z_acc{acci}+data.z_euph(1);
        P_acc.acc_conc = loc_conc{acci};

        % relative improvement over no removal fit
        P_acc.relimp_R2 = (loc_gof{acci}.rsquare - P_remin_gof.rsquare)...
                            ./P_remin_gof.rsquare;
        P_acc.relimp_rmse = (loc_gof{acci}.rmse - P_remin_gof.rmse)...
                            ./P_remin_gof.rmse;


    end

    % clean up local variables
    clear loc_fit loc_gof loc_out loc_z_acc loc_R2 loc_rmse loc_conf conf_conc
end


%% removal profile plots 

if plt.removal_prof == 1

    z_pred = 0:opts.max_fit_depth-data.z_euph(1); % depth below euph zone
    z_real = z_pred+data.z_euph(1); % actual depth
    Cd_pred = expVals(Cd_acc.fit,z_pred,Cd_inert);
    P_pred = expVals(P_acc.fit,z_pred,P_inert);
    CdP_pred = Cd_pred./P_pred;

    
    figure
    sgtitle('Refit - Points Removed','FontSize',24)
    subplot(131)
    plot(data.Cd_pM,data.Depth,'ok','MarkerSize',8,'LineWidth',2)
    set(gca,'YDir','Reverse','XAxisLocation','top','FontSize',16)
    ylim([0,opts.max_fit_depth])
    xlabel('Cd [pM]')
    ylabel('Depth [m]')
    hold on
    yline(data.z_euph(1),':k','LineWidth',1)
    if isnan(Cd_acc.z)
        Cd_mkr = '--k';
    else
        Cd_mkr = 'k';
        [~,ii,~] = intersect(data.Depth,Cd_acc.z);
        plot(data.Cd_pM(ii),Cd_acc.z,'.r','MarkerSize',20)
    end
    plot(Cd_pred,z_real,Cd_mkr,'LineWidth',2)

    subplot(132)
    plot(data.P_nM,data.Depth,'ok','MarkerSize',8,'LineWidth',2)
    set(gca,'YDir','Reverse','XAxisLocation','top','FontSize',16)
    ylim([0,opts.max_fit_depth])
    xlabel('P [nM]')
    hold on
    yline(data.z_euph(1),':k','LineWidth',1)
    if isnan(P_acc.z)
        P_mkr = '--k';
    else
        P_mkr = 'k';
        [~,ii,~] = intersect(data.Depth,P_acc.z);
        plot(data.P_nM(ii),P_acc.z,'.r','MarkerSize',10)
    end
    plot(P_pred,z_real,P_mkr,'LineWidth',2)

    subplot(133)
    plot(data.CdP_mmolmol,data.Depth,'ok','MarkerSize',8,'LineWidth',2)
    set(gca,'YDir','Reverse','XAxisLocation','top','FontSize',16)
    ylim([0,opts.max_fit_depth])
    hold on
    xlabel('Cd:P [mmol:mol]')
    yline(data.z_euph(1),':k','LineWidth',1)
    plot(CdP_pred,z_real,'k','LineWidth',2)

    clear z_pred z_real Cd_pred P_pred CdP_pred
end

%% Cd accumulation

% integrate Cd excess
Cd_xsfrac = excessInt(Cd_acc,data,Cd_inert);

% integrate P excess
P_xsfrac = excessInt(P_acc,data,P_inert);

% apply P excess threshold
if Cd_xsfrac>0 && P_xsfrac<opts.P_acc_threshold
    Cd_acc_ID = 1; % Cd acc identified
else
    Cd_acc_ID = 0; % Cd acc not identified
end

% display values
fprintf('Cd_acc_ID: %i \n',Cd_acc_ID)
fprintf('Cd_xsfrac: %.2f',Cd_xsfrac)

%% functions
function xsfrac = excessInt(acc, data, const)
    if isnan(acc.z)
        xsfrac = 0;
    else
        [~,xsi,~] = intersect(data.Depth,acc.z); % excess point depths
        intz = data.Depth(xsi(1)-1:xsi(end)+1)-data.z_euph(xsi(1)); % integration depth bounds
        total_int = trapz(intz,data.Cd_pM(xsi(1)-1:xsi(end)+1)); % total depth-integrated pCd

        remin_vals = expValsDirect(acc.fit.a,acc.fit.b,acc.fit.c,acc.fit.d,...
            intz(1):intz(end),const); % remin profile
        remin_int = trapz(intz(1):intz(end),remin_vals); % depth-integrated

        xs_int = total_int-remin_int; 

        xsfrac = xs_int/total_int;

    end

end


function conc = expVals(fit_obj, depths, const)
% calculate values of exponential fit equations using fit object

    conc = (fit_obj.a .*exp(fit_obj.b .* depths)) +...
        (fit_obj.c .*exp(fit_obj.d .* depths)) + const;

end


function [conc] = expValsDirect(a, b ,c, d, depths, const)
% calculate values of exponential fit equations using parameters directly

        conc = (a .*exp(b .* depths)) +...
            (c .*exp(d .* depths)) + const;

end