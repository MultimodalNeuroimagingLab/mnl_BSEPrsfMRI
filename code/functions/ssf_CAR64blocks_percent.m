function [signaldata, out] = ssf_CAR64blocks_percent( ...
    signaldata, ttt, good_channels, perc_channels, ...
    car_timeint, notchOpt, block64Opt)
%SSF_CAR64BLOCKS_PERCENT Common-average rereferencing for BSEP/CCEP data.
%
%   [signaldata, out] = ssf_CAR64blocks_percent( ...
%       signaldata, ttt, good_channels, perc_channels, ...
%       car_timeint, notchOpt, block64Opt)
%
% Selects a proportion of low-variance recording channels for common
% average referencing. Channel variance is estimated within a specified
% post-stimulation time interval. The common average can be calculated
% across all channels or independently within 64-channel blocks.
%
% Inputs
%   signaldata   - Channels x time x epochs data array.
%   ttt          - Time vector in seconds.
%   good_channels
%                - Channel indices eligible for inclusion in the common
%                  average. Noisy, bad, and stimulation channels should
%                  be excluded before calling this function.
%   perc_channels
%                - Proportion of the lowest-variance channels to include
%                  in the common average. For example, 0.25 includes the
%                  lowest-variance 25% of eligible channels.
%   car_timeint  - Time interval used to estimate channel variance, in
%                  seconds. Default: [0.015 0.500].
%   notchOpt     - Structure controlling notch filtering used only for
%                  selecting CAR channels:
%                    notchOpt.do    = 0 or 1
%                    notchOpt.freq  = line frequency in Hz
%                    notchOpt.srate = sampling frequency in Hz
%   block64Opt   - 0: calculate the reference across all channels.
%                  1: calculate independent references in 64-channel
%                     blocks.
%
% Outputs
%   signaldata   - Common-average-referenced data,
%                  channels x time x epochs.
%   out          - Structure containing channel groups and CAR-channel
%                  selections when block64Opt = 1.
%
% This function is adapted from the MNL CCEP common-average-reference
% implementation historically used as ccep_CAR64blocks_percent1.
%
% Original attribution:
%   DH and HH, Multimodal Neuroimaging Lab, Mayo Clinic, 2020.
%
% Pipeline adaptation:
%   Maria Guadalupe Yanez Ramos
%   Developed with scientific and technical guidance from Dora Hermes
%   and the Multimodal Neuroimaging Lab (MNL) team.
%   May 2026

if isempty(perc_channels)
    perc_channels = .30;
end

if isempty(car_timeint)
    car_timeint = [0.015 0.500];
end

if ~isempty(notchOpt)
    if notchOpt.do==1 % do notch before getting channels used in car
        disp('notch filtering data to decide which channels go in Common Average')
        signalnotch = zeros(size(signaldata));
        for kk = 1:size(signaldata,1) % channels
            signalnotch(kk,:,:) = ieeg_notch(squeeze(signaldata(kk,:,:)), notchOpt.srate, notchOpt.freq);
        end
    end
else
    notchOpt.do = 0;
end


out = [];

if block64Opt==1

    % create blocks to use for 64 channel groups
    set_nrs = 1:ceil(size(signaldata,1)/64);
    for ss = set_nrs
        set_inds = (set_nrs(ss)*64-63):set_nrs(ss)*64; % 1:64 65:128 etc...
        set_inds = set_inds(set_inds<=max(good_channels)); % make sure to stay below last good channel
        out(ss).channels_set = set_inds;
    end

    % split into blocks and get car channels
    for ss = set_nrs

        % total set of channels that are good & in the set:
        these_channel_nrs = intersect(good_channels,out(ss).channels_set);

        % get data for car in this set
        if notchOpt.do==1 % only use notch filtered channels
            these_data = signalnotch(these_channel_nrs,ttt>car_timeint(1) & ttt<car_timeint(2),:);
        elseif notchOpt.do==0 % use data, no notched data
            these_data = signaldata(these_channel_nrs,ttt>car_timeint(1) & ttt<car_timeint(2),:);
        end

        % concatinate trials (these_data will be channels X time*trials
        these_data_cat = reshape(these_data,size(these_data,1),size(these_data,2)*size(these_data,3));

         % calculate variance for each channel across time
        chan_var1 = var(these_data_cat,[],2);

        % first reference with respect to median channel, then redo
        %sort the  variance in increasing order
        srtd_var=sort(chan_var1);

        %get the median (or first higher value for even lists)
        median_sort_var = srtd_var(find(srtd_var>=median(chan_var1),1));

        %get median index
        median_var_ind=find(chan_var1==median_sort_var);

        % subtract median variance channel for CAR calculation
        these_data_cat = minus(these_data_cat,these_data_cat(median_var_ind,:));




        % calculate variance for each channel across time
        chan_var = var(these_data_cat,[],2);

        % set a threshold for which channels to reject based on variance
        var_th = quantile(chan_var,perc_channels);

        % include channels with smallest variance (out of good channels in set)
        chans_incl = setdiff(1:length(these_channel_nrs),find(chan_var>var_th));

        % original channels numbers to use for car
        out(ss).car_channels = these_channel_nrs(chans_incl);
    end


    % now do CAR
    for ss = set_nrs % run across blocks

        % average across car channels
        out(ss).car_data = mean(signaldata(out(ss).car_channels,:,:),1);

        % subtract from all other data
        signaldata(out(ss).channels_set,:,:) = signaldata(out(ss).channels_set,:,:) - repmat(out(ss).car_data,length(out(ss).channels_set),1,1);

    end

elseif block64Opt==0

    % total set of channels that are good :
    these_channel_nrs = good_channels;

    % get data for good channels
    if notchOpt.do==1 % only use notch filtered channels
        these_data = signalnotch(these_channel_nrs,ttt>car_timeint(1) & ttt<car_timeint(2),:);
    elseif notchOpt.do==0 % use data, no notched data
        these_data = signaldata(these_channel_nrs,ttt>car_timeint(1) & ttt<car_timeint(2),:);
    end

    % concatinate trials (these_data will be channels X time*trials
    these_data_cat = reshape(these_data,size(these_data,1),size(these_data,2)*size(these_data,3));

    % calculate variance for each channel across time
    chan_var1 = var(these_data_cat,[],2);

    % first reference with respect to median channel, then redo
    %sort the  variance in increasing order
    srtd_var=sort(chan_var1);

    %get the median (or first higher value for even lists)
    median_sort_var = srtd_var(find(srtd_var>=median(chan_var1),1));

    %get median index
    median_var_ind=find(chan_var1==median_sort_var);

    % subtract median variance channel for CAR calculation
    these_data_cat = minus(these_data_cat,these_data_cat(median_var_ind,:));

    % subtract median variance channel from actual data as well
    signaldata = minus(signaldata,signaldata(median_var_ind,:,:));

    % calculate variance for each channel across time, again
    chan_var = var(these_data_cat,[],2);

    % set a threshold for which channels to reject based on variance
    var_th = quantile(chan_var,perc_channels);

    % include channels with smallest variance (out of good channels in set)
    chans_incl = setdiff(1:length(these_channel_nrs),find(chan_var>var_th));

    % original channels numbers to use for car
    car_channels = these_channel_nrs(chans_incl);

    % average across car channels
    car_data = mean(signaldata(car_channels,:,:),1);

    % subtract from all other data
    signaldata = signaldata - repmat(car_data,size(signaldata,1),1,1);

end
