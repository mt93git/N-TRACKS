function all_user_params_per_channel = launch_QC_Tool_v_ACME_v4(nd2_filepath, num_series, num_channels, default_params)
    % =========================================================================
    % ACME Interactive QC Tool (v4 - Channel Selection)
    % v1.1 - Portable (Self-loading JAR)
    % =========================================================================

    % --- Auto-load local Bio-Formats library ---
    persistent bioformats_loaded_qc;
    if isempty(bioformats_loaded_qc)
        fprintf('Loading local Bio-Formats library...');
        script_path = fileparts(mfilename('fullpath'));
        jar_path = fullfile(script_path, 'lib', 'bioformats_package.jar');
        if ~isfile(jar_path), error('CRITICAL: bioformats_package.jar not found in %s', fullfile(script_path, 'lib')); end
        if ~ismember(jar_path, javaclasspath('-dynamic')), javaaddpath(jar_path); fprintf(' Done.\n');
        else, fprintf(' Already loaded.\n'); end
        bioformats_loaded_qc = true;
    end

    % --- Data Storage ---
    handles = struct();
    handles.AllSeriesData = cell(num_series, num_channels);
    handles.SeriesParams = repmat(default_params, num_series, num_channels); % [series, channel]
    handles.CurrentSeries = 1;
    handles.CurrentChannel = 1;
    handles.IsPlaying = false;
    handles.Movie = [];
    handles.FrameTimepoints = 0;
    handles.CurrentFrame = 1;
    handles.Reader = bfGetReader(nd2_filepath);
    
    fig = figure('Name', 'Interactive QC Tool v4', 'Position', [100, 100, 1200, 800], 'CloseRequestFcn', @onClosing);
    
    handles.MainAxes = axes('Position', [0.05 0.05 0.6 0.9]);
    axis(handles.MainAxes, 'off');
    
    control_panel = uipanel('Title', 'Controls', 'Position', [0.7 0.05 0.28 0.9]);
    
    % --- Top-level Controls ---
    y_pos = 0.9;
    uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Select Series to Tune:', 'Position', [10, y_pos*100+10, 150, 20], 'HorizontalAlignment', 'left');
    series_names = arrayfun(@(s) sprintf('Series %d', s), 1:num_series, 'UniformOutput', false);
    handles.SeriesDropdown = uicontrol('Parent', control_panel, 'Style', 'popupmenu', 'String', series_names, 'Position', [160, y_pos*100+10, 150, 20], 'Callback', @(s,e) seriesChanged(s,e,fig));
    
    y_pos = y_pos - 0.05;
    uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Select Channel:', 'Position', [10, y_pos*100+10, 150, 20], 'HorizontalAlignment', 'left');
    channel_names = arrayfun(@(c) sprintf('Channel %d', c), 1:num_channels, 'UniformOutput', false);
    handles.ChannelDropdown = uicontrol('Parent', control_panel, 'Style', 'popupmenu', 'String', channel_names, 'Position', [160, y_pos*100+10, 150, 20], 'Callback', @(s,e) channelChanged(s,e,fig));
    
    y_pos = y_pos - 0.08;
    handles.SaveButton = uicontrol('Parent', control_panel, 'Style', 'pushbutton', 'String', 'Save Settings for This Series/Channel', 'Position', [10, y_pos*100, 300, 30], 'Callback', @(s,e) saveSeriesParams(s,e,fig));
    handles.SaveStatus = uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Unsaved changes', 'Position', [10, y_pos*100-20, 300, 20], 'ForegroundColor', 'red');
    
    y_pos = y_pos - 0.1;
    handles.FinishButton = uicontrol('Parent', control_panel, 'Style', 'pushbutton', 'String', 'Finish & Return All Parameters', 'Position', [10, y_pos*100, 300, 40], 'BackgroundColor', [0.8 1 0.8], 'Callback', @(s,e) uiresume(fig));

    % --- Sliders ---
    slider_y_start = y_pos - 0.05;
    slider_gap = 0.07;
    handles.Sliders = struct();
    handles.SliderLabels = struct();
    
    slider_defs = {
        'threshold_multiplier', 'Threshold Multiplier', [0.1, 8.0], '%.1f'
        'use_watershed', 'Use Watershed (1=Y, 0=N)', [0, 1], '%.0f'
        'watershed_sensitivity', 'Watershed Sensitivity', [0.1, 5.0], '%.1f'
    };
    
    for i = 1:length(slider_defs)
        y_pos = slider_y_start - (i-1)*slider_gap;
        name = slider_defs{i}{1};
        label = slider_defs{i}{2};
        range = slider_defs{i}{3};
        fmt = slider_defs{i}{4};
        
        uicontrol('Parent', control_panel, 'Style', 'text', 'String', label, 'Position', [10, y_pos*100+20, 200, 20], 'HorizontalAlignment', 'left');
        handles.SliderLabels.(name) = uicontrol('Parent', control_panel, 'Style', 'text', 'String', sprintf(fmt, default_params.(name)), 'Position', [250, y_pos*100+20, 60, 20], 'HorizontalAlignment', 'right');
        handles.Sliders.(name) = uicontrol('Parent', control_panel, 'Style', 'slider', ...
            'Min', range(1), 'Max', range(2), 'Value', default_params.(name), ...
            'Position', [10, y_pos*100, 300, 20], ...
            'Callback', @(s,e) onSliderChange(s,e,fig,name,fmt));
    end
    
    % --- Playback Controls ---
    y_pos = 0.1;
    handles.PlayButton = uicontrol('Parent', control_panel, 'Style', 'pushbutton', 'String', '? Play', 'Position', [10, y_pos*100, 90, 30], 'Callback', @(s,e) togglePlay(s,e,fig), 'Enable', 'off');
    handles.FrameSlider = uicontrol('Parent', control_panel, 'Style', 'slider', 'Position', [110, y_pos*100+5, 200, 20], 'Min', 1, 'Max', 1, 'Value', 1, 'Callback', @(s,e) onFrameSlider(s,e,fig), 'Enable', 'off');
    handles.FrameLabel = uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Frame 1 / 1', 'Position', [110, y_pos*100-15, 200, 20]);

    % --- Readout Labels ---
    y_pos = 0.02;
    handles.CellCountLabel = uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Detected Cells: 0', 'Position', [10, y_pos*100+10, 150, 20], 'HorizontalAlignment', 'left');
    handles.ContactCountLabel = uicontrol('Parent', control_panel, 'Style', 'text', 'String', 'Resolved Contacts: 0', 'Position', [160, y_pos*100+10, 150, 20], 'HorizontalAlignment', 'left');

    guidata(fig, handles);
    
    % Initial data load for default series/channel
    loadAllSeriesData(fig, nd2_filepath, num_series, num_channels);
    seriesChanged(handles.SeriesDropdown, [], fig); % Trigger initial display
    
    uiwait(fig);
    
    handles = guidata(fig);
    all_user_params_per_channel = handles.SeriesParams;
    delete(fig);
    handles.Reader.close();
end

% --- CALLBACKS AND HELPER FUNCTIONS ---

function loadAllSeriesData(fig, nd2_filepath, num_series, num_channels)
    handles = guidata(fig);
    hWait = waitbar(0, 'Loading all series data. This may take a moment...');
    total_steps = num_series * num_channels;
    step = 0;
    
    for s = 1:num_series
        handles.Reader.setSeries(s-1);
        nFrames = handles.Reader.getSizeT();
        nZ = handles.Reader.getSizeZ();
        
        for c = 1:num_channels
            waitbar(step / total_steps, hWait, sprintf('Loading Series %d/%d, Channel %d/%d...', s, num_series, c, num_channels));
            
            % Create a 2D movie [Y, X, Time] using Max Intensity Projection
            movie = zeros(handles.Reader.getSizeY(), handles.Reader.getSizeX(), nFrames, 'uint16');
            for t = 1:nFrames
                vol = zeros(handles.Reader.getSizeY(), handles.Reader.getSizeX(), nZ, 'uint16');
                for z = 1:nZ
                    idx = handles.Reader.getIndex(z-1, c-1, t-1) + 1;
                    vol(:,:,z) = bfGetPlane(handles.Reader, idx);
                end
                movie(:,:,t) = max(vol, [], 3);
            end
            handles.AllSeriesData{s, c} = movie;
            step = step + 1;
        end
    end
    close(hWait);
    guidata(fig, handles);
    fprintf('   - All series and channel data loaded into memory.\n');
end

function seriesChanged(~, ~, fig)
    handles = guidata(fig);
    handles.CurrentSeries = get(handles.SeriesDropdown, 'Value');
    handles.IsPlaying = false;
    set(handles.PlayButton, 'String', '? Play');
    
    % Update data for the new series
    movie = handles.AllSeriesData{handles.CurrentSeries, handles.CurrentChannel};
    handles.Movie = movie;
    handles.FrameTimepoints = size(movie, 3);
    handles.CurrentFrame = 1;
    
    % Update frame slider
    set(handles.FrameSlider, 'Min', 1, 'Max', handles.FrameTimepoints, 'Value', 1);
    if handles.FrameTimepoints > 1
        set(handles.FrameSlider, 'Enable', 'on', 'SliderStep', [1/(handles.FrameTimepoints-1), 0.1]);
        set(handles.PlayButton, 'Enable', 'on');
    else
        set(handles.FrameSlider, 'Enable', 'off');
        set(handles.PlayButton, 'Enable', 'off');
    end
    
    % Update sliders to match stored params for this series/channel
    current_params = handles.SeriesParams(handles.CurrentSeries, handles.CurrentChannel);
    slider_names = fieldnames(handles.Sliders);
    for i = 1:length(slider_names)
        name = slider_names{i};
        if isfield(current_params, name)
            set(handles.Sliders.(name), 'Value', current_params.(name));
            fmt = get(handles.Sliders.(name), 'UserData'); % Assumes format string stored in UserData (which it isn't, BUG)
            % Let's fix this by just reading the struct
            if strcmp(name, 'use_watershed'), fmt = '%.0f'; else, fmt = '%.1f'; end
            set(handles.SliderLabels.(name), 'String', sprintf(fmt, current_params.(name)));
        end
    end
    
    guidata(fig, handles);
    set(handles.SaveStatus, 'String', 'Loaded settings', 'ForegroundColor', 'black');
    updateSegmentation(fig);
end

function channelChanged(~, ~, fig)
    handles = guidata(fig);
    handles.CurrentChannel = get(handles.ChannelDropdown, 'Value');
    guidata(fig, handles);
    % Just trigger a full series reload, which will pick up the new channel
    seriesChanged(handles.SeriesDropdown, [], fig);
end

function onSliderChange(slider, ~, fig, name, fmt)
    handles = guidata(fig);
    value = get(slider, 'Value');
    set(handles.SliderLabels.(name), 'String', sprintf(fmt, value));
    
    % Store the new value in the parameter struct
    handles.SeriesParams(handles.CurrentSeries, handles.CurrentChannel).(name) = value;
    
    set(handles.SaveStatus, 'String', 'Unsaved changes', 'ForegroundColor', 'red');
    guidata(fig, handles);
    
    if ~handles.IsPlaying
        updateSegmentation(fig);
    end
end

function onFrameSlider(slider, ~, fig)
    handles = guidata(fig);
    handles.CurrentFrame = round(get(slider, 'Value'));
    set(slider, 'Value', handles.CurrentFrame); % Snap to integer
    guidata(fig, handles);
    updateSegmentation(fig);
end

function saveSeriesParams(~, ~, fig)
    handles = guidata(fig);
    % The parameters are already updated on slider change, so we just update the status
    set(handles.SaveStatus, 'String', 'Settings SAVED', 'ForegroundColor', [0 0.6 0]);
    fprintf('Parameters saved for Series %d, Channel %d.\n', handles.CurrentSeries, handles.CurrentChannel);
    guidata(fig, handles); % Save the guidata (though technically already saved)
end

function togglePlay(~, ~, fig)
    handles = guidata(fig);
    handles.IsPlaying = ~handles.IsPlaying;
    if handles.IsPlaying
        set(handles.PlayButton, 'String', '?? Pause');
        guidata(fig, handles);
        while handles.IsPlaying
            handles = guidata(fig);
            next_frame = handles.CurrentFrame + 1;
            if next_frame > handles.FrameTimepoints
                next_frame = 1; % Loop back
            end
            handles.CurrentFrame = next_frame;
            set(handles.FrameSlider, 'Value', handles.CurrentFrame);
            guidata(fig, handles);
            
            updateSegmentation(fig);
            pause(0.1); % Adjust for playback speed
            
            % Check if fig still exists
            if ~isvalid(fig), return; end
            handles = guidata(fig); % Re-fetch handles in case user paused
        end
    else
        set(handles.PlayButton, 'String', '? Play');
        guidata(fig, handles);
    end
end

function updateSegmentation(fig)
    handles = guidata(fig);
    if isempty(handles.Movie), return; end
    
    frame = handles.Movie(:,:,handles.CurrentFrame);
    params = handles.SeriesParams(handles.CurrentSeries, handles.CurrentChannel);
    
    % --- Apply Segmentation ---
    img_double = double(frame(:));
    img_mean = mean(img_double);
    img_std = std(img_double);
    thresholdVal = img_mean + (params.threshold_multiplier * img_std);
    bw = frame > thresholdVal;
    
    contacts = 0;
    if params.use_watershed
        bw_initial_count = length(regionprops(bw, 'Area'));
        D = -bwdist(~bw);
        D = imhmin(D, params.watershed_sensitivity);
        L = watershed(D);
        bw(L == 0) = 0;
        
        CC = bwconncomp(bw);
        final_count = CC.NumObjects;
        contacts = final_count - bw_initial_count;
    end
    
    % --- Display ---
    cla(handles.MainAxes);
    imshow(imadjust(frame), 'Parent', handles.MainAxes, 'Border', 'tight');
    hold(handles.MainAxes, 'on');
    
    stats = regionprops(bw, 'Centroid', 'BoundingBox');
    cell_count = 0;
    
    if ~isempty(stats)
        visboundaries(handles.MainAxes, bw, 'Color', 'r', 'LineWidth', 0.5);
        cell_count = length(stats);
    end
    
    hold(handles.MainAxes, 'off');
    
    % Update labels
    set(handles.FrameLabel, 'String', sprintf('Frame %d / %d', handles.CurrentFrame, handles.FrameTimepoints));
    set(handles.CellCountLabel, 'String', sprintf('Detected Cells: %d', cell_count));
    set(handles.ContactCountLabel, 'String', sprintf('Resolved Contacts: %d', contacts));
end

function onClosing(src, ~)
    handles = guidata(src);
    if ~isempty(handles) && isfield(handles, 'Reader')
        handles.Reader.close();
        fprintf('Bio-Formats reader closed.\n');
    end
    delete(src);
end
