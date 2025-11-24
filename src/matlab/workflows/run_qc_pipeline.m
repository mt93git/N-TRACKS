function run_QC_centric_workflow_v6()
    % =========================================================================
    % MASTER WORKFLOW - QC-Centric Version 6 (Default Group Names)
    % v1.1 - Portable (Self-loading JAR)
    % =========================================================================

    % --- Auto-load local Bio-Formats library ---
    persistent bioformats_loaded_main;
    if isempty(bioformats_loaded_main)
        fprintf('Loading local Bio-Formats library...');
        script_path = fileparts(mfilename('fullpath'));
        jar_path = fullfile(script_path, 'lib', 'bioformats_package.jar');
        if ~isfile(jar_path), error('CRITICAL: bioformats_package.jar not found in %s', fullfile(script_path, 'lib')); end
        if ~ismember(jar_path, javaclasspath('-dynamic')), javaaddpath(jar_path); fprintf(' Done.\n');
        else, fprintf(' Already loaded.\n'); end
        bioformats_loaded_main = true;
    end

    clc; fprintf('===== STARTING ACME QC-CENTRIC WORKFLOW v6 =====\n\n');
    scripts_dir=fileparts(mfilename('fullpath')); addpath(scripts_dir);
    output_dir='C:\MATLAB_Pipeline_Final_Checkpoint\03_output'; % Point to new checkpoint output

    fprintf('Step 1: Please select an ND2 file...\n');
    [fileName,pathName]=uigetfile('*.nd2','Select an ND2 File');
    if isequal(fileName,0),fprintf('   - User cancelled. Terminated.\n');return;end
    selected_nd2_path=fullfile(pathName,fileName); fprintf('   - Selected: %s\n\n',selected_nd2_path);

    fprintf('Step 2: Inspecting file metadata...\n');
    reader=bfGetReader(selected_nd2_path);
    series_count=reader.getSeriesCount();
    channel_count = reader.getSizeC();
    reader.close();
    fprintf('   - Found %d series and %d channels.\n',series_count, channel_count);

    if channel_count >= 2
        channel_to_process = 2; fprintf('   - Automatically selecting Channel 2 for processing.\n');
    else
        channel_to_process = 1; fprintf('   - Only 1 channel found. Automatically selecting Channel 1 for processing.\n');
    end

    default_p=struct('thr',3.5,'minV',90,'maxV',1326,'minS',0.57);
    qc_defaults=struct('threshold_multiplier',default_p.thr,'use_watershed',false,'watershed_sensitivity',0.5);
    series_params_list = repmat(struct('threshold_multiplier', default_p.thr), series_count, 1);

    run_qc = input('Run Interactive QC to refine threshold? (y/n) [n]: ', 's');
    if isempty(run_qc), run_qc = 'n'; end
    if strcmpi(run_qc, 'y')
        fprintf('\n--- Starting Multi-Series/Channel QC Module ---\n');
        all_series_channel_params = ntracks_qc_gui(selected_nd2_path, series_count, channel_count, qc_defaults);
        fprintf('   - QC parameters potentially updated.\n');
        for s = 1:series_count
             channel_to_get = min(2, channel_count);
             tuned_params = all_series_channel_params(s, channel_to_get);
             series_params_list(s).threshold_multiplier = tuned_params.threshold_multiplier;
             fprintf('   - Series %d will use Threshold Multiplier: %.2f (from QC Channel %d)\n', s, series_params_list(s).threshold_multiplier, channel_to_get);
        end
    else
        fprintf('   - Skipping interactive QC. Using default parameters.\n');
    end
    fprintf('\n');

    fprintf('Step 3: Building manifest...\n');
    manifest=table();
    for s=1:series_count
        current_threshold = series_params_list(s).threshold_multiplier;
        new_row={s, selected_nd2_path, s-1, 'standard', current_threshold, default_p.minV, default_p.maxV, default_p.minS, channel_to_process};
        manifest=[manifest;new_row];
    end
    manifest.Properties.VariableNames={'experiment_id','source_nd2_path','series_index','segmentation_method','segment_threshold_multiplier','min_volume','max_volume','min_sphericity', 'channel_index'};
    fprintf('   - Manifest created.\n\n');

    group_map=containers.Map('KeyType','char','ValueType','any');
    fprintf('Step 4: Enter group names manually (Defaults to experiment ID)...\n');
    for i=1:height(manifest)
        [~,nd2_name,~]=fileparts(manifest.source_nd2_path{i});
        desc_id=string(nd2_name)+"_"+manifest.series_index(i);
        prompt=sprintf('Enter group for %s [%s]: ', desc_id, desc_id);
        user_input = input(prompt,'s');
        if isempty(user_input), group_map(char(desc_id)) = desc_id;
        else, group_map(char(desc_id)) = string(user_input); end
    end

    fprintf('\nStep 5: Executing main processing pipeline (Using v9.2 INTEGRATED Engine)...\n');
    all_results={};
    for i=1:height(manifest)
        experiment_info=manifest(i,:);
        processed_data=process_engine_v9_2_INTEGRATED(experiment_info);
        if ~isempty(processed_data)
            num_rows=height(processed_data); [~,nd2_name,~]=fileparts(experiment_info.source_nd2_path{1});
            desc_id=string(nd2_name)+"_"+experiment_info.series_index;
            group_name=group_map(char(desc_id));
            metadata=table(repmat(desc_id,num_rows,1),repmat(string(fileName),num_rows,1),repmat(group_name,num_rows,1),'VariableNames',{'experiment_id','SourceFile','group_clean'});
            all_results{end+1}=[processed_data,metadata];
        end
    end

    if isempty(all_results),fprintf('No data processed.\n');return;end
    master_table=vertcat(all_results{:});

    all_cols = master_table.Properties.VariableNames;
    metadata_cols_shiny = {'CellID','experiment_id', 'SourceFile', 'group_clean'};
    feature_cols_shiny = all_cols(~ismember(all_cols, metadata_cols_shiny));
    final_ordered_table = master_table(:, [{'CellID'}, feature_cols_shiny, {'experiment_id', 'SourceFile', 'group_clean'}]);

    [~,nd2_name,~]=fileparts(selected_nd2_path);
    output_csv_path=fullfile(output_dir,['MASTER_RESULTS_',nd2_name,'_QCCentric_v6_DefaultGroup.csv']);
    writetable(final_ordered_table,output_csv_path);
    fprintf('\n? SUCCESS: Final CSV saved to %s\n',output_csv_path);
    fprintf('===== ACME QC-CENTRIC WORKFLOW v6 COMPLETE =====\n');
end

