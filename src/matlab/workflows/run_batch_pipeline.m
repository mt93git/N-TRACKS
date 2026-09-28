function run_BATCH_workflow_v2()
    % =========================================================================
    % BATCH WORKFLOW v2 - Includes Robust Error Handling
    % v1.1 - Portable (Self-loading JAR)
    % =========================================================================

    % --- Auto-load local Bio-Formats library ---
    persistent bioformats_loaded_batch;
    if isempty(bioformats_loaded_batch)
        fprintf('Loading Bio-Formats library...');
        script_path = fileparts(mfilename('fullpath'));
        repo_root = fileparts(fileparts(fileparts(script_path)));
        jar_path = fullfile(repo_root, 'lib', 'bioformats_package.jar');
        if ~isfile(jar_path)
            fprintf('\n   Bio-Formats jar not found locally. Downloading from Open Microscopy Environment...\n');
            lib_dir = fullfile(repo_root, 'lib');
            if ~exist(lib_dir, 'dir'), mkdir(lib_dir); end
            try
                websave(jar_path, 'https://downloads.openmicroscopy.org/bio-formats/7.0.0/artifacts/bioformats_package.jar');
                fprintf('   Successfully downloaded Bio-Formats bridge.\n');
            catch ME
                error('CRITICAL: bioformats_package.jar not found in %s. Please run install_ntracks_env.R to download it.', lib_dir);
            end
        end
        if ~ismember(jar_path, javaclasspath('-dynamic')), javaaddpath(jar_path); fprintf(' Done.\n');
        else, fprintf(' Already loaded.\n'); end
        bioformats_loaded_batch = true;
    end

    clc; fprintf('===== STARTING ACME BATCH WORKFLOW v2 =====\n\n');
    scripts_dir = fileparts(mfilename('fullpath'));
    addpath(scripts_dir);
    
    % Dynamic workspace discovery: locate data directories relative to repository root
    repo_root = fileparts(fileparts(fileparts(scripts_dir)));
    output_dir = fullfile(repo_root, 'data', 'output_examples');
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    % Scan raw_nd2 directory dynamically if present
    raw_data_dir = fullfile(repo_root, 'data', 'raw_nd2');
    if exist(raw_data_dir, 'dir')
        nd2_entries = dir(fullfile(raw_data_dir, '*.nd2'));
        nd2_file_list = fullfile(raw_data_dir, {nd2_entries.name});
    else
        % Default fallback list with anonymized research identifiers
        nd2_file_list = { ...
            fullfile(repo_root, 'data', 'demo_Donor_P1_Vehicle.nd2'), ...
            fullfile(repo_root, 'data', 'demo_Donor_P1_LPS.nd2') ...
        };
    end
    fprintf('Found %d ND2 files configured for batch mode.\n\n', numel(nd2_file_list));
    default_p=struct('thr',3.5,'minV',90,'maxV',1326,'minS',0.57);
    master_results_list = {};

    for f_idx = 1:numel(nd2_file_list)
        selected_nd2_path = nd2_file_list{f_idx};
        [~, fileName, ~] = fileparts(selected_nd2_path);
        fprintf('--- Processing File %d of %d: %s ---\n', f_idx, numel(nd2_file_list), fileName);

        try
            if ~isfile(selected_nd2_path)
                fprintf('   ? ERROR: File not found. Skipping.\n\n');
                continue;
            end

            reader=bfGetReader(selected_nd2_path);
            series_count=reader.getSeriesCount();
            channel_count = reader.getSizeC();
            reader.close();

            if channel_count >= 2, channel_to_process = 2; else, channel_to_process = 1; end
            fprintf('   - Found %d series. Using Channel %d for processing.\n', series_count, channel_to_process);

            manifest=table();
            for s=1:series_count
                current_threshold = default_p.thr;
                new_row={s, selected_nd2_path, s-1, 'standard', current_threshold, default_p.minV, default_p.maxV, default_p.minS, channel_to_process};
                manifest=[manifest;new_row];
            end
            manifest.Properties.VariableNames={'experiment_id','source_nd2_path','series_index','segmentation_method','segment_threshold_multiplier','min_volume','max_volume','min_sphericity', 'channel_index'};

            file_results_list = {};
            for i=1:height(manifest)
                experiment_info=manifest(i,:);
                fprintf('   - Processing Series %d...\n', experiment_info.series_index);
                processed_data=ntracks_core_engine(experiment_info);

                if ~isempty(processed_data)
                    num_rows=height(processed_data); [~,nd2_name,~]=fileparts(experiment_info.source_nd2_path{1});
                    desc_id=string(nd2_name)+"_"+experiment_info.series_index; group_name=desc_id;
                    metadata=table(repmat(desc_id,num_rows,1),repmat(string(fileName),num_rows,1),repmat(group_name,num_rows,1),'VariableNames',{'experiment_id','SourceFile','group_clean'});
                    file_results_list{end+1}=[processed_data,metadata];
                end
            end

            if ~isempty(file_results_list)
                 master_results_list{end+1} = vertcat(file_results_list{:});
                 fprintf('   - ? Finished processing %s successfully.\n\n', fileName);
            else
                fprintf('   - ?? No tracks found in %s. Skipping aggregation for this file.\n\n', fileName);
            end
        catch ME
            fprintf('   ? ERROR processing file %s. Skipping this file.\n', fileName);
            fprintf('   Error details: %s\n\n', ME.message);
        end
    end

    if isempty(master_results_list),fprintf('? ERROR: No data processed successfully from any file.\n');return;end

    fprintf('--- Aggregating results from all successfully processed files ---\n');
    master_table=vertcat(master_results_list{:});

    all_cols = master_table.Properties.VariableNames;
    metadata_cols_shiny = {'CellID','experiment_id', 'SourceFile', 'group_clean'};
    feature_cols_shiny = all_cols(~ismember(all_cols, metadata_cols_shiny));
    final_ordered_table = master_table(:, [{'CellID'}, feature_cols_shiny, {'experiment_id', 'SourceFile', 'group_clean'}]);

    output_csv_path=fullfile(output_dir,'MASTER_RESULTS_BATCH_RUN_v2.csv');
    writetable(final_ordered_table,output_csv_path);
    fprintf('\n? SUCCESS: Final Batch CSV saved. Total tracks processed: %d\n', height(final_ordered_table));
    fprintf('   Output Location: %s\n', output_csv_path);
    fprintf('===== ACME BATCH WORKFLOW v2 COMPLETE =====\n');
end

