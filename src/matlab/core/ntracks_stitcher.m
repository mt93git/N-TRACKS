function run_interactive_stitcher()
    % =========================================================================
    % Interactive CSV Stitcher (v3 - with Metadata Editor)
    % =========================================================================
    clc;
    fprintf('===== STARTING INTERACTIVE CSV STITCHER =====\n\n');

    % --- Step 1: File Selection ---
    [fileNames, pathName] = uigetfile('*.csv', 'Select CSV files to stitch', 'MultiSelect', 'on');
    if isequal(fileNames, 0), fprintf('   - User cancelled. Terminated.\n'); return; end
    if ~iscell(fileNames), fileNames = {fileNames}; end
    fprintf('   - Selected %d file(s).\n\n', numel(fileNames));

    % --- Step 2: Stitching ---
    fprintf('Step 2: Loading and combining data...\n');
    all_tables = {};
    for i = 1:numel(fileNames)
        try
            T = readtable(fullfile(pathName, fileNames{i}), 'TextType', 'string');
            all_tables{end+1} = T;
        catch ME, warning('Could not read %s. Skipping.', fileNames{i}); end
    end
    if isempty(all_tables), fprintf('   - No data loaded. Terminated.\n'); return; end
    master_table = vertcat(all_tables{:});
    fprintf('   - Stitching complete. Total tracks: %d\n\n', height(master_table));

    % --- Step 3: Launch Interactive Editor ---
    fprintf('Step 3: Launching metadata editor...\n');
    
    % This function call opens the GUI and waits until it is closed.
    % It returns the fully updated master_table.
    updated_master_table = launch_editor_gui(master_table);

    % Check if the user closed the editor window without applying changes
    if isempty(updated_master_table)
        fprintf('   - Metadata editing cancelled. Terminated.\n');
        return;
    end
    
    fprintf('   - Metadata updates applied.\n\n');

    % --- Step 4: Save Final File ---
    fprintf('Step 4: Please choose a location to save the final stitched file...\n');
    [saveFile, savePath] = uiputfile('*.csv', 'Save Stitched CSV As', 'stitched_master_results.csv');
    if isequal(saveFile, 0), fprintf('   - User cancelled save. Terminated.\n'); return; end
    
    full_save_path = fullfile(savePath, saveFile);
    writetable(updated_master_table, full_save_path);
    
    fprintf('   - Final file saved to: %s\n', full_save_path);
    fprintf('\n===== INTERACTIVE STITCHER COMPLETE =====\n');
end

function modified_table = launch_editor_gui(original_table)
    summary_table = unique(original_table(:, {'experiment_id', 'group_clean'}));
    num_items = height(summary_table);

    fig_height = max(400, 60 + num_items * 30);
    fig = uifigure('Name', 'Interactive Metadata Editor', 'Position', [200 200 800 fig_height]);
    
    % Create a scrollable panel to hold the text boxes
    panel = uipanel(fig, 'Position', [20 70 760 fig_height-80], 'Scrollable', 'on');

    edit_handles = cell(num_items, 2);
    
    uilabel(panel, 'Text', 'Original experiment_id', 'Position', [10, num_items*30, 350, 22], 'FontWeight', 'bold');
    uilabel(panel, 'Text', 'New group_clean', 'Position', [380, num_items*30, 350, 22], 'FontWeight', 'bold');

    for i = 1:num_items
        y_pos = (num_items - i) * 30 + 10;
        
        % Original experiment_id (read-only)
        uilabel(panel, 'Text', summary_table.experiment_id(i), 'Position', [10, y_pos, 350, 22]);
        
        % Editable group_clean
        edit_handles{i, 1} = summary_table.experiment_id(i); % Store original ID for lookup
        edit_handles{i, 2} = uieditfield(panel, 'text', 'Value', summary_table.group_clean(i), 'Position', [380, y_pos, 350, 22]);
    end
    
    uibutton(fig, 'Text', 'Apply Changes & Continue', 'Position', [300, 20, 200, 30], 'ButtonPushedFcn', @applyChanges);

    modified_table = []; % Default return value
    uiwait(fig); % Pause execution until uiresume is called

    function applyChanges(~, ~)
        fprintf('   - Applying changes...\n');
        temp_table = original_table;
        for j = 1:num_items
            original_id = edit_handles{j, 1};
            new_group_name = string(edit_handles{j, 2}.Value);
            
            % Find all rows that match the original experiment_id and update their group_clean
            rows_to_update = (temp_table.experiment_id == original_id);
            temp_table.group_clean(rows_to_update) = new_group_name;
        end
        modified_table = temp_table;
        uiresume(fig);
        delete(fig);
    end
end
