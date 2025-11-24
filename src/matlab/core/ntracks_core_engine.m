function final_track_data = process_engine_v9_2_INTEGRATED(experiment_info)
    % =========================================================================
    % INTEGRATED Scientific Engine (Exact v9.2 Logic + Channel 2 Fixed)
    % v1.1 - Portable (Self-loading JAR)
    % =========================================================================

    % --- Auto-load local Bio-Formats library ---
    persistent bioformats_loaded_engine;
    if isempty(bioformats_loaded_engine)
        fprintf('   (Engine: Loading local Bio-Formats library...)');
        script_path = fileparts(mfilename('fullpath'));
        jar_path = fullfile(script_path, 'lib', 'bioformats_package.jar');
        if ~isfile(jar_path), error('CRITICAL: bioformats_package.jar not found in %s', fullfile(script_path, 'lib')); end
        if ~ismember(jar_path, javaclasspath('-dynamic')), javaaddpath(jar_path); fprintf(' Done.\n');
        else, fprintf(' Already loaded.\n'); end
        bioformats_loaded_engine = true;
    end

    exp_id_str = num2str(experiment_info.experiment_id);
    fprintf('--- Starting v9.2 INTEGRATED Engine for: %s ---\n', exp_id_str);
    nd2_file_path = char(experiment_info.source_nd2_path);

    p = struct();
    p.threshold_multiplier = experiment_info.segment_threshold_multiplier;
    p.min_volume = experiment_info.min_volume;
    p.max_volume = experiment_info.max_volume;
    p.min_sphericity = experiment_info.min_sphericity;
    p.minVoxelCount = 50;
    p.costOfNonAssign = 20;
    p.channel_index = 2; % Hardcoded to Channel 2 as requested

    reader = bfGetReader(nd2_file_path);
    reader.setSeries(experiment_info.series_index);
    nFrames = reader.getSizeT();
    image_depth_z = reader.getSizeZ();

    final_track_data = run_core_processing_v9_2(reader, p, nFrames, image_depth_z);

    reader.close();
    fprintf('Processed %d valid tracks.\n', height(final_track_data));
    fprintf('--- v9.2 INTEGRATED Engine finished for: %s ---\n\n', exp_id_str);
end

% --- Core processing logic REPLICATING v9.2 ---
function final_results_table = run_core_processing_v9_2(reader, p, nFrames, image_depth_z)
    AllCellData=struct('Cells', cell(nFrames,1)); nextAvailableID=1; prevFrameTable=[];
    for t=1:nFrames
        volData=get_volume(reader, t, p.channel_index);
        bw3D=segment_volume_v9_2_simple(volData,p);
        labelVol=bwlabeln(bw3D,26);
        stats3D=regionprops3(labelVol,im2double(volData),'Volume','SurfaceArea','Centroid','PrincipalAxisLength','MeanIntensity');
        if ~isempty(stats3D)
            sphericity=(pi^(1/3)*(6*stats3D.Volume).^(2/3))./stats3D.SurfaceArea;
            volume_filter=stats3D.Volume>=p.min_volume & stats3D.Volume<=p.max_volume;
            sphericity_filter=sphericity>=p.min_sphericity;
            stats3D=stats3D(volume_filter & sphericity_filter,:);
        end
        if isempty(stats3D), AllCellData(t).Cells=[]; continue; end
        nCells=height(stats3D); newFrameTable=zeros(nCells,7);
        for cdx=1:nCells
            volu=stats3D.Volume(cdx);sArea=stats3D.SurfaceArea(cdx);
            if sArea>0,spheric=(pi^(1/3)*((6*volu)^(2/3)))/sArea;else,spheric=0;end
            newFrameTable(cdx,:)=[-1,stats3D.Centroid(cdx,:),volu,sArea,spheric];
        end
        if t>1 && ~isempty(prevFrameTable) && ~isempty(newFrameTable)
            costM=pdist2(prevFrameTable(:,2:4),newFrameTable(:,2:4));
            [assignMap,~,unassDet]=assignDetectionsToTracks(costM,p.costOfNonAssign);
            for a=1:size(assignMap,1),newFrameTable(assignMap(a,2),1)=prevFrameTable(assignMap(a,1),1);end
            for ud=1:numel(unassDet),newFrameTable(unassDet(ud),1)=nextAvailableID;nextAvailableID=nextAvailableID+1;end
            if ~isempty(newFrameTable),maxID=max(newFrameTable(:,1));if maxID>=nextAvailableID,nextAvailableID=maxID+1;end,end
        else
            for iCell=1:nCells,newFrameTable(iCell,1)=nextAvailableID;nextAvailableID=nextAvailableID+1;end
        end
        prevFrameTable=newFrameTable;
        AllCellData(t).Cells=localMakeCellStruct3D_v9_2(newFrameTable,stats3D);
    end
    final_results_table=aggregate_all_features_DEFINITIVE(AllCellData,nFrames,image_depth_z);
end
function vol=get_volume(reader,t,c),nZ=reader.getSizeZ();vol=zeros(reader.getSizeY(),reader.getSizeX(),nZ,'uint16');for z=1:nZ,idx=reader.getIndex(z-1,c-1,t-1)+1;vol(:,:,z)=bfGetPlane(reader,idx);end,end
function bw3D = segment_volume_v9_2_simple(volData, p)
    img_double = double(volData(:)); img_mean = mean(img_double); img_std = std(img_double);
    thresholdVal = img_mean + (p.threshold_multiplier * img_std);
    bw3D = volData > thresholdVal;
    CC = bwconncomp(bw3D, 26); numPixels = cellfun(@numel, CC.PixelIdxList);
    for sdx = 1:numel(CC.PixelIdxList), if numPixels(sdx) < p.minVoxelCount, bw3D(CC.PixelIdxList{sdx}) = false; end, end
end
function cs=localMakeCellStruct3D_v9_2(t,s3D),n=size(t,1);p=struct('ID',[],'Centroid',[],'Volume',[],'SurfaceArea',[],'Sphericity',[],'PrincipalAxisLength',[],'MeanIntensity',[]);if n==0,cs=p;return;end;cs=repmat(p,n,1);for i=1:n,cs(i).ID=t(i,1);cs(i).Centroid=t(i,2:4);cs(i).Volume=t(i,5);cs(i).SurfaceArea=t(i,6);cs(i).Sphericity=t(i,7);cs(i).PrincipalAxisLength=s3D.PrincipalAxisLength(i,:);cs(i).MeanIntensity=s3D.MeanIntensity(i);end,end
function final_results_table = aggregate_all_features_DEFINITIVE(AllCellData, nFrames, image_depth_z)
    all_cell_data_structs = AllCellData(arrayfun(@(s) isstruct(s.Cells), AllCellData));
    if isempty(all_cell_data_structs), uniqueIDs = []; else, uniqueIDs = unique(cell2mat(cellfun(@(c) [c.Cells.ID], num2cell(all_cell_data_structs), 'UniformOutput', false)')); end
    results_list = {};
    for u = 1:numel(uniqueIDs)
        cell_id = uniqueIDs(u); frame_data = table();
        for f = 1:nFrames
            if isempty(AllCellData(f).Cells) || ~isstruct(AllCellData(f).Cells), continue; end
            idx = find([AllCellData(f).Cells.ID] == cell_id);
            if ~isempty(idx), s_frame = AllCellData(f).Cells(idx); s_frame.Frame = f; frame_data = [frame_data; struct2table(s_frame)]; end
        end
        if height(frame_data) < 2, continue; end
        displacements = vecnorm(diff(frame_data.Centroid), 2, 2); s = struct(); s.CellID = cell_id; s.NumTimepoints = height(frame_data);
        s.TraveledDistance = sum(displacements); s.MeanVelocity = mean(displacements); s.StdVelocity = std(displacements); s.MinVelocity = min(displacements); s.MaxVelocity = max(displacements);
        if s.MeanVelocity > 1e-6, s.CVvelocity = s.StdVelocity / s.MeanVelocity; else, s.CVvelocity = 0; end
        linear_disp_vec = frame_data.Centroid(end,:) - frame_data.Centroid(1,:); s.LinearDistance = norm(linear_disp_vec);
        if s.LinearDistance > 1e-6, s.Tortuosity = s.TraveledDistance / s.LinearDistance; else, s.Tortuosity = 1; end
        s.MeanSurfaceArea = mean(frame_data.SurfaceArea, 'omitnan'); s.StdSurfaceArea = std(frame_data.SurfaceArea, 'omitnan'); s.MeanVolume = mean(frame_data.Volume, 'omitnan'); s.StdVolume = std(frame_data.Volume, 'omitnan');
        s.MeanSphericity = mean(frame_data.Sphericity, 'omitnan'); s.StdSphericity = std(frame_data.Sphericity, 'omitnan');
        axes = frame_data.PrincipalAxisLength; s.MeanMajorAxis = mean(axes(:,1), 'omitnan'); s.StdMajorAxis = std(axes(:,1), 'omitnan');
        s.MeanMidAxis = mean(axes(:,2), 'omitnan'); s.StdMidAxis = std(axes(:,2), 'omitnan'); s.MeanMinAxis = mean(axes(:,3), 'omitnan'); s.StdMinAxis = std(axes(:,3), 'omitnan');
        proE = (axes(:,1) - axes(:,2)) ./ axes(:,1); oblE = (axes(:,2) - axes(:,3)) ./ axes(:,2);
        ar1 = axes(:,1) ./ axes(:,2); ar2 = axes(:,1) ./ axes(:,3); ar3 = axes(:,2) ./ axes(:,3);
        s.MeanProlateEllipt = mean(proE, 'omitnan'); s.StdProlateEllipt = std(proE, 'omitnan'); s.MeanOblateEllipt = mean(oblE, 'omitnan'); s.StdOblateEllipt = std(oblE, 'omitnan');
        s.MeanAR1 = mean(ar1, 'omitnan'); s.StdAR1 = std(ar1, 'omitnan'); s.MeanAR2 = mean(ar2, 'omitnan'); s.StdAR2 = std(ar2, 'omitnan'); s.MeanAR3 = mean(ar3, 'omitnan'); s.StdAR3 = std(ar3, 'omitnan');
        z_coords = frame_data.Centroid(:,3); dist_to_surface = min(z_coords, image_depth_z - z_coords);
        s.DistanceToSurface_Mean = mean(dist_to_surface, 'omitnan'); s.DistanceToSurface_Std = std(dist_to_surface, 'omitnan');
        disp_vectors = diff(frame_data.Centroid); angles = atan2d(disp_vectors(:,2), disp_vectors(:,1));
        s.AngleOfDirection_Sum = sum(abs(angles), 'omitnan'); s.AngleOfDirection_Std = std(angles, 'omitnan');
        results_list{end+1} = s;
    end
    if isempty(results_list), final_results_table = table(); else, final_results_table = struct2table(vertcat(results_list{:})); end
end
