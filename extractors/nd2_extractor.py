import os
import sys
try:
    import pandas as pd
except ImportError:
    pd = None
import nd2
import tifffile
import numpy as np
import subprocess
import shutil

# --- 1. Path Configurations (Dynamic & Cluster-Portable) ---
script_dir = os.path.dirname(os.path.abspath(__file__))
# Infer module root: .../src/pipeline/ntracks_unified
module_root = os.path.abspath(os.path.join(script_dir, "../../../.."))
if not os.path.exists(os.path.join(module_root, "configs")):
    # Fallback to current working directory
    module_root = os.getcwd()

project_root = os.environ.get("NTRACKS_ROOT", module_root)

manifest_path = os.environ.get(
    "NTRACKS_MANIFEST",
    os.path.join(project_root, "configs", "batch_manifest_golden.csv")
)
if not os.path.exists(manifest_path):
    manifest_path = os.path.join(project_root, "manifests", "batch_manifest.csv")

raw_nd2_dir = os.environ.get(
    "NTRACKS_RAW_DIR",
    os.path.join(project_root, "User_Downloaded_ND2_Folder")
)
if not os.path.exists(raw_nd2_dir):
    # Portable fallback search paths
    fallback_paths = [
        os.path.join(project_root, "data", "raw_nd2"),
        os.path.join(project_root, "data", "raw"),
        os.path.join(project_root, "raw_nd2")
    ]
    for p in fallback_paths:
        if os.path.exists(p):
            raw_nd2_dir = p
            break

default_temp_tiff = os.path.join(project_root, "data", "temp_tiff")
if os.environ.get("SLURM_TMPDIR"):
    default_temp_tiff = os.path.join(os.environ["SLURM_TMPDIR"], f"ntracks_tiff_{os.environ.get('SLURM_JOB_ID', 'job')}")
elif sys.platform.startswith("linux") and os.path.exists("/tmp"):
    default_temp_tiff = f"/tmp/ntracks_tiff_{os.getpid()}"

temp_tiff_root = os.environ.get("NTRACKS_TEMP_TIFF", default_temp_tiff)
staging_dir = os.path.join(project_root, "output", "staging")
output_dir = os.path.join(project_root, "output")

matlab_path = os.environ.get(
    "MATLAB_BIN",
    shutil.which("matlab") or "/Applications/MATLAB_R2024b.app/bin/matlab"
)
rscript_path = os.environ.get(
    "RSCRIPT_BIN",
    shutil.which("Rscript") or "/usr/local/bin/Rscript"
)

os.makedirs(temp_tiff_root, exist_ok=True)
os.makedirs(staging_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

print("==========================================================")
print("N-TRACKS: COHORT-WIDE BATCH PROCESSING AND TRANSITION ANALYSIS")
print("==========================================================")
print(f"Project Root : {project_root}")
print(f"Raw Reservoir: {raw_nd2_dir}")
print(f"MATLAB Exec  : {matlab_path}")
print(f"Rscript Exec : {rscript_path}")

# --- 2. Load Manifest ---
if not os.path.exists(manifest_path):
    print(f"✖ ERROR: Manifest file not found at: {manifest_path}")
    sys.exit(1)

manifest_df = pd.read_csv(manifest_path)
print(f"Loaded manifest: {manifest_path} ({len(manifest_df)} runs configured)")

# --- 3. Multi-Channel ND2 to TIFF Conversion Helper ---
def convert_nd2_series_to_tiff(nd2_path, series_idx, output_dir, channel_idx=1):
    """
    Extracts a specific series and channel from a multi-dimensional ND2 file
    into a calibrated 3D TIFF stack per frame (ImageJ compatible).
    Bypasses OS-specific Bio-Formats memory constraints and handles arbitrary channel order.
    """
    os.makedirs(output_dir, exist_ok=True)
    with nd2.ND2File(nd2_path) as f:
        sizes = f.sizes
        axes = list(f.sizes.keys())
        data = f.to_dask()
        axis_map = {axis: idx for idx, axis in enumerate(axes)}
        
        nT = sizes.get('T', 1)
        nZ = sizes.get('Z', 1)
        nP = sizes.get('P', 1)
        nC = sizes.get('C', 1)
        
        # 1-based channel_idx converted safely to 0-based index
        c_target = max(0, min(int(channel_idx) - 1, nC - 1))
        
        for t in range(nT):
            sys.stdout.write(f"\r   > Extracting frame {t+1}/{nT} (Series {series_idx}, Chan {channel_idx} [idx {c_target}])...")
            sys.stdout.flush()
            
            index_list = [slice(None)] * len(axes)
            if 'T' in axis_map:
                index_list[axis_map['T']] = t
            if 'P' in axis_map:
                index_list[axis_map['P']] = min(series_idx, nP - 1)
            if 'C' in axis_map:
                index_list[axis_map['C']] = c_target
                
            # Load target volume slice from Dask
            volume = data[tuple(index_list)].compute()
            
            remaining_axes = [ax for ax in axes if ax not in ('T', 'P', 'C')]
            rem_axis_map = {ax: idx for idx, axis in enumerate(remaining_axes)}
            
            z_idx = rem_axis_map.get('Z', 0)
            y_idx = rem_axis_map.get('Y', 1)
            x_idx = rem_axis_map.get('X', 2)
            
            if volume.ndim == 3:
                volume_zyx = np.transpose(volume, (z_idx, y_idx, x_idx))
            elif volume.ndim == 2:
                volume_zyx = np.expand_dims(volume, axis=0)
            else:
                volume_squeezed = np.squeeze(volume)
                if volume_squeezed.ndim == 3:
                    volume_zyx = volume_squeezed
                else:
                    volume_zyx = volume
            
            # Save specifically indexed channel TIFF so MATLAB get_volume_tiff reads exact channel
            tiff_filename = f"frame_{t+1:03d}_chan_{channel_idx:03d}.tif"
            tifffile.imwrite(os.path.join(output_dir, tiff_filename), volume_zyx, imagej=True)
        print("")

# --- 4. Iterative File Processing ---
for idx, row in manifest_df.iterrows():
    nd2_name = os.path.basename(row['Filepath'])
    nd2_filepath = os.path.join(raw_nd2_dir, nd2_name)
    
    # Check recursive search in raw_nd2_dir if not in immediate root
    if not os.path.exists(nd2_filepath):
        found = False
        for r, _, files in os.walk(raw_nd2_dir):
            if nd2_name in files:
                nd2_filepath = os.path.join(r, nd2_name)
                found = True
                break
        if not found:
            print(f"Skipping run {idx+1}/{len(manifest_df)}: file {nd2_name} not found.")
            continue
        
    sample_id = row['SampleID']
    condition = row['Condition']
    series_idx = int(row['series_index'])
    chan_idx = int(row['channel_index']) if not pd.isna(row['channel_index']) else 1
    segment_thr = float(row['segment_threshold_multiplier']) if not pd.isna(row['segment_threshold_multiplier']) else 3.5
    min_vol = float(row['min_volume']) if not pd.isna(row['min_volume']) else 90.0
    
    # Generate signature for restartability
    run_sig = f"{sample_id}_series{series_idx}"
    dest_raw_csv = os.path.join(staging_dir, f"raw_{run_sig}.csv")
    dest_frames_csv = os.path.join(staging_dir, f"frames_{run_sig}.csv")
    
    # Skip if already processed in a previous attempt
    if os.path.exists(dest_raw_csv) and os.path.exists(dest_frames_csv):
        print(f"-> Run {idx+1}/{len(manifest_df)} [{run_sig}] already processed. Skipping.")
        continue
        
    print(f"\nProcessing run {idx+1}/{len(manifest_df)}: {nd2_name} (Series: {series_idx}, Channel: {chan_idx})")
    
    # Extract to temporary TIFF folder
    temp_tiff_dir = os.path.join(temp_tiff_root, f"run_{run_sig}")
    print(f"   > Extracting ND2 series to temp TIFF folder...")
    try:
        convert_nd2_series_to_tiff(nd2_filepath, series_idx, temp_tiff_dir, channel_idx=chan_idx)
    except Exception as e:
        print(f"   ✖ ERROR converting ND2 to TIFF: {e}")
        continue
        
    # Write temporary MATLAB execution script
    m_script_name = f"run_temp_{run_sig}"
    m_script_path = os.path.join(project_root, f"{m_script_name}.m")
    with open(m_script_path, "w") as f:
        f.write(f"""
cd '{project_root}'
addpath(genpath('src'))
test_info = struct();
test_info.source_nd2_path = '{temp_tiff_dir}';
test_info.series_index = 0; % Extracted single-series TIFF directory
test_info.channel_index = {chan_idx};
test_info.segment_threshold_multiplier = {segment_thr};
test_info.min_volume = {min_vol};
test_info.max_volume = 1326;
test_info.min_sphericity = 0.57;

try
    [tracks, frames] = ntracks_core_engine(test_info);
    if ~isempty(tracks)
        writetable(tracks, '{dest_raw_csv}');
        writetable(frames, '{dest_frames_csv}');
        disp('TRACKING_SUCCESS');
    else
        disp('TRACKING_EMPTY');
    end
catch ME
    disp(['TRACKING_ERROR: ', ME.message]);
end
exit;
""")

    # Execute MATLAB headless batch command
    print(f"   > Running MATLAB core tracking engine...")
    try:
        result = subprocess.run(
            [matlab_path, "-batch", m_script_name],
            capture_output=True, text=True, check=True
        )
        output_log = result.stdout + result.stderr
        if "TRACKING_SUCCESS" in output_log:
            print(f"   ✔ SUCCESS: Tracking completed successfully for {run_sig}.")
        elif "TRACKING_EMPTY" in output_log:
            print(f"   > WARNING: Tracking yielded 0 cells for {run_sig}.")
        else:
            print(f"   ✖ ERROR in MATLAB tracking run:\n{output_log}")
    except Exception as e:
        print(f"   ✖ FATAL calling MATLAB batch execution: {e}")
        
    # Clean up temporary folders immediately
    if os.path.exists(temp_tiff_dir):
        shutil.rmtree(temp_tiff_dir)
    if os.path.exists(m_script_path):
        os.remove(m_script_path)

# --- 5. Cohort Stitching and Assembly ---
print("\nStitching all staging outputs together...")
raw_files = [f for f in os.listdir(staging_dir) if f.startswith("raw_") and f.endswith(".csv")]
frames_files = [f for f in os.listdir(staging_dir) if f.startswith("frames_") and f.endswith(".csv")]

meta_lookup = {}
for idx, row in manifest_df.iterrows():
    nd2_name = os.path.basename(row['Filepath'])
    sample_id = row['SampleID']
    condition = row['Condition']
    series_idx = int(row['series_index'])
    run_sig = f"{sample_id}_series{series_idx}"
    meta_lookup[run_sig] = {
        'SampleID': sample_id,
        'Condition': condition,
        'SourceFile': nd2_name
    }

# Stitch Raw cell-averaged trajectories
stitched_raw_dfs = []
for f in raw_files:
    run_sig = f[4:-4]
    meta = meta_lookup.get(run_sig)
    if not meta:
        continue
    df = pd.read_csv(os.path.join(staging_dir, f))
    df['SampleID'] = meta['SampleID']
    df['Condition'] = meta['Condition']
    df['SourceFile'] = meta['SourceFile']
    
    if 'ID' in df.columns:
        df = df.rename(columns={'ID': 'CellID'})
    stitched_raw_dfs.append(df)

if stitched_raw_dfs:
    final_raw = pd.concat(stitched_raw_dfs, ignore_index=True)
    priority = ['SampleID', 'Condition', 'CellID', 'SourceFile']
    others = [c for c in final_raw.columns if c not in priority]
    final_raw = final_raw[priority + others]
    final_raw.to_csv(os.path.join(output_dir, "MASTER_RESULTS_RAW.csv"), index=False)
    print(f"✔ Saved cohort averages to output/MASTER_RESULTS_RAW.csv ({len(final_raw)} cells)")

# Stitch Frame-by-frame cell coordinates
stitched_frame_dfs = []
for f in frames_files:
    run_sig = f[7:-4]
    meta = meta_lookup.get(run_sig)
    if not meta:
        continue
    df = pd.read_csv(os.path.join(staging_dir, f))
    df['SampleID'] = meta['SampleID']
    df['Condition'] = meta['Condition']
    df['SourceFile'] = meta['SourceFile']
    
    if 'ID' in df.columns:
        df = df.rename(columns={'ID': 'CellID'})
    stitched_frame_dfs.append(df)

if stitched_frame_dfs:
    final_frames = pd.concat(stitched_frame_dfs, ignore_index=True)
    priority = ['SampleID', 'Condition', 'CellID', 'Frame', 'SourceFile']
    others = [c for c in final_frames.columns if c not in priority]
    final_frames = final_frames[priority + others]
    final_frames.to_csv(os.path.join(output_dir, "MASTER_RESULTS_FRAME_BY_FRAME.csv"), index=False)
    print(f"✔ Saved cohort frame-by-frame tracks to output/MASTER_RESULTS_FRAME_BY_FRAME.csv ({len(final_frames)} coordinates)")

# --- 6. Run Category-Wise Markov Transition Analysis ---
print("\nExecuting R Markov behavior transition analysis...")
r_script_candidates = [
    os.path.join(project_root, "src/ntracks/pipeline2_mac_intel/r_markov/markov_transition_cohort_analysis.R"),
    os.path.join(project_root, "src/r_integration/markov_transition_cohort_analysis.R"),
    os.path.join(project_root, "r_markov/markov_transition_cohort_analysis.R")
]
r_script_path = None
for candidate in r_script_candidates:
    if os.path.exists(candidate):
        r_script_path = candidate
        break

if r_script_path:
    try:
        r_result = subprocess.run(
            [rscript_path, "--no-init-file", r_script_path],
            capture_output=True, text=True, check=True
        )
        print(r_result.stdout)
        print("✔ SUCCESS: Category transition matrices and plots generated.")
    except Exception as e:
        print(f"✖ ERROR executing category transition analysis: {e}")
        if hasattr(e, 'stderr') and e.stderr:
            print(e.stderr)
else:
    print("ℹ R Markov script not located in standard search paths. Skipping R step.")

print("\nCohort processing completed successfully.")
print("==========================================================")
