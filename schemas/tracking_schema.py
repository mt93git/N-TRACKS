from pydantic import BaseModel, Field, ConfigDict

class TrackingCalibrationParams(BaseModel):
    """Phase 8 Forensic Calibration Constants."""
    model_config = ConfigDict(frozen=True, extra="forbid")

    k_vol: float = Field(default=1.107220, description="3D Volume calibration factor")
    k_vel: float = Field(default=1.288988, description="Velocity calibration factor")
    k_sph: float = Field(default=1.004726, description="Sphericity calibration factor")

class SegmentationThresholds(BaseModel):
    """3D Volumetric Segmentation & Association Parameters."""
    model_config = ConfigDict(frozen=True, extra="forbid")

    segment_threshold_multiplier: float = Field(default=3.5, ge=0.5, le=10.0)
    min_volume: float = Field(default=90.0, ge=10.0, le=500.0)
    max_volume: float = Field(default=1326.0, ge=500.0, le=10000.0)
    min_sphericity: float = Field(default=0.57, ge=0.0, le=1.0)
    min_voxel_count: int = Field(default=50, ge=5)
    cost_of_non_assign: float = Field(default=20.0, ge=1.0)
    default_channel: int = Field(default=1, ge=1, le=4)

class TrackingConfiguration(BaseModel):
    """Unified Tracking & Calibration Configuration."""
    model_config = ConfigDict(frozen=True, extra="forbid")

    k_vol: float = Field(default=1.107220, description="3D Volume calibration factor")
    k_vel: float = Field(default=1.288988, description="Velocity calibration factor")
    k_sph: float = Field(default=1.004726, description="Sphericity calibration factor")
    segment_threshold_multiplier: float = Field(default=3.5, ge=0.5, le=10.0)
    min_volume: float = Field(default=90.0, ge=10.0, le=500.0)
    max_volume: float = Field(default=1326.0, ge=500.0, le=10000.0)
    min_sphericity: float = Field(default=0.57, ge=0.0, le=1.0)
    min_voxel_count: int = Field(default=50, ge=5)
    cost_of_non_assign: float = Field(default=20.0, ge=1.0)
    default_channel: int = Field(default=1, ge=1, le=4)

class FrameByFrameRecord(BaseModel):
    """Instantaneous frame coordinate schema for MASTER_RESULTS_FRAME_BY_FRAME.csv."""
    model_config = ConfigDict(frozen=True, extra="ignore")

    SampleID: str
    Condition: str
    CellID: int
    Frame: int
    SourceFile: str
    Centroid_X: float = Field(description="Raw unscaled voxel coordinate")
    Centroid_Y: float = Field(description="Raw unscaled voxel coordinate")
    Centroid_Z: float = Field(description="Raw unscaled slice coordinate")
    Volume: float = Field(description="Calibrated 3D volume (K_VOL=1.107220 pre-applied in engine)")
    SurfaceArea: float = Field(description="Mesh surface area in um^2")
    Sphericity: float = Field(description="Intermediate sphericity from calibrated volume")
    AxisLength_1: float
    AxisLength_2: float
    AxisLength_3: float
    MeanIntensity: float

