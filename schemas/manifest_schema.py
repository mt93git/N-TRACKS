from typing import Optional, List
from pydantic import BaseModel, Field, ConfigDict

class BatchRunEntry(BaseModel):
    """Single run specification within execution manifest."""
    model_config = ConfigDict(frozen=True, extra="forbid")

    filepath: str = Field(description="Relative or absolute path to the raw ND2 file")
    sample_id: str = Field(description="Authoritative sample identifier (e.g. P22_20251013_Septic_Shock_Control)")
    condition: str = Field(description="High-level condition group (Donor or Patient)")
    series_index: int = Field(default=0, ge=0, description="ND2 multiposition series index")
    channel_index: int = Field(default=1, ge=1, description="Acquisition fluorescence channel")
    segment_threshold_multiplier: float = Field(default=3.5)
    min_volume: float = Field(default=90.0)
    manual_check: Optional[str] = Field(default=None, description="Quality control or audit note")

class BatchManifest(BaseModel):
    """Complete manifest collection."""
    model_config = ConfigDict(frozen=True, extra="forbid")

    runs: List[BatchRunEntry]
