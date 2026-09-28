from typing import Optional
from pydantic import BaseModel, Field, ConfigDict

class PatientClinicalRecord(BaseModel):
    """Immutable clinical demographic and in-depth etiology record."""
    model_config = ConfigDict(frozen=True, extra="forbid")

    patient_id: str = Field(description="Patient identifier (e.g., P1, P22)")
    mrn: Optional[str] = Field(default=None, description="Hospital Medical Record Number")
    age: Optional[float] = Field(default=None, ge=0, le=120)
    gender: Optional[str] = Field(default=None)
    global_condition: str = Field(description="Global category: Hospitalized_Control, Infection_Wo_Sepsis, Severe_Sepsis, Septic_Shock, FluA")
    specific_etiology: str = Field(description="Detailed clinical diagnosis/etiology (e.g. C. diff colitis, DKA, Cholangitis)")
    ex_vivo_treatment: str = Field(description="Ex vivo stimulation or control: Baseline, Vehicle_PBS, TL02-59, RvD4")
    raw_file: Optional[str] = Field(default=None, description="Source raw ND2 filename")
    file_status: str = Field(description="Physical file presence on disk: VERIFIED_LOCAL or MISSING_ON_DISK")
    clinical_note: Optional[str] = Field(default=None)

class DonorTreatmentRecord(BaseModel):
    """Immutable healthy donor experimental record."""
    model_config = ConfigDict(frozen=True, extra="forbid")

    sample_id: str = Field(description="Unique run sample ID (e.g., D16_20250930_LPS)")
    donor_id: str = Field(description="Donor number (e.g., D0, D16)")
    date: Optional[str] = Field(default=None)
    treatment_molecule: str = Field(description="Treatment condition: LPS, LTB4, RvD4, PGE2, PMA, TL0259, or Control")
    category: str = Field(description="High-level category: Healthy_Untreated, Healthy_Treated, Healthy_Treated_Recovery")
    source_nd2_file: str = Field(description="Source ND2 filename")
    series_index: int = Field(default=0)
    channel_index: int = Field(default=1)
