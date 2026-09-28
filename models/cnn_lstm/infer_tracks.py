"""
N-TRACKS CNN-LSTM Inference & Tracking Guidance Runner
======================================================
Author: Maxence Tricaud
License: MIT License

Executes forward inference on 4D time-series microscopy data, predicting:
1. 3D velocity vectors (dx, dy, dz) for Hungarian LAP tracking integration.
2. Phenotypic morphokinetic state sequence per cell.
3. Cohort-level Markov behavioral transition probabilities.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import torch
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    from .architecture import SpatiotemporalCNNLSTM
except (ImportError, ValueError):
    from architecture import SpatiotemporalCNNLSTM

STATE_NAMES = [
    "Stationary_Arrested",
    "Directed_Chemotactic",
    "Meandering_Exploratory",
    "Swarming_Cluster",
    "Extravasating_Active"
]


def generate_synthetic_demo_patch(seq_len: int = 10, d: int = 16, h: int = 32, w: int = 32) -> torch.Tensor:
    """Generates synthetic 4D leukocyte volumetric patch for demonstration and verification."""
    # (B=1, T, C=1, D, H, W)
    patch = torch.randn(1, seq_len, 1, d, h, w) * 0.1
    # Simulate bright cell center
    for t in range(seq_len):
        cz, cy, cx = d // 2, h // 2, w // 2
        patch[0, t, 0, cz-2:cz+3, cy-4:cy+5, cx-4:cx+5] += 2.0 + 0.1 * t
    return patch


def run_inference(
    model: SpatiotemporalCNNLSTM,
    patches: torch.Tensor,
    device: str = "cpu"
) -> dict:
    model.eval()
    model.to(device)
    patches = patches.to(device)

    with torch.no_grad():
        out = model(patches)

    displacements = out["displacement"].cpu().tolist()
    state_logits = out["state_logits"].cpu()
    state_preds = torch.argmax(state_logits, dim=-1).tolist()
    trans_matrix = out["transition_matrix"].cpu().tolist()

    return {
        "displacements": displacements,
        "state_indices": state_preds,
        "state_names": [[STATE_NAMES[idx] for idx in seq] for seq in state_preds],
        "transition_matrix": trans_matrix
    }


def main():
    parser = argparse.ArgumentParser(description="N-TRACKS CNN-LSTM Spatiotemporal Inference")
    parser.add_argument("--seq-len", type=int, default=10, help="Temporal sequence length")
    parser.add_argument("--output-csv", type=str, default="cnn_lstm_predictions.csv", help="Output CSV path")
    args = parser.parse_args()

    print("==================================================================")
    print("      N-TRACKS :: Spatiotemporal CNN-LSTM Inference Engine        ")
    print("==================================================================")
    print(f"Initializing SpatiotemporalCNNLSTM (Sequence Length: {args.seq_len})...")

    model = SpatiotemporalCNNLSTM(
        in_channels=1,
        feature_dim=128,
        hidden_dim=256,
        num_lstm_layers=2,
        num_states=5,
        bidirectional=True
    )

    print("Synthesizing multi-frame 4D volumetric test patch...")
    test_patches = generate_synthetic_demo_patch(seq_len=args.seq_len)

    print("Executing forward pass through 3D-CNN encoder and Bi-LSTM...")
    results = run_inference(model, test_patches)

    print("\n--- Inference Results Summary ---")
    print(f"Predicted Displacement Sequence Length: {len(results['displacements'][0])}")
    print(f"Sample First-Step (dx, dy, dz): {[round(x, 4) for x in results['displacements'][0][0]]}")
    print("\nPredicted Phenotypic State Sequence across frames:")
    for t, st in enumerate(results["state_names"][0]):
        print(f"  Frame {t+1:02d}: {st}")

    print("\nDeep-Estimated Markov Transition Matrix P_ij:")
    df_trans = pd.DataFrame(results["transition_matrix"][0], index=STATE_NAMES, columns=STATE_NAMES)
    print(df_trans.round(3))

    # Save summary dataframe
    pred_records = []
    for t in range(args.seq_len):
        pred_records.append({
            "Frame": t + 1,
            "Predicted_dx": float(results["displacements"][0][t][0]),
            "Predicted_dy": float(results["displacements"][0][t][1]),
            "Predicted_dz": float(results["displacements"][0][t][2]),
            "Predicted_State": results["state_names"][0][t]
        })
    df_out = pd.DataFrame(pred_records)
    df_out.to_csv(args.output_csv, index=False)
    print(f"\nSaved predictions to {args.output_csv}")
    print("==================================================================")


if __name__ == "__main__":
    main()
