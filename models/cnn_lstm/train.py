"""
N-TRACKS CNN-LSTM Training Pipeline
===================================
Author: Maxence Tricaud
License: MIT License

Trains the Spatiotemporal CNN-LSTM architecture on 4D microscopy sequences
coupling volumetric spatial feature extraction, displacement regression,
and Markovian state transition supervision.
"""

import os
import argparse
import torch
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import sys
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    from .architecture import SpatiotemporalCNNLSTM, MorphokineticLoss
except (ImportError, ValueError):
    from architecture import SpatiotemporalCNNLSTM, MorphokineticLoss

class Synthetic4DTrackDataset(Dataset):
    """Synthetic dataset simulating 4D leukocyte volumetric time-series for model verification."""

    def __init__(self, num_samples: int = 40, seq_len: int = 8, d: int = 16, h: int = 32, w: int = 32):
        self.num_samples = num_samples
        self.seq_len = seq_len
        self.d = d
        self.h = h
        self.w = w

    def __len__(self) -> int:
        return self.num_samples

    def __getitem__(self, idx: int):
        # Generate simulated 4D patch (T, C=1, D, H, W)
        patches = torch.randn(self.seq_len, 1, self.d, self.h, self.w) * 0.05
        # Simulate moving cell
        target_disp = torch.zeros(self.seq_len, 3)
        target_states = torch.zeros(self.seq_len, dtype=torch.long)
        
        state_id = idx % 5
        for t in range(self.seq_len):
            dx, dy, dz = (t * 0.1, math.sin(t * 0.5) * 0.2, 0.05)
            target_disp[t] = torch.tensor([dx, dy, dz])
            target_states[t] = state_id

        return patches, target_disp, target_states


import math


def train_model(epochs: int = 3, batch_size: int = 2, seq_len: int = 8, lr: float = 1e-3):
    print("Initializing N-TRACKS CNN-LSTM Training Harness...")
    dataset = Synthetic4DTrackDataset(num_samples=20, seq_len=seq_len)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    model = SpatiotemporalCNNLSTM(
        in_channels=1,
        feature_dim=64,
        hidden_dim=128,
        num_lstm_layers=1,
        num_states=5,
        bidirectional=True
    )
    criterion = MorphokineticLoss(alpha_disp=1.0, beta_state=1.0)
    optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)

    model.train()
    for epoch in range(1, epochs + 1):
        total_loss = 0.0
        for patches, disp_targets, state_targets in loader:
            optimizer.zero_grad()
            preds = model(patches)
            losses = criterion(preds, disp_targets, state_targets)
            losses["total_loss"].backward()
            optimizer.step()
            total_loss += losses["total_loss"].item()

        avg_loss = total_loss / len(loader)
        print(f"Epoch [{epoch}/{epochs}] - Loss: {avg_loss:.4f} (Disp Loss: {losses['displacement_loss'].item():.4f}, State Loss: {losses['state_loss'].item():.4f})")

    print("Training loop successfully verified!")
    return model


if __name__ == "__main__":
    train_model(epochs=2, batch_size=2, seq_len=4)
