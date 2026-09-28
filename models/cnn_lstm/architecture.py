"""
N-TRACKS Spatiotemporal CNN-LSTM Deep Learning Architecture
============================================================
Author: Maxence Tricaud
License: MIT License

Combines 3D Convolutional Neural Networks (3D-CNN) for volumetric spatial feature
extraction with Bidirectional Long Short-Term Memory (Bi-LSTM) networks for temporal
kinetic tracking, velocity forecasting, and Markovian behavioral state transition modeling
in 3D/4D confocal live-cell imaging of motile and resident leukocytes.
"""

import math
from typing import Dict, Tuple, Optional
import torch
import torch.nn as nn
import torch.nn.functional as F


class VolumetricSpatialEncoder(nn.Module):
    """
    3D-CNN Backbone for extracting spatial and morphological features
    from 3D voxel patches (D x H x W) across time frames.
    """

    def __init__(self, in_channels: int = 1, base_filters: int = 32, feature_dim: int = 128):
        super().__init__()
        self.conv1 = nn.Conv3d(in_channels, base_filters, kernel_size=3, padding=1)
        self.bn1 = nn.BatchNorm3d(base_filters)
        
        self.conv2 = nn.Conv3d(base_filters, base_filters * 2, kernel_size=3, padding=1)
        self.bn2 = nn.BatchNorm3d(base_filters * 2)
        
        self.conv3 = nn.Conv3d(base_filters * 2, base_filters * 4, kernel_size=3, padding=1)
        self.bn3 = nn.BatchNorm3d(base_filters * 4)
        
        self.pool = nn.MaxPool3d(kernel_size=2, stride=2)
        self.adaptive_pool = nn.AdaptiveAvgPool3d((2, 2, 2))
        self.proj = nn.Linear(base_filters * 4 * 8, feature_dim)
        self.dropout = nn.Dropout3d(0.15)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        x: Tensor of shape (B, C, D, H, W)
        Returns: Tensor of shape (B, feature_dim)
        """
        out = F.leaky_relu(self.bn1(self.conv1(x)), negative_slope=0.1)
        out = self.pool(out)
        out = F.leaky_relu(self.bn2(self.conv2(out)), negative_slope=0.1)
        out = self.pool(out)
        out = F.leaky_relu(self.bn3(self.conv3(out)), negative_slope=0.1)
        out = self.adaptive_pool(out)
        out = self.dropout(out)
        out = torch.flatten(out, start_dim=1)
        features = self.proj(out)
        return features


class SpatiotemporalCNNLSTM(nn.Module):
    """
    Hybrid 3D-CNN + Bi-LSTM Deep Learning Pipeline for 4D Cell Kinetics.
    
    1. 3D-CNN extracts volumetric shape & texture embeddings per frame t.
    2. Bi-LSTM integrates temporal dynamics across sequence t=1..T.
    3. Multi-task heads:
       - Displacement Regression Head (dx, dy, dz for LAP tracking guidance)
       - Morphokinetic State Classifier Head (Arrested, Motile, Meandering, Swarming, Extravasating)
       - Markov State Transition Probability Matrix Head P_ij
    """

    def __init__(
        self,
        in_channels: int = 1,
        feature_dim: int = 128,
        hidden_dim: int = 256,
        num_lstm_layers: int = 2,
        num_states: int = 5,
        bidirectional: bool = True
    ):
        super().__init__()
        self.num_states = num_states
        self.hidden_dim = hidden_dim
        self.bidirectional = bidirectional
        num_directions = 2 if bidirectional else 1

        # Spatial 3D-CNN encoder
        self.spatial_encoder = VolumetricSpatialEncoder(
            in_channels=in_channels,
            base_filters=32,
            feature_dim=feature_dim
        )

        # Temporal Bi-LSTM
        self.lstm = nn.LSTM(
            input_size=feature_dim,
            hidden_size=hidden_dim,
            num_layers=num_lstm_layers,
            batch_first=True,
            bidirectional=bidirectional,
            dropout=0.2 if num_lstm_layers > 1 else 0.0
        )

        lstm_out_dim = hidden_dim * num_directions

        # Head 1: 3D Velocity / Displacement Forecast (dx, dy, dz)
        self.displacement_head = nn.Sequential(
            nn.Linear(lstm_out_dim, 128),
            nn.LeakyReLU(0.1),
            nn.Linear(128, 3)
        )

        # Head 2: Morphokinetic State Classification
        # (0: Quiescent, 1: Directed-Chemotactic, 2: Meandering, 3: Swarming, 4: Extravasating)
        self.state_head = nn.Sequential(
            nn.Linear(lstm_out_dim, 128),
            nn.LeakyReLU(0.1),
            nn.Linear(128, num_states)
        )

        # Head 3: Markov Transition Matrix Prediction P_ij
        self.transition_head = nn.Sequential(
            nn.Linear(lstm_out_dim, num_states * num_states)
        )

    def forward(self, x_seq: torch.Tensor) -> Dict[str, torch.Tensor]:
        """
        x_seq: Tensor of shape (B, T, C, D, H, W) representing 4D temporal volumetric patches.
        Returns:
            {
                'displacement': (B, T, 3) predicted delta x, y, z
                'state_logits': (B, T, num_states) predicted state logits
                'transition_matrix': (B, num_states, num_states) row-stochastic Markov transition matrix
                'latent_features': (B, T, hidden_dim * num_directions)
            }
        """
        batch_size, seq_len, c, d, h, w = x_seq.shape

        # Fold time into batch for 3D-CNN encoding
        x_flat = x_seq.view(batch_size * seq_len, c, d, h, w)
        spatial_feats = self.spatial_encoder(x_flat)  # (B * T, feature_dim)

        # Unfold sequence for Bi-LSTM
        lstm_input = spatial_feats.view(batch_size, seq_len, -1)
        lstm_out, _ = self.lstm(lstm_input)  # (B, T, hidden_dim * directions)

        # Predict displacements & states per frame
        displacements = self.displacement_head(lstm_out)
        state_logits = self.state_head(lstm_out)

        # Global sequence embedding for Markov transition matrix
        global_repr = torch.mean(lstm_out, dim=1)  # (B, hidden_dim * directions)
        trans_raw = self.transition_head(global_repr).view(batch_size, self.num_states, self.num_states)
        transition_matrix = F.softmax(trans_raw, dim=-1)  # row-stochastic: sum_j P_ij = 1

        return {
            "displacement": displacements,
            "state_logits": state_logits,
            "transition_matrix": transition_matrix,
            "latent_features": lstm_out
        }


class MorphokineticLoss(nn.Module):
    """
    Multi-task loss combining:
    1. Smooth L1 loss for 3D displacement forecasting.
    2. Cross-entropy loss for morphokinetic state prediction.
    3. Markov consistency loss enforcing stationary state preservation.
    """

    def __init__(self, alpha_disp: float = 1.0, beta_state: float = 1.0, gamma_markov: float = 0.5):
        super().__init__()
        self.alpha_disp = alpha_disp
        self.beta_state = beta_state
        self.gamma_markov = gamma_markov
        self.l1_loss = nn.SmoothL1Loss()
        self.ce_loss = nn.CrossEntropyLoss()

    def forward(
        self,
        predictions: Dict[str, torch.Tensor],
        target_displacement: torch.Tensor,
        target_states: torch.Tensor,
        target_transitions: Optional[torch.Tensor] = None
    ) -> Dict[str, torch.Tensor]:
        loss_disp = self.l1_loss(predictions["displacement"], target_displacement)
        
        # Flatten for cross-entropy
        b, t, num_states = predictions["state_logits"].shape
        logits_flat = predictions["state_logits"].view(b * t, num_states)
        targets_flat = target_states.view(b * t)
        loss_state = self.ce_loss(logits_flat, targets_flat)

        loss_markov = torch.tensor(0.0, device=predictions["displacement"].device)
        if target_transitions is not None:
            loss_markov = F.kl_div(
                predictions["transition_matrix"].log(),
                target_transitions,
                reduction="batchmean"
            )

        total_loss = (self.alpha_disp * loss_disp +
                      self.beta_state * loss_state +
                      self.gamma_markov * loss_markov)

        return {
            "total_loss": total_loss,
            "displacement_loss": loss_disp,
            "state_loss": loss_state,
            "markov_loss": loss_markov
        }
