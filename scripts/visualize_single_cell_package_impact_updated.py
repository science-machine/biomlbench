#!/usr/bin/env python3
"""
Visualize Single Cell Omics performance comparing agents with and without pre-installed scanpy.

Updated version that works with the combined v1/v2 analysis approach.
The "with scanpy" data comes from the combined analysis (which uses v2 for single cell).
The "without scanpy" data needs to be from a v1-only analysis.

Usage:
    python scripts/visualize_single_cell_package_impact_updated.py \
        --v1-only-dir analysis_results_v1_only \
        --combined-dir analysis_results_v1v2_combined \
        --output-dir single_cell_figures_updated
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple
import warnings

warnings.filterwarnings('ignore')

# Set publication-ready style
plt.style.use('seaborn-v0_8-whitegrid')

# Use Set2 palette colors in the exact order of AGENT_ORDER from original
# Order in original: ['dummy', 'mlagentbench', 'aide', 'stella', 'biomni']
BASE_COLORS = {
    'mlagentbench': '#fc8d62', # Set2[1] - Orange
    'aide': '#8da0cb',        # Set2[2] - Blue  
    'stella': '#e78ac3',      # Set2[3] - Pink
    'biomni': '#a6d854'       # Set2[4] - Green
}

# Agent display names - Simple labels since we have group headers
AGENT_DISPLAY_NAMES = {
    'mlagentbench': '−scanpy',
    'mlagentbench_scanpy': '+scanpy',
    'aide': '−scanpy',
    'aide_scanpy': '+scanpy',
    'stella': '−scanpy',
    'stella_scanpy': '+scanpy',
    'biomni': '−scanpy',
    'biomni_scanpy': '+scanpy'
}

# Agent ordering - Updated: scanpy version first, then without
AGENT_ORDER = [
    'mlagentbench_scanpy', 'mlagentbench',
    'aide_scanpy', 'aide',
    'stella_scanpy', 'stella',
    'biomni_scanpy', 'biomni'
]

def log(message: str):
    """Simple logging with timestamp."""
    import time
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def calculate_sem(values):
    """Calculate standard error of the mean."""
    if len(values) == 0:
        return 0
    return np.std(values, ddof=1) / np.sqrt(len(values))

def load_and_prepare_data(v1_only_dir: Path, combined_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load data from v1-only analysis and combined analysis."""
    
    log("Loading data from v1-only analysis (without scanpy)...")
    v1_replicates = pd.read_csv(v1_only_dir / 'all_replicates.csv')
    
    log("Loading data from combined analysis (v2 for single cell = with scanpy)...")
    combined_replicates = pd.read_csv(combined_dir / 'all_replicates.csv')
    
    # Filter to Single Cell Omics tasks only
    v1_single_cell = v1_replicates[v1_replicates['task_id'].str.startswith('manual/')].copy()
    v2_single_cell = combined_replicates[combined_replicates['task_id'].str.startswith('manual/')].copy()
    
    # Remove dummy agent from both
    v1_single_cell = v1_single_cell[v1_single_cell['agent'] != 'dummy']
    v2_single_cell = v2_single_cell[v2_single_cell['agent'] != 'dummy']
    
    log(f"V1 (without scanpy): {len(v1_single_cell)} Single Cell Omics replicates")
    log(f"V2 (with scanpy): {len(v2_single_cell)} Single Cell Omics replicates")
    
    return v1_single_cell, v2_single_cell

def calculate_statistics(v1_data: pd.DataFrame, v2_data: pd.DataFrame) -> pd.DataFrame:
    """Calculate statistics for plotting."""
    
    results = []
    
    # Get unique agents (excluding dummy)
    agents = sorted(set(v1_data['agent'].unique()) | set(v2_data['agent'].unique()))
    agents = [a for a in agents if a != 'dummy']
    
    # Process each agent
    for agent in agents:
        # V1 data (without scanpy)
        v1_agent_data = v1_data[v1_data['agent'] == agent]
        if len(v1_agent_data) > 0:
            # Calculate task-level means first
            task_means = v1_agent_data.groupby('task_id')['leaderboard_percentile'].mean()
            
            # Calculate overall statistics
            mean_percentile = task_means.mean()
            sem_percentile = calculate_sem(task_means)
            
            results.append({
                'agent': agent,
                'agent_label': agent,
                'has_scanpy': False,
                'mean_percentile': mean_percentile,
                'sem_percentile': sem_percentile,
                'n_tasks': len(task_means),
                'n_replicates': len(v1_agent_data)
            })
        
        # V2 data (with scanpy)
        v2_agent_data = v2_data[v2_data['agent'] == agent]
        if len(v2_agent_data) > 0:
            # Calculate task-level means first
            task_means = v2_agent_data.groupby('task_id')['leaderboard_percentile'].mean()
            
            # Calculate overall statistics
            mean_percentile = task_means.mean()
            sem_percentile = calculate_sem(task_means)
            
            results.append({
                'agent': agent,
                'agent_label': f"{agent}_scanpy",
                'has_scanpy': True,
                'mean_percentile': mean_percentile,
                'sem_percentile': sem_percentile,
                'n_tasks': len(task_means),
                'n_replicates': len(v2_agent_data)
            })
    
    return pd.DataFrame(results)

def create_single_cell_performance_plot_with_points(v1_data: pd.DataFrame, v2_data: pd.DataFrame, 
                                                   stats_df: pd.DataFrame, output_dir: Path):
    """Create the plot with individual task points overlaid."""
    
    log("Creating Single Cell Omics performance plot with task-level points...")
    
    # Create figure - match visualize_results.py sizing
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Sort data by agent order
    plot_data = []
    for agent_label in AGENT_ORDER:
        row = stats_df[stats_df['agent_label'] == agent_label]
        if len(row) > 0:
            plot_data.append(row.iloc[0])
    
    plot_df = pd.DataFrame(plot_data)
    
    # Prepare bar positions with grouped spacing
    # Each group has 2 bars (with/without scanpy) with small gap within group, larger gap between groups
    n_groups = 4  # MLAgentBench, AIDE, STELLA, Biomni
    within_group_spacing = 0.2  # Small gap within each agent group
    between_group_spacing = 0.8  # Smaller gap between agent groups
    width = 0.4  # Slightly wider bars
    
    # Calculate x positions for each bar
    x_positions = []
    for group_idx in range(n_groups):
        base_x = group_idx * (2 * width + within_group_spacing + between_group_spacing)
        x_positions.extend([base_x, base_x + width + within_group_spacing])
    
    x = np.array(x_positions)
    
    # Create bars and overlay points
    for i, (idx, row) in enumerate(plot_df.iterrows()):
        base_agent = row['agent'].lower()
        color = BASE_COLORS[base_agent]
        
        # Get task-level data for points
        if row['has_scanpy']:
            agent_task_data = v2_data[v2_data['agent'] == row['agent']].groupby('task_id')['leaderboard_percentile'].mean()
        else:
            agent_task_data = v1_data[v1_data['agent'] == row['agent']].groupby('task_id')['leaderboard_percentile'].mean()
        
        # Create bar - swap hatching: without scanpy gets the texture
        if row['has_scanpy']:
            bar = ax.bar(x[i], row['mean_percentile'], width,
                         yerr=row['sem_percentile'],
                         color=color, alpha=0.8,
                         edgecolor='black', linewidth=2,
                         capsize=5, error_kw={'linewidth': 2})
        else:
            bar = ax.bar(x[i], row['mean_percentile'], width, 
                         yerr=row['sem_percentile'],
                         color=color, alpha=0.6,
                         edgecolor='black', linewidth=2,
                         hatch='///',
                         capsize=5, error_kw={'linewidth': 2})
        
        # Overlay individual task points
        if len(agent_task_data) > 0:
            # Add jitter to x position for better visibility
            jitter = np.random.normal(0, 0.08, len(agent_task_data))
            x_positions_scatter = np.full(len(agent_task_data), x[i]) + jitter
            
            ax.scatter(x_positions_scatter, agent_task_data.values, 
                      color='black', s=25, alpha=0.7, 
                      edgecolors='white', linewidth=1, zorder=10)
    
    # Customize plot - match visualize_results.py styling
    ax.set_ylabel('Leaderboard Percentile', fontsize=15)
    ax.set_title('Single Cell Omics: Impact of Pre-installed scanpy', fontsize=20)
    
    # Set x-axis with multi-level structure
    # Primary labels (with/without scanpy)
    labels = [AGENT_DISPLAY_NAMES[row['agent_label']] for _, row in plot_df.iterrows()]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='center', fontsize=12)
    
    # Add agent group labels below
    agent_groups = ['MLAgentBench', 'AIDE', 'STELLA', 'Biomni']
    
    # Calculate center positions for each group
    group_positions = []
    for group_idx in range(n_groups):
        # Get the x positions of the two bars in this group
        bar1_x = x[group_idx * 2]
        bar2_x = x[group_idx * 2 + 1]
        # Center position is the midpoint between the centers of the two bars
        center_x = (bar1_x + bar2_x) / 2
        group_positions.append(center_x)
    
    # Add vertical lines to separate agent groups (between groups)
    for i in range(1, len(agent_groups)):
        separator_x = (group_positions[i-1] + group_positions[i]) / 2
        ax.axvline(x=separator_x, color='gray', linestyle='--', alpha=0.3)
    
    # Add agent group labels below the bars - moved down to compensate for rotation
    for group, pos in zip(agent_groups, group_positions):
        ax.text(pos, -0.20, group, transform=ax.get_xaxis_transform(),
                ha='center', va='top', fontsize=13)
    
    # Add "Agent" label below the group labels - moved down to compensate
    center_x = (group_positions[0] + group_positions[-1]) / 2

    
    # Set y-axis properties
    ax.set_ylim(0, 100)
    ax.tick_params(axis='y', labelsize=13)
    ax.tick_params(axis='x', labelsize=13)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add legend - updated to match new scheme
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    legend_elements = [
        Patch(facecolor='gray', edgecolor='black', label='With scanpy', alpha=0.8, linewidth=2),
        Patch(facecolor='gray', edgecolor='black', hatch='///', label='Without scanpy (−scanpy)', alpha=0.6, linewidth=2),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='black', 
               markersize=8, label='Individual tasks', markeredgecolor='white', markeredgewidth=1)
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=12, title_fontsize=13)
    
    # Set x-axis limits to center the groups nicely
    ax.set_xlim(-0.5, max(x) + 0.5)
    
    # Adjust layout to accommodate multi-level x-axis and rotated labels
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.22)  # Make room for the rotated labels and agent group labels
    plt.savefig(output_dir / 'single_cell_omics_scanpy_comparison_with_points.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    log(f"Saved plot to {output_dir / 'single_cell_omics_scanpy_comparison_with_points.png'}")

def main():
    parser = argparse.ArgumentParser(
        description="Visualize Single Cell Omics performance with scanpy comparison"
    )
    parser.add_argument(
        "--v1-only-dir",
        required=True,
        help="Directory for v1-only analysis results (without scanpy)"
    )
    parser.add_argument(
        "--combined-dir",
        default="analysis_results_v1v2_combined",
        help="Directory for combined v1/v2 results (with scanpy for single cell)"
    )
    parser.add_argument(
        "--output-dir",
        default="single_cell_figures_updated",
        help="Directory to save outputs"
    )
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Check if v1-only directory exists
    v1_only_dir = Path(args.v1_only_dir)
    if not v1_only_dir.exists():
        log(f"ERROR: v1-only directory not found: {v1_only_dir}")
        log("\nTo create v1-only results, run:")
        log(f"  python scripts/analyze_results.py --output-dir {args.v1_only_dir} --prefix v1/artifacts")
        return 1
    
    # Load and prepare data
    v1_data, v2_data = load_and_prepare_data(v1_only_dir, Path(args.combined_dir))
    
    # Calculate statistics
    log("Calculating statistics...")
    stats_df = calculate_statistics(v1_data, v2_data)
    
    # Save statistics
    stats_df.to_csv(output_dir / 'single_cell_statistics.csv', index=False)
    log(f"Saved statistics to {output_dir / 'single_cell_statistics.csv'}")
    
    # Create plots
    create_single_cell_performance_plot_with_points(v1_data, v2_data, stats_df, output_dir)
    
    # Print summary
    log("\n" + "="*60)
    log("SINGLE CELL OMICS ANALYSIS COMPLETE")
    log("="*60)
    
    log("\nPerformance Summary:")
    for _, row in stats_df.iterrows():
        agent_name = AGENT_DISPLAY_NAMES[row['agent_label']]
        log(f"  {agent_name}: {row['mean_percentile']:.1f} ± {row['sem_percentile']:.1f}")
    
    # Print comparison
    log("\nImpact of pre-installed scanpy:")
    for agent in ['mlagentbench', 'aide', 'stella', 'biomni']:
        without = stats_df[(stats_df['agent'] == agent) & (~stats_df['has_scanpy'])]
        with_scanpy = stats_df[(stats_df['agent'] == agent) & (stats_df['has_scanpy'])]
        
        if len(without) > 0 and len(with_scanpy) > 0:
            delta = with_scanpy.iloc[0]['mean_percentile'] - without.iloc[0]['mean_percentile']
            log(f"  {agent}: {delta:+.1f} percentile points")

if __name__ == "__main__":
    main() 