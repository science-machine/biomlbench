#!/usr/bin/env python3
"""
Visualize Single Cell Omics performance comparing agents with and without pre-installed scanpy.

Creates a bar plot showing leaderboard performance for Single Cell Omics tasks,
with agents displayed as: Dummy, MLAgentBench, MLAgentBench (+scanpy), AIDE, AIDE (+scanpy), 
STELLA, STELLA (+scanpy), Biomni, Biomni (+scanpy).

Usage:
    python scripts/visualize_single_cell_package_impact.py \
        --set1-dir analysis_results \
        --set2-dir analysis_resultsv2 \
        --output-dir single_cell_figures
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

# Agent display names
AGENT_DISPLAY_NAMES = {
    'mlagentbench': 'MLAgentBench',
    'mlagentbench_scanpy': 'MLAgentBench (+scanpy)',
    'aide': 'AIDE',
    'aide_scanpy': 'AIDE (+scanpy)',
    'stella': 'STELLA',
    'stella_scanpy': 'STELLA (+scanpy)',
    'biomni': 'Biomni',
    'biomni_scanpy': 'Biomni (+scanpy)'
}

# Agent ordering
AGENT_ORDER = [
    'mlagentbench', 'mlagentbench_scanpy',
    'aide', 'aide_scanpy',
    'stella', 'stella_scanpy',
    'biomni', 'biomni_scanpy'
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

def load_and_prepare_data(set1_dir: Path, set2_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load data from both sets and filter to Single Cell Omics tasks only."""
    
    log("Loading data from Set 1 (without scanpy)...")
    set1_replicates = pd.read_csv(set1_dir / 'all_replicates.csv')
    
    log("Loading data from Set 2 (with scanpy)...")
    set2_replicates = pd.read_csv(set2_dir / 'all_replicates.csv')
    
    # Filter to Single Cell Omics tasks only
    set1_single_cell = set1_replicates[set1_replicates['task_id'].str.startswith('manual/')].copy()
    set2_single_cell = set2_replicates[set2_replicates['task_id'].str.startswith('manual/')].copy()
    
    log(f"Set 1: {len(set1_single_cell)} Single Cell Omics replicates")
    log(f"Set 2: {len(set2_single_cell)} Single Cell Omics replicates")
    
    return set1_single_cell, set2_single_cell

def calculate_statistics(set1_data: pd.DataFrame, set2_data: pd.DataFrame) -> pd.DataFrame:
    """Calculate statistics for plotting."""
    
    results = []
    
    # Process Set 1 (without scanpy)
    for agent in set1_data['agent'].unique():
        agent_data = set1_data[set1_data['agent'] == agent]
        
        # Calculate task-level means first
        task_means = agent_data.groupby('task_id')['leaderboard_percentile'].mean()
        
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
            'n_replicates': len(agent_data)
        })
    
    # Process Set 2 (with scanpy) - exclude dummy
    for agent in set2_data['agent'].unique():
        if agent.lower() == 'dummy':
            continue  # Skip dummy for +scanpy version
            
        agent_data = set2_data[set2_data['agent'] == agent]
        
        # Calculate task-level means first
        task_means = agent_data.groupby('task_id')['leaderboard_percentile'].mean()
        
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
            'n_replicates': len(agent_data)
        })
    
    return pd.DataFrame(results)


def create_single_cell_performance_plot_with_points(set1_data: pd.DataFrame, set2_data: pd.DataFrame, 
                                                   stats_df: pd.DataFrame, output_dir: Path):
    """Create the plot with individual task points overlaid."""
    
    log("Creating Single Cell Omics performance plot with task-level points...")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Sort data by agent order
    plot_data = []
    for agent_label in AGENT_ORDER:
        row = stats_df[stats_df['agent_label'] == agent_label]
        if len(row) > 0:
            plot_data.append(row.iloc[0])
    
    plot_df = pd.DataFrame(plot_data)
    
    # Prepare bar positions
    x = np.arange(len(plot_df))
    width = 0.8
    
    # Create bars and overlay points
    for idx, row in plot_df.iterrows():
        base_agent = row['agent'].lower()
        color = BASE_COLORS[base_agent]

        
        # Get task-level data for points
        if row['has_scanpy']:
            agent_task_data = set2_data[set2_data['agent'] == row['agent']].groupby('task_id')['leaderboard_percentile'].mean()
        else:
            agent_task_data = set1_data[set1_data['agent'] == row['agent']].groupby('task_id')['leaderboard_percentile'].mean()
        
        # Create bar - use position index instead of dataframe index
        bar_position = list(plot_df.index).index(idx)
        if row['has_scanpy']:
            bar = ax.bar(bar_position, row['mean_percentile'], width, 
                         yerr=row['sem_percentile'],
                         color=color, alpha=0.6,
                         edgecolor='black', linewidth=2,
                         hatch='///',
                         capsize=5, error_kw={'linewidth': 2})
        else:
            bar = ax.bar(bar_position, row['mean_percentile'], width,
                         yerr=row['sem_percentile'],
                         color=color, alpha=0.8,
                         edgecolor='black', linewidth=2,
                         capsize=5, error_kw={'linewidth': 2})

        
        # Overlay individual task points
        if len(agent_task_data) > 0:
            # Add jitter to x position for better visibility
            jitter = np.random.normal(0, 0.1, len(agent_task_data))
            x_positions = np.full(len(agent_task_data), bar_position) + jitter
            
            ax.scatter(x_positions, agent_task_data.values, 
                      color='black', s=50, alpha=0.7, 
                      edgecolors='white', linewidth=1, zorder=10)
    
    # Customize plot
    ax.set_xlabel('Agent', fontsize=16)
    ax.set_ylabel('Leaderboard Percentile', fontsize=16)
    ax.set_title('Single Cell Omics Performance: Impact of Pre-installed scanpy', fontsize=18, fontweight='bold')
    
    # Set x-axis labels
    labels = [AGENT_DISPLAY_NAMES[row['agent_label']] for _, row in plot_df.iterrows()]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=12)
    
    # Set y-axis properties
    ax.set_ylim(0, 100)
    ax.tick_params(axis='y', labelsize=12)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add legend
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    legend_elements = [
        Patch(facecolor='gray', edgecolor='black', label='Without scanpy', alpha=0.8),
        Patch(facecolor='gray', edgecolor='black', hatch='///', label='With scanpy (+scanpy)', alpha=0.6),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='black', 
               markersize=8, label='Individual tasks', markeredgecolor='white')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'single_cell_omics_scanpy_comparison_with_points.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description="Visualize Single Cell Omics performance with scanpy comparison"
    )
    parser.add_argument(
        "--set1-dir",
        default="analysis_results",
        help="Directory for Set 1 results (without scanpy)"
    )
    parser.add_argument(
        "--set2-dir",
        default="analysis_resultsv2",
        help="Directory for Set 2 results (with scanpy)"
    )
    parser.add_argument(
        "--output-dir",
        default="single_cell_figures",
        help="Directory to save outputs"
    )
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Load and prepare data
    set1_data, set2_data = load_and_prepare_data(Path(args.set1_dir), Path(args.set2_dir))
    
    # Calculate statistics
    log("Calculating statistics...")
    if "dummy" in set1_data['agent'].unique():
        set1_data = set1_data[set1_data['agent'] != "dummy"]
    if "dummy" in set2_data['agent'].unique():
        set2_data = set2_data[set2_data['agent'] != "dummy"]
    stats_df = calculate_statistics(set1_data, set2_data)
    
    # Save statistics
    stats_df.to_csv(output_dir / 'single_cell_statistics.csv', index=False)
    log(f"Saved statistics to {output_dir / 'single_cell_statistics.csv'}")
    
    # Create plots
    create_single_cell_performance_plot_with_points(set1_data, set2_data, stats_df, output_dir)
    
    # Print summary
    log("\n" + "="*60)
    log("SINGLE CELL OMICS ANALYSIS COMPLETE")
    log("="*60)
    
    log("\nPerformance Summary:")
    for _, row in stats_df.iterrows():
        agent_name = AGENT_DISPLAY_NAMES[row['agent_label']]
        log(f"  {agent_name}: {row['mean_percentile']:.1f} ± {row['sem_percentile']:.1f}")
    
    log("\nOutputs generated:")
    log(f"  - {output_dir / 'single_cell_statistics.csv'}: Statistics summary")
    log(f"  - {output_dir / 'single_cell_omics_scanpy_comparison.png'}: Main comparison plot")
    log(f"  - {output_dir / 'single_cell_omics_scanpy_comparison_clean.png'}: Clean version")
    log(f"  - {output_dir / 'single_cell_omics_scanpy_comparison_with_points.png'}: With task points")

if __name__ == "__main__":
    main() 