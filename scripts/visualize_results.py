#!/usr/bin/env python3
"""
BioML-bench Results Visualization Script

Takes raw replicate-level results and creates:
1. Per-task aggregation with side-by-side agent comparison plots
2. Per-task-group aggregation with domain-level analysis
3. Publication-ready figures and statistical summaries

Usage:
  python scripts/visualize_results.py --input analysis_results/all_replicates.csv --output-dir figures/
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List
import warnings

# Set publication-ready style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("Set2")
warnings.filterwarnings('ignore')

def log(message: str):
    """Simple logging with timestamp."""
    import time
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def categorize_tasks(df: pd.DataFrame) -> pd.DataFrame:
    """Add task_group column based on task_id patterns."""
    
    def get_task_group(task_id: str) -> str:
        if task_id.startswith('proteingym-dms/'):
            return 'Protein Engineering'
        elif task_id.startswith('polarishub/'):
            return 'Drug Discovery'
        elif task_id.startswith('kaggle/'):
            return 'Biomedical Imaging'
        elif task_id.startswith('manual/'):
            return 'Single Cell Omics'
        else:
            raise ValueError(f"Unknown task category for {task_id}")
    
    df['task_group'] = df['task_id'].apply(get_task_group)
    return df

def aggregate_per_task(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate results per task (across replicates)."""
    
    # Group by agent and task, calculate statistics
    agg_funcs = {
        'score': ['mean', 'std', 'median'],
        'leaderboard_percentile': ['mean', 'std', 'median'],
        'above_median': 'mean',  # Fraction above median
        'any_medal': 'mean'      # Fraction with medals
    }
    
    task_agg = df.groupby(['agent', 'task_id', 'task_group']).agg(agg_funcs).reset_index()
    
    # Flatten column names
    task_agg.columns = ['agent', 'task_id', 'task_group'] + [
        f"{col[0]}_{col[1]}" if col[1] != '' else col[0] 
        for col in task_agg.columns[3:]
    ]
    
    return task_agg

def aggregate_per_task_group(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate results per task group (across tasks and replicates)."""
    
    # Group by agent and task_group, calculate statistics across ALL replicates in that group
    agg_funcs = {
        'score': ['mean', 'std', 'median'],
        'leaderboard_percentile': ['mean', 'std', 'median'], 
        'above_median': 'mean',
        'any_medal': 'mean'
    }
    
    group_agg = df.groupby(['agent', 'task_group']).agg(agg_funcs).reset_index()
    
    # Flatten column names
    group_agg.columns = ['agent', 'task_group'] + [
        f"{col[0]}_{col[1]}" if col[1] != '' else col[0]
        for col in group_agg.columns[2:]
    ]
    
    # Add replicate count
    replicate_counts = df.groupby(['agent', 'task_group']).size().reset_index(name='total_replicates')
    group_agg = group_agg.merge(replicate_counts, on=['agent', 'task_group'])
    
    # Add unique task count
    task_counts = df.groupby(['agent', 'task_group'])['task_id'].nunique().reset_index(name='unique_tasks')
    group_agg = group_agg.merge(task_counts, on=['agent', 'task_group'])
    
    return group_agg

def plot_per_task_comparison(df: pd.DataFrame, output_dir: Path):
    """Create side-by-side plots for each task showing agent replicates."""
    
    tasks = df['task_id'].unique()
    agents = df['agent'].unique()
    
    for task in tasks:
        task_data = df[df['task_id'] == task]
        
        if len(task_data) == 0:
            continue
            
        # Create figure with subplots for different metrics
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f'Task: {task}', fontsize=16, fontweight='bold')
        
        # Plot 1: Scores
        ax1 = axes[0, 0]
        sns.boxplot(data=task_data, x='agent', y='score', ax=ax1)
        sns.swarmplot(data=task_data, x='agent', y='score', ax=ax1, color='black', size=6)
        ax1.set_title('Scores Across Replicates')
        ax1.set_ylabel('Score')
        
        # Plot 2: Leaderboard Percentiles
        ax2 = axes[0, 1]
        sns.boxplot(data=task_data, x='agent', y='leaderboard_percentile', ax=ax2)
        sns.swarmplot(data=task_data, x='agent', y='leaderboard_percentile', ax=ax2, color='black', size=6)
        ax2.set_title('Leaderboard Percentiles')
        ax2.set_ylabel('Percentile')
        ax2.set_ylim(0, 100)
        
        # Plot 3: Medal rates
        ax3 = axes[1, 0]
        medal_rates = task_data.groupby('agent')['any_medal'].mean()
        medal_rates.plot(kind='bar', ax=ax3, color='gold')
        ax3.set_title('Medal Rate')
        ax3.set_ylabel('Fraction with Medals')
        ax3.set_ylim(0, 1)
        ax3.tick_params(axis='x', rotation=45)
        
        # Plot 4: Above median rates
        ax4 = axes[1, 1] 
        above_median_rates = task_data.groupby('agent')['above_median'].mean()
        above_median_rates.plot(kind='bar', ax=ax4, color='lightblue')
        ax4.set_title('Above Median Rate')
        ax4.set_ylabel('Fraction Above Median')
        ax4.set_ylim(0, 1)
        ax4.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        # Save with safe filename
        safe_task_name = task.replace('/', '_').replace('-', '_')
        plt.savefig(output_dir / f'task_{safe_task_name}.png', dpi=300, bbox_inches='tight')
        plt.close()

def plot_per_task_group_comparison(df: pd.DataFrame, task_group_agg: pd.DataFrame, output_dir: Path):
    """Create plots comparing agents across task groups."""
    
    task_groups = df['task_group'].unique()
    
    # Overall comparison plot
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Agent Performance Across Task Groups', fontsize=16, fontweight='bold')
    
    # Plot 1: Average percentiles by task group
    ax1 = axes[0, 0]
    pivot_data = task_group_agg.pivot(index='task_group', columns='agent', values='leaderboard_percentile_mean')
    pivot_data.plot(kind='bar', ax=ax1, width=0.8)
    ax1.set_title('Average Leaderboard Percentile by Task Group')
    ax1.set_ylabel('Average Percentile')
    ax1.legend(title='Agent')
    ax1.tick_params(axis='x', rotation=45)
    
    # Plot 2: Medal rates by task group
    ax2 = axes[0, 1]
    medal_pivot = task_group_agg.pivot(index='task_group', columns='agent', values='any_medal_mean')
    medal_pivot.plot(kind='bar', ax=ax2, width=0.8)
    ax2.set_title('Medal Rate by Task Group')
    ax2.set_ylabel('Fraction with Medals')
    ax2.legend(title='Agent')
    ax2.tick_params(axis='x', rotation=45)
    
    # Plot 3: Score distributions by task group
    ax3 = axes[1, 0]
    sns.boxplot(data=df, x='task_group', y='score', hue='agent', ax=ax3)
    ax3.set_title('Score Distributions by Task Group')
    ax3.set_ylabel('Score')
    ax3.tick_params(axis='x', rotation=45)
    
    # Plot 4: Above median rates
    ax4 = axes[1, 1]
    above_median_pivot = task_group_agg.pivot(index='task_group', columns='agent', values='above_median_mean')
    above_median_pivot.plot(kind='bar', ax=ax4, width=0.8)
    ax4.set_title('Above Median Rate by Task Group')
    ax4.set_ylabel('Fraction Above Median')
    ax4.legend(title='Agent')
    ax4.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'task_group_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Individual task group detailed plots
    for task_group in task_groups:
        group_data = df[df['task_group'] == task_group]
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f'{task_group} - Detailed Analysis', fontsize=16, fontweight='bold')
        
        # Plot 1: Per-task percentiles in this group
        ax1 = axes[0, 0]
        sns.boxplot(data=group_data, x='task_id', y='leaderboard_percentile', hue='agent', ax=ax1)
        ax1.set_title('Percentiles by Task')
        ax1.set_ylabel('Leaderboard Percentile')
        ax1.tick_params(axis='x', rotation=45)
        
        # Plot 2: Per-task scores in this group  
        ax2 = axes[0, 1]
        sns.boxplot(data=group_data, x='task_id', y='score', hue='agent', ax=ax2)
        ax2.set_title('Scores by Task')
        ax2.set_ylabel('Score')
        ax2.tick_params(axis='x', rotation=45)
        
        # Plot 3: Medal distribution
        ax3 = axes[1, 0]
        medal_data = group_data.groupby(['agent', 'task_id'])['any_medal'].mean().reset_index()
        sns.barplot(data=medal_data, x='task_id', y='any_medal', hue='agent', ax=ax3)
        ax3.set_title('Medal Rate by Task')
        ax3.set_ylabel('Fraction with Medals')
        ax3.tick_params(axis='x', rotation=45)
        
        # Plot 4: Above median distribution
        ax4 = axes[1, 1]
        median_data = group_data.groupby(['agent', 'task_id'])['above_median'].mean().reset_index()
        sns.barplot(data=median_data, x='task_id', y='above_median', hue='agent', ax=ax4)
        ax4.set_title('Above Median Rate by Task')
        ax4.set_ylabel('Fraction Above Median')
        ax4.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        # Save with safe filename
        safe_group_name = task_group.replace(' ', '_').replace('-', '_')
        plt.savefig(output_dir / f'group_{safe_group_name.lower()}_detailed.png', dpi=300, bbox_inches='tight')
        plt.close()

def create_summary_tables(df: pd.DataFrame, task_agg: pd.DataFrame, group_agg: pd.DataFrame, output_dir: Path):
    """Create summary tables with statistics."""
    
    log("Creating summary tables...")
    
    # Per-task summary table
    task_summary = task_agg.copy()
    task_summary = task_summary.round(3)
    task_summary.to_csv(output_dir / 'per_task_summary.csv', index=False)
    
    # Per-task-group summary table  
    group_summary = group_agg.copy()
    group_summary = group_summary.round(3)
    group_summary.to_csv(output_dir / 'per_task_group_summary.csv', index=False)
    
    # Replicate count verification table
    replicate_counts = df.groupby(['agent', 'task_id', 'task_group']).size().reset_index(name='replicate_count')
    replicate_counts.to_csv(output_dir / 'replicate_counts.csv', index=False)
    
    # Statistical significance tests (if multiple agents)
    agents = df['agent'].unique()
    if len(agents) > 1:
        from scipy import stats
        
        significance_results = []
        
        for task_group in df['task_group'].unique():
            group_data = df[df['task_group'] == task_group]
            
            # Pairwise t-tests between agents for this task group
            for i, agent1 in enumerate(agents):
                for j, agent2 in enumerate(agents):
                    if i >= j:  # Only test each pair once
                        continue
                    
                    agent1_scores = group_data[group_data['agent'] == agent1]['leaderboard_percentile']
                    agent2_scores = group_data[group_data['agent'] == agent2]['leaderboard_percentile']
                    
                    if len(agent1_scores) > 0 and len(agent2_scores) > 0:
                        t_stat, p_value = stats.ttest_ind(agent1_scores, agent2_scores)
                        
                        significance_results.append({
                            'task_group': task_group,
                            'agent1': agent1,
                            'agent2': agent2,
                            'agent1_mean': agent1_scores.mean(),
                            'agent2_mean': agent2_scores.mean(),
                            't_statistic': t_stat,
                            'p_value': p_value,
                            'significant': p_value < 0.05
                        })
        
        if significance_results:
            sig_df = pd.DataFrame(significance_results)
            sig_df = sig_df.round(4)
            sig_df.to_csv(output_dir / 'statistical_significance.csv', index=False)

def create_publication_plots(df: pd.DataFrame, task_agg: pd.DataFrame, group_agg: pd.DataFrame, output_dir: Path):
    """Create publication-ready plots."""
    
    log("Creating publication-ready plots...")
    
    # Main figure: Task group performance overview
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('BioML-bench: Agent Performance Across Biomedical Domains', 
                 fontsize=18, fontweight='bold')
    
    # Plot 1: Average percentile by task group (with error bars)
    ax1 = axes[0, 0]
    
    # Calculate means and standard errors
    group_means = group_agg.pivot(index='task_group', columns='agent', values='leaderboard_percentile_mean')
    group_stds = group_agg.pivot(index='task_group', columns='agent', values='leaderboard_percentile_std')
    
    # Create bar plot with error bars
    x_pos = np.arange(len(group_means.index))
    width = 0.25
    
    for i, agent in enumerate(group_means.columns):
        if agent in group_means.columns:
            means = group_means[agent].fillna(0)
            stds = group_stds[agent].fillna(0)
            ax1.bar(x_pos + i*width, means, width, yerr=stds, 
                   label=agent, capsize=5, alpha=0.8)
    
    ax1.set_xlabel('Task Group')
    ax1.set_ylabel('Average Leaderboard Percentile')
    ax1.set_title('Performance by Biomedical Domain')
    ax1.set_xticks(x_pos + width)
    ax1.set_xticklabels(group_means.index, rotation=45)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Medal rates by task group
    ax2 = axes[0, 1]
    medal_pivot = group_agg.pivot(index='task_group', columns='agent', values='any_medal_mean')
    medal_pivot.plot(kind='bar', ax=ax2, width=0.8)
    ax2.set_title('Medal Achievement Rate by Domain')
    ax2.set_ylabel('Fraction of Runs with Medals')
    ax2.legend(title='Agent')
    ax2.tick_params(axis='x', rotation=45)
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Score distributions 
    ax3 = axes[1, 0]
    sns.violinplot(data=df, x='task_group', y='score', hue='agent', ax=ax3)
    ax3.set_title('Score Distributions by Domain')
    ax3.set_ylabel('Task Score')
    ax3.tick_params(axis='x', rotation=45)
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Above median rates
    ax4 = axes[1, 1]
    above_median_pivot = group_agg.pivot(index='task_group', columns='agent', values='above_median_mean')
    above_median_pivot.plot(kind='bar', ax=ax4, width=0.8)
    ax4.set_title('Above Median Achievement Rate')
    ax4.set_ylabel('Fraction Above Median')
    ax4.legend(title='Agent')
    ax4.tick_params(axis='x', rotation=45)
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'publication_main_figure.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Individual task group plots
    plot_per_task_group_comparison(df, group_agg, output_dir)

def plot_per_task_group_comparison(df: pd.DataFrame, group_agg: pd.DataFrame, output_dir: Path):
    """Create detailed plots for each task group."""
    
    for task_group in df['task_group'].unique():
        group_data = df[df['task_group'] == task_group]
        
        # Create comprehensive figure for this task group
        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        fig.suptitle(f'{task_group} - Comprehensive Analysis', fontsize=16, fontweight='bold')
        
        # Plot 1: Individual task performance
        ax1 = axes[0, 0]
        task_means = group_data.groupby(['agent', 'task_id'])['leaderboard_percentile'].mean().reset_index()
        sns.barplot(data=task_means, x='task_id', y='leaderboard_percentile', hue='agent', ax=ax1)
        ax1.set_title('Average Percentile by Task')
        ax1.set_ylabel('Leaderboard Percentile')
        ax1.tick_params(axis='x', rotation=45)
        
        # Plot 2: Score consistency (std dev)
        ax2 = axes[0, 1]
        score_stds = group_data.groupby(['agent', 'task_id'])['score'].std().reset_index()
        sns.barplot(data=score_stds, x='task_id', y='score', hue='agent', ax=ax2)
        ax2.set_title('Score Variability (Std Dev)')
        ax2.set_ylabel('Standard Deviation')
        ax2.tick_params(axis='x', rotation=45)
        
        # Plot 3: Percentile consistency
        ax3 = axes[0, 2]
        perc_stds = group_data.groupby(['agent', 'task_id'])['leaderboard_percentile'].std().reset_index()
        sns.barplot(data=perc_stds, x='task_id', y='leaderboard_percentile', hue='agent', ax=ax3)
        ax3.set_title('Percentile Variability (Std Dev)')
        ax3.set_ylabel('Standard Deviation')
        ax3.tick_params(axis='x', rotation=45)
        
        # Plot 4: All replicates scatter
        ax4 = axes[1, 0]
        for agent in group_data['agent'].unique():
            agent_data = group_data[group_data['agent'] == agent]
            ax4.scatter(agent_data['score'], agent_data['leaderboard_percentile'], 
                       label=agent, alpha=0.7, s=50)
        ax4.set_xlabel('Score')
        ax4.set_ylabel('Leaderboard Percentile')
        ax4.set_title('Score vs Percentile (All Replicates)')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        # Plot 5: Replicate reliability
        ax5 = axes[1, 1]
        reliability = group_data.groupby(['agent', 'task_id']).agg({
            'leaderboard_percentile': 'std',
            'score': 'std'
        }).reset_index()
        sns.scatterplot(data=reliability, x='score', y='leaderboard_percentile', 
                       hue=reliability['agent'], ax=ax5, s=100)
        ax5.set_xlabel('Score Std Dev')
        ax5.set_ylabel('Percentile Std Dev')
        ax5.set_title('Replicate Consistency')
        ax5.grid(True, alpha=0.3)
        
        # Plot 6: Summary statistics table as plot
        ax6 = axes[1, 2]
        ax6.axis('off')
        
        # Create summary stats for this group
        group_stats = group_agg[group_agg['task_group'] == task_group].copy()
        
        # Format for display
        summary_text = f"Summary Statistics - {task_group}\n\n"
        for _, row in group_stats.iterrows():
            agent = row['agent']
            summary_text += f"{agent}:\n"
            summary_text += f"  Avg Percentile: {row['leaderboard_percentile_mean']:.1f} ± {row['leaderboard_percentile_std']:.1f}\n"
            summary_text += f"  Tasks: {row['unique_tasks']}, Replicates: {row['total_replicates']}\n"
            summary_text += f"  Medal Rate: {row['any_medal_mean']:.1%}\n"
            summary_text += f"  Above Median: {row['above_median_mean']:.1%}\n\n"
        
        ax6.text(0.1, 0.9, summary_text, transform=ax6.transAxes, fontsize=11,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        plt.tight_layout()
        
        # Save with safe filename
        safe_group_name = task_group.replace(' ', '_').replace('-', '_')
        plt.savefig(output_dir / f'group_{safe_group_name.lower()}_detailed.png', dpi=300, bbox_inches='tight')
        plt.close()

def main():
    parser = argparse.ArgumentParser(
        description="Aggregate and visualize BioML-bench replicate results"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to all_replicates.csv file"
    )
    parser.add_argument(
        "--output-dir",
        default="publication_figures",
        help="Directory to save plots and tables (default: publication_figures)"
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for saved figures (default: 300)"
    )
    
    args = parser.parse_args()
    
    # Load data
    input_file = Path(args.input)
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    log(f"Loading data from {input_file}")
    df = pd.read_csv(input_file)
    
    if len(df) == 0:
        raise ValueError("Empty input file!")
    
    log(f"Loaded {len(df)} individual results")
    
    # Validate required columns
    required_cols = ['agent', 'task_id', 'score', 'leaderboard_percentile', 'above_median', 'any_medal']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")
    
    # Add task groupings
    df = categorize_tasks(df)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    log(f"Unique agents: {sorted(df['agent'].unique())}")
    log(f"Unique task groups: {sorted(df['task_group'].unique())}")
    log(f"Unique tasks: {len(df['task_id'].unique())}")
    
    # Replicate count validation
    replicate_counts = df.groupby(['agent', 'task_id']).size()
    log("Replicate counts per agent-task:")
    for (agent, task), count in replicate_counts.items():
        log(f"  {agent} x {task}: {count} replicates")
    
    # Aggregate at different levels
    log("Aggregating per-task statistics...")
    task_agg = aggregate_per_task(df)
    
    log("Aggregating per-task-group statistics...")
    group_agg = aggregate_per_task_group(df)
    
    # Create summary tables
    create_summary_tables(df, task_agg, group_agg, output_dir)
    
    # Create plots
    log("Creating per-task comparison plots...")
    plot_per_task_comparison(df, output_dir)
    
    log("Creating publication-ready figures...")
    create_publication_plots(df, task_agg, group_agg, output_dir)
    
    log("Analysis complete!")
    log(f"📊 Tables saved to: {output_dir}")
    log(f"📈 Figures saved to: {output_dir}")
    log("Key outputs:")
    log(f"  - per_task_summary.csv: Statistics aggregated per task")
    log(f"  - per_task_group_summary.csv: Statistics aggregated per domain") 
    log(f"  - publication_main_figure.png: Main results figure")
    log(f"  - task_*.png: Individual task comparison plots")
    log(f"  - group_*_detailed.png: Detailed domain analysis plots")

if __name__ == "__main__":
    main() 