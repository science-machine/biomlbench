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

# Define consistent agent ordering
AGENT_ORDER = ['dummy', 'mlagentbench', 'aide', 'stella', 'biomni']
AGENT_DISPLAY_NAMES = {
    'dummy': 'Dummy',
    'mlagentbench': 'MLAgentBench', 
    'aide': 'AIDE',
    'stella': 'STELLA',
    'biomni': 'Biomni'
}

def get_agent_order(agent: str) -> int:
    """Get sort order for agents."""
    try:
        return AGENT_ORDER.index(agent.lower())
    except ValueError:
        return 999

def format_agent_name(agent: str) -> str:
    """Format agent names for display."""
    return AGENT_DISPLAY_NAMES.get(agent.lower(), agent)

def log(message: str):
    """Simple logging with timestamp."""
    import time
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def sort_by_agent_order(df: pd.DataFrame) -> pd.DataFrame:
    """Sort dataframe by agent order."""
    df = df.copy()
    df['agent_order'] = df['agent'].apply(get_agent_order)
    df = df.sort_values(['agent_order', 'task_id'])
    df = df.drop('agent_order', axis=1)
    return df

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

def aggregate_per_task_group(df: pd.DataFrame, completion_rates_df: pd.DataFrame = None) -> pd.DataFrame:
    """Aggregate results per task group (across tasks and replicates)."""
    
    # Group by agent and task_group, calculate statistics across ALL replicates in that group  
    # NOTE: We do NOT aggregate scores across tasks since they are task-specific with different scales
    agg_funcs = {
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
    
    # Add completion rate analysis if available
    if completion_rates_df is not None:
        # Add task_group to completion rates
        completion_rates_with_groups = completion_rates_df.copy()
        completion_rates_with_groups['task_group'] = completion_rates_with_groups['task_id'].apply(get_task_group)
        
        # Aggregate completion rates by task group
        completion_agg = completion_rates_with_groups.groupby(['agent', 'task_group']).agg({
            'completion_rate': ['mean', 'std', 'median'],
            'total_attempts': 'sum',
            'successful_completions': 'sum'
        }).reset_index()
        
        # Flatten completion rate column names
        completion_agg.columns = ['agent', 'task_group'] + [
            f"completion_{col[0]}_{col[1]}" if col[1] != '' else f"completion_{col[0]}"
            for col in completion_agg.columns[2:]
        ]
        
        # Calculate overall completion rate per task group
        # Column names are flattened as completion_total_attempts_sum, etc.
        total_attempts_col = [col for col in completion_agg.columns if 'total_attempts' in col and 'sum' in col][0]
        successful_completions_col = [col for col in completion_agg.columns if 'successful_completions' in col and 'sum' in col][0]
        
        completion_agg['completion_rate_overall'] = (
            completion_agg[successful_completions_col] / 
            completion_agg[total_attempts_col]
        ).fillna(0.0)
        
        # Merge with main aggregation
        group_agg = group_agg.merge(completion_agg, on=['agent', 'task_group'], how='left')
    
    return group_agg

def get_task_group(task_id: str) -> str:
    """Helper function to categorize task by ID pattern."""
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

def plot_per_task_comparison(df: pd.DataFrame, output_dir: Path):
    """Create side-by-side plots for each task showing agent replicates."""
    
    tasks = df['task_id'].unique()
    agents = sorted(df['agent'].unique(), key=get_agent_order)
    
    for task in tasks:
        task_data = df[df['task_id'] == task]
        
        if len(task_data) == 0:
            continue
            
        # Create figure with subplots for different metrics
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle(f'Task: {task}', fontsize=16, fontweight='bold')
        
        # Plot 1: Scores
        ax1 = axes[0, 0]
        sns.boxplot(data=task_data, x='agent', y='score', ax=ax1, order=agents)
        sns.swarmplot(data=task_data, x='agent', y='score', ax=ax1, color='black', size=6, order=agents)
        ax1.set_title('Scores Across Replicates')
        ax1.set_ylabel('Score')
        ax1.set_xticklabels([format_agent_name(a) for a in agents])
        
        # Plot 2: Leaderboard Percentiles
        ax2 = axes[0, 1]
        sns.boxplot(data=task_data, x='agent', y='leaderboard_percentile', ax=ax2, order=agents)
        sns.swarmplot(data=task_data, x='agent', y='leaderboard_percentile', ax=ax2, color='black', size=6, order=agents)
        ax2.set_title('Leaderboard Percentiles')
        ax2.set_ylabel('Percentile')
        ax2.set_ylim(0, 100)
        ax2.set_xticklabels([format_agent_name(a) for a in agents])
        
        # Plot 3: Medal rates
        ax3 = axes[1, 0]
        medal_rates = task_data.groupby('agent')['any_medal'].mean()
        medal_rates = medal_rates.reindex(agents)
        medal_rates.plot(kind='bar', ax=ax3, color='gold')
        ax3.set_title('Medal Rate')
        ax3.set_ylabel('Fraction with Medals')
        ax3.set_ylim(0, 1)
        ax3.tick_params(axis='x', rotation=45)
        ax3.set_xticklabels([format_agent_name(a) for a in agents])
        
        # Plot 4: Above median rates
        ax4 = axes[1, 1] 
        above_median_rates = task_data.groupby('agent')['above_median'].mean()
        above_median_rates = above_median_rates.reindex(agents)
        above_median_rates.plot(kind='bar', ax=ax4, color='lightblue')
        ax4.set_title('Above Median Rate')
        ax4.set_ylabel('Fraction Above Median')
        ax4.set_ylim(0, 1)
        ax4.tick_params(axis='x', rotation=45)
        ax4.set_xticklabels([format_agent_name(a) for a in agents])
        
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

def create_summary_tables(df: pd.DataFrame, task_agg: pd.DataFrame, group_agg: pd.DataFrame, completion_rates_df: pd.DataFrame, output_dir: Path):
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
    
    # Completion rates table (if available)
    if completion_rates_df is not None:
        completion_summary = completion_rates_df.copy()
        completion_summary = completion_summary.round(3)
        completion_summary.to_csv(output_dir / 'completion_rates_summary.csv', index=False)
        
        # Overall completion rate summary by agent
        agent_completion_summary = completion_rates_df.groupby('agent').agg({
            'total_attempts': 'sum',
            'successful_completions': 'sum',
            'completion_rate': 'mean'
        }).reset_index()
        agent_completion_summary['overall_completion_rate'] = (
            agent_completion_summary['successful_completions'] / 
            agent_completion_summary['total_attempts']
        )
        agent_completion_summary = agent_completion_summary.round(3)
        agent_completion_summary.to_csv(output_dir / 'agent_completion_summary.csv', index=False)
    
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
    
    # Helper function for SEM calculation
    def calculate_sem(x):
        return np.std(x) / np.sqrt(len(x)) if len(x) > 0 else 0
    
    # Get sorted agent list
    agents = sorted(df['agent'].unique(), key=get_agent_order)
    agent_names = [format_agent_name(a) for a in agents]
    
    # Read pre-calculated statistics if available
    latex_stats_file = output_dir.parent / 'latex_tables' / 'subdomain_statistics.csv'
    if latex_stats_file.exists():
        log(f"Using pre-calculated statistics from {latex_stats_file}")
        subdomain_stats = pd.read_csv(latex_stats_file)
        # Sort by agent order
        subdomain_stats['agent_order'] = subdomain_stats['agent'].apply(get_agent_order)
        subdomain_stats = subdomain_stats.sort_values(['task_group', 'agent_order'])
    else:
        log("Warning: Pre-calculated statistics not found, will use group_agg")
        subdomain_stats = None
    
    # Main figure: Task group performance overview
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('BioML-bench: Agent Performance Across Biomedical Domains', 
                 fontsize=18, fontweight='bold')
    
    # Setup for bar plots
    x_labels = sorted(df['task_group'].unique())
    x_pos = np.arange(len(x_labels))
    width = 0.15
    
    # Plot 1: Leaderboard percentile - bar plot with error bars + jittered points
    ax1 = axes[0, 0]
    
    # Calculate task-level means for scatter points
    task_means = df.groupby(['agent', 'task_id', 'task_group'])['leaderboard_percentile'].mean().reset_index()
    
    # Prepare data for bar plot
    if subdomain_stats is not None:
        # Use pre-calculated statistics
        bar_data = subdomain_stats[['agent', 'task_group', 'leaderboard_mean', 'leaderboard_sem']].copy()
        bar_data.columns = ['agent', 'task_group', 'mean', 'sem']
    else:
        # Fallback: calculate from task means
        bar_data = task_means.groupby(['agent', 'task_group']).agg({
            'leaderboard_percentile': ['mean', calculate_sem]
        }).reset_index()
        bar_data.columns = ['agent', 'task_group', 'mean', 'sem']
    
    # Create grouped bar plot using seaborn
    sns.barplot(data=bar_data, x='task_group', y='mean', hue='agent', 
                ax=ax1, alpha=0.8, hue_order=agents, order=x_labels,
                capsize=0.1, errwidth=2, ci=None)
    
    # Add error bars manually
    bars = ax1.patches
    n_groups = len(x_labels)
    n_agents = len(agents)
    for i, (_, row) in enumerate(bar_data.iterrows()):
        bar_idx = (i % n_agents) * n_groups + (i // n_agents)
        if bar_idx < len(bars):
            bar = bars[bar_idx]
            ax1.errorbar(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                        yerr=row['sem'], fmt='none', c='black', capsize=5)
    
    # Overlay strip plot - this will now align correctly
    sns.stripplot(data=task_means, x='task_group', y='leaderboard_percentile', hue='agent',
                  ax=ax1, dodge=True, size=5, alpha=0.7, edgecolor='black', linewidth=0.5, 
                  hue_order=agents, order=x_labels, legend=False)
    
    ax1.set_xlabel('Task Domain', fontsize=16)
    ax1.set_ylabel('Leaderboard Percentile', fontsize=20)
    ax1.set_title('Performance Distribution by Domain', fontsize=24)
    ax1.tick_params(axis='x', rotation=45, labelsize=12)
    ax1.tick_params(axis='y', labelsize=12)
    
    # Fix legend to show proper names
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles[:len(agents)], agent_names, title='Agent', loc='upper left', fontsize=12, title_fontsize=13)
    
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.set_ylim(0, 100)
    
    # Save individual panel
    fig_individual = plt.figure(figsize=(8, 6))
    ax_individual = fig_individual.add_subplot(111)
    
    sns.barplot(data=bar_data, x='task_group', y='mean', hue='agent', 
                ax=ax_individual, alpha=0.8, hue_order=agents, order=x_labels,
                capsize=0.1, errwidth=2, ci=None)
    
    # Add error bars to individual panel
    bars = ax_individual.patches
    for i, (_, row) in enumerate(bar_data.iterrows()):
        bar_idx = (i % n_agents) * n_groups + (i // n_agents)
        if bar_idx < len(bars):
            bar = bars[bar_idx]
            ax_individual.errorbar(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                                  yerr=row['sem'], fmt='none', c='black', capsize=5)
    
    sns.stripplot(data=task_means, x='task_group', y='leaderboard_percentile', hue='agent',
                  ax=ax_individual, dodge=True, size=5, alpha=0.7, edgecolor='black', linewidth=0.5, 
                  hue_order=agents, order=x_labels, legend=False)
    
    ax_individual.set_xlabel('Task Domain', fontsize=15)
    ax_individual.set_ylabel('Leaderboard Percentile', fontsize=15)
    ax_individual.set_title('Performance Distribution by Domain', fontsize=20)
    ax_individual.tick_params(axis='x', rotation=45, labelsize=13)
    ax_individual.tick_params(axis='y', labelsize=13)
    handles, labels = ax_individual.get_legend_handles_labels()
    ax_individual.legend(handles[:len(agents)], agent_names, title='Agent', loc='upper left', fontsize=12, title_fontsize=13)
    ax_individual.grid(True, alpha=0.3, axis='y')
    ax_individual.set_ylim(0, 100)
    plt.tight_layout()
    plt.savefig(output_dir / 'performance_by_domain.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Medal rates - bar plot with jittered points
    ax2 = axes[0, 1]
    
    # Calculate task-level medal rates for scatter points
    task_medal_rates = df.groupby(['agent', 'task_id', 'task_group'])['any_medal'].mean().reset_index()
    task_medal_rates['any_medal_pct'] = task_medal_rates['any_medal'] * 100
    
    # Prepare data for bar plot
    if subdomain_stats is not None:
        bar_data = subdomain_stats[['agent', 'task_group', 'any_medal_mean', 'any_medal_sem']].copy()
        bar_data.columns = ['agent', 'task_group', 'mean', 'sem']
    else:
        bar_data = task_medal_rates.groupby(['agent', 'task_group']).agg({
            'any_medal_pct': ['mean', calculate_sem]
        }).reset_index()
        bar_data.columns = ['agent', 'task_group', 'mean', 'sem']
    
    # Create grouped bar plot
    sns.barplot(data=bar_data, x='task_group', y='mean', hue='agent',
                ax=ax2, alpha=0.8, hue_order=agents, order=x_labels,
                capsize=0.1, errwidth=2, ci=None)
    
    # Add error bars manually
    bars = ax2.patches
    for i, (_, row) in enumerate(bar_data.iterrows()):
        bar_idx = (i % n_agents) * n_groups + (i // n_agents)
        if bar_idx < len(bars):
            bar = bars[bar_idx]
            ax2.errorbar(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                        yerr=row['sem'], fmt='none', c='black', capsize=5)
    
    # Overlay strip plot
    sns.stripplot(data=task_medal_rates, x='task_group', y='any_medal_pct', hue='agent',
                  ax=ax2, dodge=True, size=5, alpha=0.7, edgecolor='black', linewidth=0.5, 
                  hue_order=agents, order=x_labels, legend=False)
    
    ax2.set_xlabel('Task Group')
    ax2.set_ylabel('Medal Rate (%)')
    ax2.set_title('Medal Achievement by Domain')
    ax2.tick_params(axis='x', rotation=45)
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles[:len(agents)], agent_names, title='Agent', loc='upper left')
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_ylim(-5, 105)
    
    # Save individual panel
    fig_individual = plt.figure(figsize=(8, 6))
    ax_individual = fig_individual.add_subplot(111)
    
    sns.barplot(data=bar_data, x='task_group', y='mean', hue='agent',
                ax=ax_individual, alpha=0.8, hue_order=agents, order=x_labels,
                capsize=0.1, errwidth=2, ci=None)
    
    bars = ax_individual.patches
    for i, (_, row) in enumerate(bar_data.iterrows()):
        bar_idx = (i % n_agents) * n_groups + (i // n_agents)
        if bar_idx < len(bars):
            bar = bars[bar_idx]
            ax_individual.errorbar(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                                  yerr=row['sem'], fmt='none', c='black', capsize=5)
    
    sns.stripplot(data=task_medal_rates, x='task_group', y='any_medal_pct', hue='agent',
                  ax=ax_individual, dodge=True, size=5, alpha=0.7, edgecolor='black', linewidth=0.5, 
                  hue_order=agents, order=x_labels, legend=False)
    
    ax_individual.set_xlabel('Task Group')
    ax_individual.set_ylabel('Medal Rate (%)')
    ax_individual.set_title('Medal Achievement by Domain')
    ax_individual.tick_params(axis='x', rotation=45)
    handles, labels = ax_individual.get_legend_handles_labels()
    ax_individual.legend(handles[:len(agents)], agent_names, title='Agent', loc='upper left')
    ax_individual.grid(True, alpha=0.3, axis='y')
    ax_individual.set_ylim(-5, 105)
    plt.tight_layout()
    plt.savefig(output_dir / 'medal_rates_by_domain.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: Completion rates by task group (bar plot)
    ax3 = axes[1, 0]
    if subdomain_stats is not None and 'completion_rate' in subdomain_stats.columns:
        # Prepare completion rate data
        comp_data = subdomain_stats[['agent', 'task_group', 'completion_rate']].copy()
        comp_data['completion_rate'] = comp_data['completion_rate'] / 100.0  # Convert to fraction
        
        sns.barplot(data=comp_data, x='task_group', y='completion_rate', hue='agent',
                    ax=ax3, alpha=0.8, hue_order=agents, order=x_labels)
        
        ax3.set_title('Task Completion Rate by Domain')
        ax3.set_ylabel('Completion Rate')
        ax3.set_ylim(0, 1.05)
        
        # Save individual panel
        fig_individual = plt.figure(figsize=(8, 6))
        ax_individual = fig_individual.add_subplot(111)
        sns.barplot(data=comp_data, x='task_group', y='completion_rate', hue='agent',
                    ax=ax_individual, alpha=0.8, hue_order=agents, order=x_labels)
        ax_individual.set_title('Task Completion Rate by Domain')
        ax_individual.set_ylabel('Completion Rate')
        ax_individual.set_ylim(0, 1.05)
        ax_individual.set_xlabel('Task Group')
        ax_individual.tick_params(axis='x', rotation=45)
        handles, labels = ax_individual.get_legend_handles_labels()
        ax_individual.legend(handles[:len(agents)], agent_names, title='Agent', loc='upper left')
        ax_individual.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(output_dir / 'completion_rates_by_domain.png', dpi=300, bbox_inches='tight')
        plt.close()
    else:
        # Fallback to group_agg completion rates
        if 'completion_rate_overall' in group_agg.columns:
            comp_data = group_agg[['agent', 'task_group', 'completion_rate_overall']].copy()
            sns.barplot(data=comp_data, x='task_group', y='completion_rate_overall', hue='agent',
                        ax=ax3, alpha=0.8, hue_order=agents, order=x_labels)
        
        ax3.set_title('Task Completion Rate by Domain')
        ax3.set_ylabel('Completion Rate')
        ax3.set_ylim(0, 1.05)
    
    ax3.set_xlabel('Task Group')
    ax3.tick_params(axis='x', rotation=45)
    handles, labels = ax3.get_legend_handles_labels()
    ax3.legend(handles[:len(agents)], agent_names, title='Agent', loc='upper left')
    ax3.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Above median rates - bar plot with error bars + jittered points
    ax4 = axes[1, 1]
    
    # Calculate task-level above median rates for scatter points
    task_above_median = df.groupby(['agent', 'task_id', 'task_group'])['above_median'].mean().reset_index()
    task_above_median['above_median_pct'] = task_above_median['above_median'] * 100
    
    # Prepare data for bar plot
    if subdomain_stats is not None:
        bar_data = subdomain_stats[['agent', 'task_group', 'above_median_mean', 'above_median_sem']].copy()
        bar_data.columns = ['agent', 'task_group', 'mean', 'sem']
    else:
        bar_data = task_above_median.groupby(['agent', 'task_group']).agg({
            'above_median_pct': ['mean', calculate_sem]
        }).reset_index()
        bar_data.columns = ['agent', 'task_group', 'mean', 'sem']
    
    # Create grouped bar plot
    sns.barplot(data=bar_data, x='task_group', y='mean', hue='agent',
                ax=ax4, alpha=0.8, hue_order=agents, order=x_labels,
                capsize=0.1, errwidth=2, ci=None)
    
    # Add error bars manually
    bars = ax4.patches
    for i, (_, row) in enumerate(bar_data.iterrows()):
        bar_idx = (i % n_agents) * n_groups + (i // n_agents)
        if bar_idx < len(bars):
            bar = bars[bar_idx]
            ax4.errorbar(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                        yerr=row['sem'], fmt='none', c='black', capsize=5)
    
    # Overlay strip plot
    sns.stripplot(data=task_above_median, x='task_group', y='above_median_pct', hue='agent',
                  ax=ax4, dodge=True, size=5, alpha=0.7, edgecolor='black', linewidth=0.5, 
                  hue_order=agents, order=x_labels, legend=False)
    
    ax4.set_xlabel('Task Group')
    ax4.set_ylabel('Above Median Rate (%)')
    ax4.set_title('Above Median Achievement by Domain')
    ax4.tick_params(axis='x', rotation=45)
    handles, labels = ax4.get_legend_handles_labels()
    ax4.legend(handles[:len(agents)], agent_names, title='Agent', loc='upper left')
    ax4.grid(True, alpha=0.3, axis='y')
    ax4.set_ylim(-5, 105)
    
    # Save individual panel
    fig_individual = plt.figure(figsize=(8, 6))
    ax_individual = fig_individual.add_subplot(111)
    
    sns.barplot(data=bar_data, x='task_group', y='mean', hue='agent',
                ax=ax_individual, alpha=0.8, hue_order=agents, order=x_labels,
                capsize=0.1, errwidth=2, ci=None)
    
    bars = ax_individual.patches
    for i, (_, row) in enumerate(bar_data.iterrows()):
        bar_idx = (i % n_agents) * n_groups + (i // n_agents)
        if bar_idx < len(bars):
            bar = bars[bar_idx]
            ax_individual.errorbar(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                                  yerr=row['sem'], fmt='none', c='black', capsize=5)
    
    sns.stripplot(data=task_above_median, x='task_group', y='above_median_pct', hue='agent',
                  ax=ax_individual, dodge=True, size=5, alpha=0.7, edgecolor='black', linewidth=0.5, 
                  hue_order=agents, order=x_labels, legend=False)
    
    ax_individual.set_xlabel('Task Group')
    ax_individual.set_ylabel('Above Median Rate (%)')
    ax_individual.set_title('Above Median Achievement by Domain')
    ax_individual.tick_params(axis='x', rotation=45)
    handles, labels = ax_individual.get_legend_handles_labels()
    ax_individual.legend(handles[:len(agents)], agent_names, title='Agent', loc='upper left')
    ax_individual.grid(True, alpha=0.3, axis='y')
    ax_individual.set_ylim(-5, 105)
    plt.tight_layout()
    plt.savefig(output_dir / 'above_median_by_domain.png', dpi=300, bbox_inches='tight')
    plt.close()
    
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

def create_capability_reliability_plot(df: pd.DataFrame, completion_rates_df: pd.DataFrame, output_dir: Path):
    """Create capability vs reliability scatter plot at per-task level."""
    
    log("Creating capability vs reliability scatter plot (per-task level)...")
    
    # Calculate per-task averages
    task_means = df.groupby(['agent', 'task_id', 'task_group'])['leaderboard_percentile'].mean().reset_index()
    
    # Merge with completion rates
    if completion_rates_df is not None:
        task_data = task_means.merge(
            completion_rates_df[['agent', 'task_id', 'completion_rate']], 
            on=['agent', 'task_id'], 
            how='left'
        )
        task_data['completion_rate'] = task_data['completion_rate'] * 100  # Convert to percentage
    else:
        # If no completion rates, assume 100% for all
        task_data = task_means.copy()
        task_data['completion_rate'] = 100.0
    
    # Create figure with extra space for legends
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Define colors for agents
    colors = plt.cm.tab10(np.linspace(0, 1, len(task_data['agent'].unique())))
    color_map = dict(zip(sorted(task_data['agent'].unique(), key=get_agent_order), colors))
    
    # Define markers for task groups
    markers = ['o', 's', '^', 'D']
    marker_map = dict(zip(sorted(task_data['task_group'].unique()), markers[:len(task_data['task_group'].unique())]))
    
    # Plot each agent-task pair
    for _, row in task_data.iterrows():
        ax.scatter(row['completion_rate'], row['leaderboard_percentile'], 
                  s=200,  # Fixed size for all points
                  c=[color_map[row['agent']]],
                  marker=marker_map[row['task_group']],
                  alpha=0.7,
                  edgecolors='black',
                  linewidth=2)
    
    # Create custom legend for agents
    agent_handles = []
    for agent in sorted(task_data['agent'].unique(), key=get_agent_order):
        agent_handles.append(plt.scatter([], [], c=color_map[agent], s=200, 
                                       label=format_agent_name(agent), alpha=0.7, 
                                       edgecolors='black', linewidth=2))
    
    # Create custom legend for domains
    domain_handles = []
    for domain in sorted(task_data['task_group'].unique()):
        domain_handles.append(plt.scatter([], [], c='gray', marker=marker_map[domain], 
                                        s=200, label=domain, alpha=0.7,
                                        edgecolors='black', linewidth=2))
    
    # Add legends outside the plot area
    legend1 = ax.legend(handles=agent_handles, title='Agent', loc='center left', 
                       bbox_to_anchor=(1.05, 0.7), fontsize=16, title_fontsize=18)
    ax.add_artist(legend1)
    
    legend2 = ax.legend(handles=domain_handles, title='Domain', loc='center left',
                       bbox_to_anchor=(1.05, 0.3), fontsize=16, title_fontsize=18)
    
    # Set labels and title
    ax.set_xlabel('Task Completion Rate (%)', fontsize=24)
    ax.set_ylabel('Task Leaderboard Percentile', fontsize=24)
    ax.set_title('Agent Capability vs Reliability', fontsize=28)
    
    # Set axis properties
    ax.set_xlim(-5, 105)
    ax.set_ylim(-5, 105)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='both', labelsize=16)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'capability_vs_reliability.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    log("Capability vs reliability plot saved to capability_vs_reliability.png")

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
        "--completion-rates",
        help="Path to completion_rates.csv file (optional, will auto-detect from input directory)"
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
    
    # Load replicate data
    input_file = Path(args.input)
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    log(f"Loading replicate data from {input_file}")
    df = pd.read_csv(input_file)
    
    if len(df) == 0:
        raise ValueError("Empty input file!")
    
    log(f"Loaded {len(df)} individual results")
    
    # Sort by agent order
    df = sort_by_agent_order(df)
    
    # Load completion rates data
    completion_rates_df = None
    if args.completion_rates:
        completion_rates_file = Path(args.completion_rates)
    else:
        # Auto-detect completion rates file in same directory
        completion_rates_file = input_file.parent / 'completion_rates.csv'
    
    if completion_rates_file.exists():
        log(f"Loading completion rates from {completion_rates_file}")
        completion_rates_df = pd.read_csv(completion_rates_file)
        log(f"Loaded completion rates for {len(completion_rates_df)} agent-task combinations")
    else:
        log("⚠️  No completion rates data found - will skip completion rate analysis")
    
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
    group_agg = aggregate_per_task_group(df, completion_rates_df)
    
    # Create summary tables
    create_summary_tables(df, task_agg, group_agg, completion_rates_df, output_dir)
    
    # Create plots
    log("Creating per-task comparison plots...")
    plot_per_task_comparison(df, output_dir)
    
    log("Creating publication-ready figures...")
    create_publication_plots(df, task_agg, group_agg, output_dir)

    # Create capability vs reliability plot
    create_capability_reliability_plot(df, completion_rates_df, output_dir)
    
    log("Analysis complete!")
    log(f"📊 Tables saved to: {output_dir}")
    log(f"📈 Figures saved to: {output_dir}")
    log("Key outputs:")
    log(f"  - per_task_summary.csv: Statistics aggregated per task")
    log(f"  - per_task_group_summary.csv: Statistics aggregated per domain") 
    log(f"  - publication_main_figure.png: Main results figure")
    log(f"  - capability_vs_reliability.png: Capability vs reliability scatter plot")
    log(f"  - task_*.png: Individual task comparison plots")
    log(f"  - group_*_detailed.png: Detailed domain analysis plots")

if __name__ == "__main__":
    main() 