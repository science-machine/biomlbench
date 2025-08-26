#!/usr/bin/env python3
"""
Compare BioML-bench results with and without pre-downloaded packages.

This script analyzes the impact of providing pre-downloaded packages to agents
by comparing results from two different runs focusing on:
- Protein Engineering (proteingym-dms/)
- Single Cell Omics (manual/)
- Drug Discovery (polarishub/)

Usage:
    python scripts/compare_package_impact.py \
        --set1-dir analysis_results \
        --set2-dir analysis_resultsv2 \
        --output-dir package_impact_analysis
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple
from scipy import stats
import warnings

warnings.filterwarnings('ignore')

# Set publication-ready style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("Set2")

# Define categories to analyze
CATEGORIES_TO_ANALYZE = {
    'Protein Engineering': 'proteingym-dms/',
    'Single Cell Omics': 'manual/',
    'Drug Discovery': 'polarishub/'
}

# Agent display names
AGENT_DISPLAY_NAMES = {
    'dummy': 'Dummy',
    'mlagentbench': 'MLAgentBench',
    'aide': 'AIDE',
    'stella': 'STELLA',
    'biomni': 'Biomni'
}

# Agent ordering
AGENT_ORDER = ['dummy', 'mlagentbench', 'aide', 'stella', 'biomni']

def log(message: str):
    """Simple logging with timestamp."""
    import time
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def get_task_category(task_id: str) -> str:
    """Get category for a task ID."""
    for category, prefix in CATEGORIES_TO_ANALYZE.items():
        if task_id.startswith(prefix):
            return category
    return None

def filter_to_target_categories(df: pd.DataFrame) -> pd.DataFrame:
    """Filter dataframe to only include target categories."""
    df['category'] = df['task_id'].apply(get_task_category)
    return df[df['category'].notna()].copy()

def calculate_sem(values):
    """Calculate standard error of the mean."""
    if len(values) == 0:
        return 0
    return np.std(values, ddof=1) / np.sqrt(len(values))

def load_and_prepare_data(set1_dir: Path, set2_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load data from both sets and filter to target categories."""
    
    log("Loading data from Set 1 (without packages)...")
    set1_replicates = pd.read_csv(set1_dir / 'all_replicates.csv')
    set1_completion = pd.read_csv(set1_dir / 'completion_rates.csv')
    
    log("Loading data from Set 2 (with packages)...")
    set2_replicates = pd.read_csv(set2_dir / 'all_replicates.csv')
    set2_completion = pd.read_csv(set2_dir / 'completion_rates.csv')
    
    # Filter to target categories
    set1_replicates = filter_to_target_categories(set1_replicates)
    set2_replicates = filter_to_target_categories(set2_replicates)
    
    set1_completion['category'] = set1_completion['task_id'].apply(get_task_category)
    set2_completion['category'] = set2_completion['task_id'].apply(get_task_category)
    set1_completion = set1_completion[set1_completion['category'].notna()]
    set2_completion = set2_completion[set2_completion['category'].notna()]
    
    log(f"Set 1: {len(set1_replicates)} replicates, {len(set1_completion)} completion records")
    log(f"Set 2: {len(set2_replicates)} replicates, {len(set2_completion)} completion records")
    
    return set1_replicates, set1_completion, set2_replicates, set2_completion

def calculate_comparative_statistics(set1_replicates, set1_completion, set2_replicates, set2_completion):
    """Calculate comparative statistics between the two sets."""
    
    results = []
    
    # Get all unique agent-category combinations
    agent_categories = set()
    for df in [set1_replicates, set2_replicates]:
        for agent in df['agent'].unique():
            for category in df['category'].unique():
                agent_categories.add((agent, category))
    
    for agent, category in sorted(agent_categories):
        # Filter data for this agent-category
        s1_data = set1_replicates[(set1_replicates['agent'] == agent) & 
                                  (set1_replicates['category'] == category)]
        s2_data = set2_replicates[(set2_replicates['agent'] == agent) & 
                                  (set2_replicates['category'] == category)]
        
        # Skip if no data in either set
        if len(s1_data) == 0 and len(s2_data) == 0:
            continue
        
        # Calculate task-level means for more stable comparisons
        s1_task_means = s1_data.groupby('task_id').agg({
            'leaderboard_percentile': 'mean',
            'above_median': 'mean',
            'any_medal': 'mean'
        }).reset_index() if len(s1_data) > 0 else pd.DataFrame()
        
        s2_task_means = s2_data.groupby('task_id').agg({
            'leaderboard_percentile': 'mean',
            'above_median': 'mean',
            'any_medal': 'mean'
        }).reset_index() if len(s2_data) > 0 else pd.DataFrame()
        
        # Calculate metrics for set 1
        if len(s1_task_means) > 0:
            s1_lb_mean = s1_task_means['leaderboard_percentile'].mean()
            s1_lb_sem = calculate_sem(s1_task_means['leaderboard_percentile'])
            s1_above_median = s1_task_means['above_median'].mean() * 100
            s1_above_median_sem = calculate_sem(s1_task_means['above_median'] * 100)
            s1_medal = s1_task_means['any_medal'].mean() * 100
            s1_medal_sem = calculate_sem(s1_task_means['any_medal'] * 100)
            s1_n_tasks = len(s1_task_means)
            s1_n_replicates = len(s1_data)
        else:
            s1_lb_mean = s1_lb_sem = s1_above_median = s1_above_median_sem = 0
            s1_medal = s1_medal_sem = s1_n_tasks = s1_n_replicates = 0
        
        # Calculate metrics for set 2
        if len(s2_task_means) > 0:
            s2_lb_mean = s2_task_means['leaderboard_percentile'].mean()
            s2_lb_sem = calculate_sem(s2_task_means['leaderboard_percentile'])
            s2_above_median = s2_task_means['above_median'].mean() * 100
            s2_above_median_sem = calculate_sem(s2_task_means['above_median'] * 100)
            s2_medal = s2_task_means['any_medal'].mean() * 100
            s2_medal_sem = calculate_sem(s2_task_means['any_medal'] * 100)
            s2_n_tasks = len(s2_task_means)
            s2_n_replicates = len(s2_data)
        else:
            s2_lb_mean = s2_lb_sem = s2_above_median = s2_above_median_sem = 0
            s2_medal = s2_medal_sem = s2_n_tasks = s2_n_replicates = 0
        
        # Get completion rates
        s1_comp = set1_completion[(set1_completion['agent'] == agent) & 
                                  (set1_completion['category'] == category)]
        s2_comp = set2_completion[(set2_completion['agent'] == agent) & 
                                  (set2_completion['category'] == category)]
        
        if len(s1_comp) > 0:
            s1_total_attempts = s1_comp['total_attempts'].sum()
            s1_successful = s1_comp['successful_completions'].sum()
            s1_completion_rate = (s1_successful / s1_total_attempts * 100) if s1_total_attempts > 0 else 0
        else:
            s1_completion_rate = 0
            s1_total_attempts = 0
            s1_successful = 0
        
        if len(s2_comp) > 0:
            s2_total_attempts = s2_comp['total_attempts'].sum()
            s2_successful = s2_comp['successful_completions'].sum()
            s2_completion_rate = (s2_successful / s2_total_attempts * 100) if s2_total_attempts > 0 else 0
        else:
            s2_completion_rate = 0
            s2_total_attempts = 0
            s2_successful = 0
        
        # Calculate changes
        delta_lb = s2_lb_mean - s1_lb_mean
        delta_above_median = s2_above_median - s1_above_median
        delta_medal = s2_medal - s1_medal
        delta_completion = s2_completion_rate - s1_completion_rate
        
        # Perform statistical tests if we have data in both sets
        p_value_lb = None
        p_value_above_median = None
        if len(s1_task_means) > 0 and len(s2_task_means) > 0:
            # Find common tasks
            common_tasks = set(s1_task_means['task_id']) & set(s2_task_means['task_id'])
            if len(common_tasks) > 0:
                s1_common = s1_task_means[s1_task_means['task_id'].isin(common_tasks)]
                s2_common = s2_task_means[s2_task_means['task_id'].isin(common_tasks)]
                
                # Paired t-test for leaderboard percentile
                if len(s1_common) > 1:
                    _, p_value_lb = stats.ttest_rel(
                        s2_common['leaderboard_percentile'].values,
                        s1_common['leaderboard_percentile'].values
                    )
                    _, p_value_above_median = stats.ttest_rel(
                        s2_common['above_median'].values,
                        s1_common['above_median'].values
                    )
        
        results.append({
            'agent': agent,
            'category': category,
            # Set 1 (without packages)
            's1_lb_mean': s1_lb_mean,
            's1_lb_sem': s1_lb_sem,
            's1_above_median': s1_above_median,
            's1_above_median_sem': s1_above_median_sem,
            's1_medal': s1_medal,
            's1_medal_sem': s1_medal_sem,
            's1_completion_rate': s1_completion_rate,
            's1_n_tasks': s1_n_tasks,
            's1_n_replicates': s1_n_replicates,
            's1_total_attempts': s1_total_attempts,
            's1_successful': s1_successful,
            # Set 2 (with packages)
            's2_lb_mean': s2_lb_mean,
            's2_lb_sem': s2_lb_sem,
            's2_above_median': s2_above_median,
            's2_above_median_sem': s2_above_median_sem,
            's2_medal': s2_medal,
            's2_medal_sem': s2_medal_sem,
            's2_completion_rate': s2_completion_rate,
            's2_n_tasks': s2_n_tasks,
            's2_n_replicates': s2_n_replicates,
            's2_total_attempts': s2_total_attempts,
            's2_successful': s2_successful,
            # Deltas
            'delta_lb': delta_lb,
            'delta_above_median': delta_above_median,
            'delta_medal': delta_medal,
            'delta_completion': delta_completion,
            # Statistical tests
            'p_value_lb': p_value_lb,
            'p_value_above_median': p_value_above_median
        })
    
    return pd.DataFrame(results)

def create_comparison_visualizations(stats_df: pd.DataFrame, output_dir: Path):
    """Create comparison visualizations."""
    
    log("Creating comparison visualizations...")
    
    # Filter out dummy agent for cleaner plots (optional)
    plot_df = stats_df[stats_df['agent'] != 'dummy'].copy()
    
    # Main comparison figure
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Impact of Pre-downloaded Packages on BioML-bench Performance', 
                 fontsize=18, fontweight='bold')
    
    categories = sorted(plot_df['category'].unique())
    agents = [a for a in AGENT_ORDER if a in plot_df['agent'].unique() and a != 'dummy']
    
    # Plot 1-3: Leaderboard percentile by category
    for i, category in enumerate(categories):
        ax = axes[0, i]
        cat_data = plot_df[plot_df['category'] == category]
        
        x = np.arange(len(agents))
        width = 0.35
        
        for j, agent in enumerate(agents):
            agent_data = cat_data[cat_data['agent'] == agent]
            if len(agent_data) > 0:
                row = agent_data.iloc[0]
                
                # Without packages
                ax.bar(j - width/2, row['s1_lb_mean'], width, 
                      yerr=row['s1_lb_sem'], label='Without packages' if j == 0 else '',
                      color='lightcoral', alpha=0.8, capsize=5)
                
                # With packages
                ax.bar(j + width/2, row['s2_lb_mean'], width,
                      yerr=row['s2_lb_sem'], label='With packages' if j == 0 else '',
                      color='lightgreen', alpha=0.8, capsize=5)
                
                # Add significance stars
                if row['p_value_lb'] is not None and row['p_value_lb'] < 0.05:
                    y_pos = max(row['s1_lb_mean'] + row['s1_lb_sem'], 
                               row['s2_lb_mean'] + row['s2_lb_sem']) + 2
                    ax.text(j, y_pos, '*', ha='center', va='bottom', fontsize=16)
        
        ax.set_xlabel('Agent')
        ax.set_ylabel('Leaderboard Percentile')
        ax.set_title(f'{category}')
        ax.set_xticks(x)
        ax.set_xticklabels([AGENT_DISPLAY_NAMES[a] for a in agents], rotation=45)
        ax.set_ylim(0, 100)
        ax.grid(True, alpha=0.3, axis='y')
        if i == 0:
            ax.legend()
    
    # Plot 4-6: Completion rates by category
    for i, category in enumerate(categories):
        ax = axes[1, i]
        cat_data = plot_df[plot_df['category'] == category]
        
        x = np.arange(len(agents))
        width = 0.35
        
        for j, agent in enumerate(agents):
            agent_data = cat_data[cat_data['agent'] == agent]
            if len(agent_data) > 0:
                row = agent_data.iloc[0]
                
                # Without packages
                ax.bar(j - width/2, row['s1_completion_rate'], width,
                      label='Without packages' if j == 0 else '',
                      color='lightcoral', alpha=0.8)
                
                # With packages
                ax.bar(j + width/2, row['s2_completion_rate'], width,
                      label='With packages' if j == 0 else '',
                      color='lightgreen', alpha=0.8)
        
        ax.set_xlabel('Agent')
        ax.set_ylabel('Completion Rate (%)')
        ax.set_title(f'{category} - Completion Rates')
        ax.set_xticks(x)
        ax.set_xticklabels([AGENT_DISPLAY_NAMES[a] for a in agents], rotation=45)
        ax.set_ylim(0, 105)
        ax.grid(True, alpha=0.3, axis='y')
        if i == 0:
            ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'package_impact_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create delta plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle('Change in Performance with Pre-downloaded Packages', 
                 fontsize=16, fontweight='bold')
    
    for i, category in enumerate(categories):
        ax = axes[i]
        cat_data = plot_df[plot_df['category'] == category]
        
        # Prepare data for grouped bar plot
        metrics = ['Percentile', 'Above Median', 'Medal Rate', 'Completion']
        x = np.arange(len(agents))
        width = 0.2
        
        for j, agent in enumerate(agents):
            agent_data = cat_data[cat_data['agent'] == agent]
            if len(agent_data) > 0:
                row = agent_data.iloc[0]
                
                deltas = [
                    row['delta_lb'],
                    row['delta_above_median'],
                    row['delta_medal'],
                    row['delta_completion']
                ]
                
                positions = x[j] + np.arange(len(metrics)) * width - (len(metrics)-1) * width / 2
                colors = ['green' if d > 0 else 'red' for d in deltas]
                ax.bar(positions, deltas, width, color=colors, alpha=0.7,
                      label=AGENT_DISPLAY_NAMES[agent])
        
        ax.set_xlabel('Metric')
        ax.set_ylabel('Change (percentage points)')
        ax.set_title(f'{category}')
        ax.set_xticks(x)
        ax.set_xticklabels([AGENT_DISPLAY_NAMES[a] for a in agents])
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        ax.grid(True, alpha=0.3, axis='y')
    
    # Add custom legend for metrics
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='gray', alpha=0.7, label=m) 
        for m in metrics
    ]
    fig.legend(handles=legend_elements, loc='center', bbox_to_anchor=(0.5, -0.05), ncol=4)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'package_impact_deltas.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create scatter plot - capability vs reliability
    fig, ax = plt.subplots(figsize=(10, 8))
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(agents)))
    color_map = dict(zip(agents, colors))
    
    for agent in agents:
        agent_data = plot_df[plot_df['agent'] == agent]
        
        # Plot points and connect with arrows
        for _, row in agent_data.iterrows():
            # Without packages (start point)
            ax.scatter(row['s1_completion_rate'], row['s1_lb_mean'], 
                      s=200, c=[color_map[agent]], marker='o', alpha=0.5,
                      edgecolors='black', linewidth=2)
            
            # With packages (end point)
            ax.scatter(row['s2_completion_rate'], row['s2_lb_mean'],
                      s=200, c=[color_map[agent]], marker='s', alpha=0.8,
                      edgecolors='black', linewidth=2,
                      label=f"{AGENT_DISPLAY_NAMES[agent]} - {row['category']}")
            
            # Arrow from without to with
            ax.annotate('', xy=(row['s2_completion_rate'], row['s2_lb_mean']),
                       xytext=(row['s1_completion_rate'], row['s1_lb_mean']),
                       arrowprops=dict(arrowstyle='->', color=color_map[agent],
                                     alpha=0.6, lw=2))
    
    ax.set_xlabel('Completion Rate (%)', fontsize=14)
    ax.set_ylabel('Mean Leaderboard Percentile', fontsize=14)
    ax.set_title('Impact of Pre-downloaded Packages: Capability vs Reliability', fontsize=16)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-5, 105)
    ax.set_ylim(-5, 105)
    
    # Add legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    
    # Add annotation
    ax.text(0.02, 0.98, 'Circles: Without packages\nSquares: With packages\nArrows: Change direction',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_dir / 'capability_reliability_change.png', dpi=300, bbox_inches='tight')
    plt.close()

def generate_latex_tables(stats_df: pd.DataFrame, output_dir: Path):
    """Generate LaTeX tables for the comparison."""
    
    log("Generating LaTeX tables...")
    
    # Sort by agent order and category
    stats_df['agent_order'] = stats_df['agent'].apply(lambda x: AGENT_ORDER.index(x) if x in AGENT_ORDER else 999)
    stats_df = stats_df.sort_values(['category', 'agent_order'])
    
    # Table 1: Detailed comparison by domain
    latex = []
    latex.append("\\begin{table*}[htbp]")
    latex.append("\\centering")
    latex.append("\\caption{Impact of Pre-downloaded Packages on BioML-bench Performance}")
    latex.append("\\label{tab:package_impact}")
    latex.append("\\resizebox{\\textwidth}{!}{%")
    latex.append("\\begin{tabular}{llcccccc}")
    latex.append("\\toprule")
    latex.append("\\multirow{2}{*}{Domain} & \\multirow{2}{*}{Agent} & \\multicolumn{2}{c}{Leaderboard Percentile} & \\multicolumn{2}{c}{Completion Rate (\\%)} & \\multicolumn{2}{c}{Change} \\\\")
    latex.append("\\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}")
    latex.append(" & & Without Packages & With Packages & Without Packages & With Packages & $\\Delta$ Percentile & $\\Delta$ Completion \\\\")
    latex.append("\\midrule")
    
    for i, category in enumerate(sorted(stats_df['category'].unique())):
        if i > 0:
            latex.append("\\midrule")
        
        cat_data = stats_df[stats_df['category'] == category]
        
        for j, (_, row) in enumerate(cat_data.iterrows()):
            # Skip dummy agent
            if row['agent'] == 'dummy':
                continue
            
            # Format values
            s1_lb = f"{row['s1_lb_mean']:.1f} ± {row['s1_lb_sem']:.1f}"
            s2_lb = f"{row['s2_lb_mean']:.1f} ± {row['s2_lb_sem']:.1f}"
            
            # Bold if improvement is significant
            if row['p_value_lb'] is not None and row['p_value_lb'] < 0.05 and row['delta_lb'] > 0:
                s2_lb = f"\\textbf{{{s2_lb}}}"
            
            s1_comp = f"{row['s1_completion_rate']:.1f}"
            s2_comp = f"{row['s2_completion_rate']:.1f}"
            
            # Bold if completion rate improved substantially (>5%)
            if row['delta_completion'] > 5:
                s2_comp = f"\\textbf{{{s2_comp}}}"
            
            delta_lb = f"{row['delta_lb']:+.1f}"
            if row['p_value_lb'] is not None and row['p_value_lb'] < 0.05:
                delta_lb += "*"
            
            delta_comp = f"{row['delta_completion']:+.1f}"
            
            # Color code deltas
            if row['delta_lb'] > 0:
                delta_lb = f"\\textcolor{{darkgreen}}{{{delta_lb}}}"
            elif row['delta_lb'] < 0:
                delta_lb = f"\\textcolor{{red}}{{{delta_lb}}}"
            
            if row['delta_completion'] > 0:
                delta_comp = f"\\textcolor{{darkgreen}}{{{delta_comp}}}"
            elif row['delta_completion'] < 0:
                delta_comp = f"\\textcolor{{red}}{{{delta_comp}}}"
            
            # Domain name only on first row
            domain_str = category if j == 0 else ""
            
            agent_display = AGENT_DISPLAY_NAMES[row['agent']]
            
            latex.append(f"{domain_str} & {agent_display} & {s1_lb} & {s2_lb} & {s1_comp} & {s2_comp} & {delta_lb} & {delta_comp} \\\\")
    
    latex.append("\\bottomrule")
    latex.append("\\end{tabular}%")
    latex.append("}")
    latex.append("\\vspace{2mm}")
    latex.append("\\caption*{Note: * indicates statistical significance (p < 0.05) using paired t-test on task-level means. Green indicates improvement, red indicates degradation.}")
    latex.append("\\end{table*}")
    
    # Save detailed table
    with open(output_dir / 'table_package_impact_detailed.tex', 'w') as f:
        f.write('\n'.join(latex))
    
    # Table 2: Summary statistics
    latex_summary = []
    latex_summary.append("\\begin{table}[htbp]")
    latex_summary.append("\\centering")
    latex_summary.append("\\caption{Summary of Package Impact Across All Domains}")
    latex_summary.append("\\label{tab:package_impact_summary}")
    latex_summary.append("\\begin{tabular}{lccc}")
    latex_summary.append("\\toprule")
    latex_summary.append("Agent & Avg $\\Delta$ Percentile & Avg $\\Delta$ Completion (\\%) & Domains Improved \\\\")
    latex_summary.append("\\midrule")
    
    # Calculate summary stats per agent
    agent_summary = stats_df[stats_df['agent'] != 'dummy'].groupby('agent').agg({
        'delta_lb': 'mean',
        'delta_completion': 'mean',
        'category': 'count'
    }).reset_index()
    
    # Count domains with improvement
    improvements = stats_df[stats_df['agent'] != 'dummy'].groupby('agent').apply(
        lambda x: (x['delta_lb'] > 0).sum()
    ).reset_index(name='domains_improved')
    
    agent_summary = agent_summary.merge(improvements, on='agent')
    
    # Sort by agent order
    agent_summary['agent_order'] = agent_summary['agent'].apply(
        lambda x: AGENT_ORDER.index(x) if x in AGENT_ORDER else 999
    )
    agent_summary = agent_summary.sort_values('agent_order')
    
    for _, row in agent_summary.iterrows():
        agent_display = AGENT_DISPLAY_NAMES[row['agent']]
        
        delta_lb_str = f"{row['delta_lb']:+.1f}"
        delta_comp_str = f"{row['delta_completion']:+.1f}"
        
        # Color code
        if row['delta_lb'] > 0:
            delta_lb_str = f"\\textcolor{{darkgreen}}{{{delta_lb_str}}}"
        elif row['delta_lb'] < 0:
            delta_lb_str = f"\\textcolor{{red}}{{{delta_lb_str}}}"
        
        if row['delta_completion'] > 0:
            delta_comp_str = f"\\textcolor{{darkgreen}}{{{delta_comp_str}}}"
        elif row['delta_completion'] < 0:
            delta_comp_str = f"\\textcolor{{red}}{{{delta_comp_str}}}"
        
        domains_str = f"{row['domains_improved']}/{row['category']}"
        
        latex_summary.append(f"{agent_display} & {delta_lb_str} & {delta_comp_str} & {domains_str} \\\\")
    
    latex_summary.append("\\bottomrule")
    latex_summary.append("\\end{tabular}")
    latex_summary.append("\\end{table}")
    
    # Save summary table
    with open(output_dir / 'table_package_impact_summary.tex', 'w') as f:
        f.write('\n'.join(latex_summary))
    
    # Create standalone document
    standalone = f"""\\documentclass{{article}}
\\usepackage{{booktabs}}
\\usepackage{{multirow}}
\\usepackage{{graphicx}}
\\usepackage{{xcolor}}
\\definecolor{{darkgreen}}{{rgb}}{{0,0.5,0}}

\\begin{{document}}

{''.join(latex)}

\\vspace{{1cm}}

{''.join(latex_summary)}

\\end{{document}}
"""
    
    with open(output_dir / 'tables_standalone.tex', 'w') as f:
        f.write(standalone)

def main():
    parser = argparse.ArgumentParser(
        description="Compare BioML-bench results with and without pre-downloaded packages"
    )
    parser.add_argument(
        "--set1-dir",
        default="analysis_results",
        help="Directory for Set 1 results (without packages)"
    )
    parser.add_argument(
        "--set2-dir", 
        default="analysis_resultsv2",
        help="Directory for Set 2 results (with packages)"
    )
    parser.add_argument(
        "--output-dir",
        default="package_impact_analysis",
        help="Directory to save comparison outputs"
    )
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Load and prepare data
    set1_replicates, set1_completion, set2_replicates, set2_completion = load_and_prepare_data(
        Path(args.set1_dir), Path(args.set2_dir)
    )
    
    # Calculate comparative statistics
    log("Calculating comparative statistics...")
    stats_df = calculate_comparative_statistics(
        set1_replicates, set1_completion, set2_replicates, set2_completion
    )
    
    # Save statistics
    stats_df.to_csv(output_dir / 'comparison_statistics.csv', index=False)
    log(f"Saved comparison statistics to {output_dir / 'comparison_statistics.csv'}")
    
    # Create visualizations
    create_comparison_visualizations(stats_df, output_dir)
    log(f"Saved visualizations to {output_dir}")
    
    # Generate LaTeX tables
    generate_latex_tables(stats_df, output_dir)
    log(f"Saved LaTeX tables to {output_dir}")
    
    # Print summary
    log("\n" + "="*60)
    log("PACKAGE IMPACT ANALYSIS COMPLETE")
    log("="*60)
    
    # Overall statistics
    overall_stats = stats_df[stats_df['agent'] != 'dummy'].groupby('agent').agg({
        'delta_lb': 'mean',
        'delta_completion': 'mean'
    }).round(1)
    
    log("\nOverall Impact by Agent:")
    for agent, row in overall_stats.iterrows():
        log(f"  {AGENT_DISPLAY_NAMES[agent]}: "
            f"Δ Percentile = {row['delta_lb']:+.1f}, "
            f"Δ Completion = {row['delta_completion']:+.1f}%")
    
    log("\nOutputs generated:")
    log(f"  - {output_dir / 'comparison_statistics.csv'}: Detailed statistics")
    log(f"  - {output_dir / 'package_impact_comparison.png'}: Main comparison figure")
    log(f"  - {output_dir / 'package_impact_deltas.png'}: Change visualization")
    log(f"  - {output_dir / 'capability_reliability_change.png'}: Scatter plot")
    log(f"  - {output_dir / 'table_package_impact_detailed.tex'}: Detailed LaTeX table")
    log(f"  - {output_dir / 'table_package_impact_summary.tex'}: Summary LaTeX table")

if __name__ == "__main__":
    main() 