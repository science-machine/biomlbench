#!/usr/bin/env python3
"""
Generate LaTeX summary tables for BioML-bench results.

Creates two tables:
1. Results broken down by subdomain (task group)
2. Overall results across all tasks (MLE-bench style)

Usage:
    python scripts/generate_latex_tables.py --data-dir analysis_results/ --output-dir latex_tables/
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Tuple

def log(message: str):
    """Simple logging with timestamp."""
    import time
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}")

def get_task_group(task_id: str) -> str:
    """Categorize task by ID pattern."""
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

def calculate_sem(values):
    """Calculate standard error of the mean."""
    if len(values) == 0:
        return 0
    return np.std(values, ddof=1) / np.sqrt(len(values))

def format_value_with_sem(mean, sem, is_percentage=False):
    """Format value as mean ± SEM."""
    if is_percentage:
        return f"{mean:.1f} ± {sem:.1f}"
    else:
        return f"{mean:.2f} ± {sem:.2f}"

def bold_if_best(value, is_best):
    """Wrap value in LaTeX bold if it's the best."""
    if is_best:
        return f"\\textbf{{{value}}}"
    return value

def format_agent_name(agent: str) -> str:
    """Format agent names for display."""
    name_map = {
        'mlagentbench': 'MLAgentBench',
        'stella': 'STELLA',
        'aide': 'AIDE',
        'biomni': 'Biomni',
        'dummy': 'Dummy'
    }
    return name_map.get(agent.lower(), agent)

def get_agent_order(agent: str) -> int:
    """Get sort order for agents."""
    order_map = {
        'dummy': 0,
        'mlagentbench': 1,
        'aide': 2,
        'stella': 3,
        'biomni': 4
    }
    return order_map.get(agent.lower(), 999)

def calculate_subdomain_statistics(df: pd.DataFrame, completion_df: pd.DataFrame = None) -> pd.DataFrame:
    """Calculate statistics for each agent/subdomain combination."""
    
    # Add task groups
    df['task_group'] = df['task_id'].apply(get_task_group)
    
    results = []
    
    # Get unique agents and sort them according to our custom order
    agents = sorted(df['agent'].unique(), key=get_agent_order)
    
    for task_group in sorted(df['task_group'].unique()):
        group_data = df[df['task_group'] == task_group]
        
        for agent in agents:
            agent_data = group_data[group_data['agent'] == agent]
            
            if len(agent_data) == 0:
                continue
            
            # Calculate task-level means first
            agg_dict = {
                'leaderboard_percentile': 'mean',
                'above_median': 'mean',
                'any_medal': 'mean'
            }
            
            # Add rank if available
            if 'rank' in agent_data.columns:
                agg_dict['rank'] = 'mean'
            
            task_means = agent_data.groupby('task_id').agg(agg_dict).reset_index()
            
            # Calculate statistics across tasks
            lb_mean = task_means['leaderboard_percentile'].mean()
            lb_sem = calculate_sem(task_means['leaderboard_percentile'])
            
            above_median_mean = task_means['above_median'].mean() * 100  # Convert to percentage
            above_median_sem = calculate_sem(task_means['above_median'] * 100)
            
            any_medal_mean = task_means['any_medal'].mean() * 100  # Convert to percentage
            any_medal_sem = calculate_sem(task_means['any_medal'] * 100)
            
            # Calculate rank statistics if available
            rank_mean = None
            rank_sem = None
            if 'rank' in task_means.columns:
                rank_mean = task_means['rank'].mean()
                rank_sem = calculate_sem(task_means['rank'])
            
            # Get completion rate if available
            completion_rate = None
            if completion_df is not None:
                comp_data = completion_df[
                    (completion_df['agent'] == agent) & 
                    (completion_df['task_id'].apply(get_task_group) == task_group)
                ]
                if len(comp_data) > 0:
                    total_attempts = comp_data['total_attempts'].sum()
                    successful_completions = comp_data['successful_completions'].sum()
                    completion_rate = (successful_completions / total_attempts * 100) if total_attempts > 0 else 0
            
            results.append({
                'agent': agent,
                'task_group': task_group,
                'n_tasks': len(task_means),
                'n_replicates': len(agent_data),
                'leaderboard_mean': lb_mean,
                'leaderboard_sem': lb_sem,
                'above_median_mean': above_median_mean,
                'above_median_sem': above_median_sem,
                'any_medal_mean': any_medal_mean,
                'any_medal_sem': any_medal_sem,
                'rank_mean': rank_mean,
                'rank_sem': rank_sem,
                'completion_rate': completion_rate
            })
    
    return pd.DataFrame(results)

def calculate_overall_statistics(df: pd.DataFrame, completion_df: pd.DataFrame = None) -> pd.DataFrame:
    """Calculate overall statistics across all tasks."""
    
    results = []
    
    # Get unique agents and sort them according to our custom order
    agents = sorted(df['agent'].unique(), key=get_agent_order)
    
    for agent in agents:
        agent_data = df[df['agent'] == agent]
        
        # Calculate task-level means first
        agg_dict = {
            'leaderboard_percentile': 'mean',
            'above_median': 'mean',
            'any_medal': 'mean'
        }
        
        # Add rank if available
        if 'rank' in agent_data.columns:
            agg_dict['rank'] = 'mean'
        
        task_means = agent_data.groupby('task_id').agg(agg_dict).reset_index()
        
        # Calculate statistics across all tasks
        lb_mean = task_means['leaderboard_percentile'].mean()
        lb_sem = calculate_sem(task_means['leaderboard_percentile'])
        
        above_median_mean = task_means['above_median'].mean() * 100
        above_median_sem = calculate_sem(task_means['above_median'] * 100)
        
        any_medal_mean = task_means['any_medal'].mean() * 100
        any_medal_sem = calculate_sem(task_means['any_medal'] * 100)
        
        # Calculate rank statistics if available
        rank_mean = None
        rank_sem = None
        if 'rank' in task_means.columns:
            rank_mean = task_means['rank'].mean()
            rank_sem = calculate_sem(task_means['rank'])
        
        # Get overall completion rate
        completion_rate = None
        if completion_df is not None:
            comp_data = completion_df[completion_df['agent'] == agent]
            if len(comp_data) > 0:
                total_attempts = comp_data['total_attempts'].sum()
                successful_completions = comp_data['successful_completions'].sum()
                completion_rate = (successful_completions / total_attempts * 100) if total_attempts > 0 else 0
        
        results.append({
            'agent': agent,
            'n_tasks': len(task_means),
            'n_replicates': len(agent_data),
            'leaderboard_mean': lb_mean,
            'leaderboard_sem': lb_sem,
            'above_median_mean': above_median_mean,
            'above_median_sem': above_median_sem,
            'any_medal_mean': any_medal_mean,
            'any_medal_sem': any_medal_sem,
            'rank_mean': rank_mean,
            'rank_sem': rank_sem,
            'completion_rate': completion_rate
        })
    
    return pd.DataFrame(results)

def generate_subdomain_latex_table(stats_df: pd.DataFrame) -> str:
    """Generate LaTeX table for subdomain results."""
    
    latex = []
    latex.append("\\begin{table}[htbp]")
    latex.append("\\centering")
    latex.append("\\caption{BioML-bench Results by Domain}")
    latex.append("\\label{tab:bioml_results_by_domain}")
    latex.append("\\resizebox{\\textwidth}{!}{%")
    # Check if rank data is available
    has_rank = 'rank_mean' in stats_df.columns and stats_df['rank_mean'].notna().any()
    
    if has_rank:
        latex.append("\\begin{tabular}{llccccc}")
        latex.append("\\toprule")
        latex.append("Domain & Agent & Leaderboard Percentile & Mean Rank & Above Median (\\%) & Any Medal (\\%) & Completion Rate (\\%) \\\\")
    else:
        latex.append("\\begin{tabular}{llcccc}")
        latex.append("\\toprule")
        latex.append("Domain & Agent & Leaderboard Percentile & Above Median (\\%) & Any Medal (\\%) & Completion Rate (\\%) \\\\")
    latex.append("\\midrule")
    
    # Process each task group
    for i, task_group in enumerate(sorted(stats_df['task_group'].unique())):
        if i > 0:
            latex.append("\\midrule")
        
        group_data = stats_df[stats_df['task_group'] == task_group]
        
        # Find best values for bolding (excluding dummy for completion rate)
        best_lb = group_data['leaderboard_mean'].max()
        best_above_median = group_data['above_median_mean'].max()
        best_medal = group_data['any_medal_mean'].max()
        
        # For rank, lower is better
        best_rank = None
        if has_rank and 'rank_mean' in group_data.columns:
            rank_values = group_data['rank_mean'].dropna()
            if len(rank_values) > 0:
                best_rank = rank_values.min()
        
        # For completion rate, exclude dummy agent when finding best
        non_dummy_data = group_data[group_data['agent'].str.lower() != 'dummy']
        best_completion = non_dummy_data['completion_rate'].max() if len(non_dummy_data) > 0 and 'completion_rate' in non_dummy_data.columns else None
        
        for j, (_, row) in enumerate(group_data.iterrows()):
            # Format values
            lb_str = format_value_with_sem(row['leaderboard_mean'], row['leaderboard_sem'])
            above_median_str = format_value_with_sem(row['above_median_mean'], row['above_median_sem'], is_percentage=True)
            any_medal_str = format_value_with_sem(row['any_medal_mean'], row['any_medal_sem'], is_percentage=True)
            
            # Bold if best
            lb_str = bold_if_best(lb_str, abs(row['leaderboard_mean'] - best_lb) < 0.01)
            above_median_str = bold_if_best(above_median_str, abs(row['above_median_mean'] - best_above_median) < 0.01)
            any_medal_str = bold_if_best(any_medal_str, abs(row['any_medal_mean'] - best_medal) < 0.01)
            
            # Completion rate - special handling for dummy
            if row['agent'].lower() == 'dummy':
                comp_str = "NA"
            elif row['completion_rate'] is not None:
                comp_str = f"{row['completion_rate']:.1f}"
                if best_completion is not None:
                    comp_str = bold_if_best(comp_str, abs(row['completion_rate'] - best_completion) < 0.01)
            else:
                comp_str = "--"
            
            # Domain name only on first row
            domain_str = task_group if j == 0 else ""
            
            # Format agent name
            agent_display = format_agent_name(row['agent'])
            
            # Format rank if available
            if has_rank and row.get('rank_mean') is not None and row.get('rank_sem') is not None:
                rank_str = format_value_with_sem(row['rank_mean'], row['rank_sem'])
                if best_rank is not None:
                    rank_str = bold_if_best(rank_str, abs(row['rank_mean'] - best_rank) < 0.01)
                latex.append(f"{domain_str} & {agent_display} & {lb_str} & {rank_str} & {above_median_str} & {any_medal_str} & {comp_str} \\\\")
            else:
                latex.append(f"{domain_str} & {agent_display} & {lb_str} & {above_median_str} & {any_medal_str} & {comp_str} \\\\")
    
    latex.append("\\bottomrule")
    latex.append("\\end{tabular}%")
    latex.append("}")
    latex.append("\\end{table}")
    
    return '\n'.join(latex)

def generate_overall_latex_table(stats_df: pd.DataFrame) -> str:
    """Generate LaTeX table for overall results."""
    
    latex = []
    latex.append("\\begin{table}[htbp]")
    latex.append("\\centering")
    latex.append("\\caption{BioML-bench Overall Results}")
    latex.append("\\label{tab:bioml_results_overall}")
    latex.append("\\resizebox{\\textwidth}{!}{%")
    
    # Check if rank data is available
    has_rank = 'rank_mean' in stats_df.columns and stats_df['rank_mean'].notna().any()
    
    if has_rank:
        latex.append("\\begin{tabular}{lccccc}")
        latex.append("\\toprule")
        latex.append("Agent & Leaderboard Percentile & Mean Rank & Above Median (\\%) & Any Medal (\\%) & Completion Rate (\\%) \\\\")
    else:
        latex.append("\\begin{tabular}{lcccc}")
        latex.append("\\toprule")
        latex.append("Agent & Leaderboard Percentile & Above Median (\\%) & Any Medal (\\%) & Completion Rate (\\%) \\\\")
    latex.append("\\midrule")
    
    # Find best values for bolding
    best_lb = stats_df['leaderboard_mean'].max()
    best_above_median = stats_df['above_median_mean'].max()
    best_medal = stats_df['any_medal_mean'].max()
    
    # For rank, lower is better
    best_rank = None
    if has_rank and 'rank_mean' in stats_df.columns:
        rank_values = stats_df['rank_mean'].dropna()
        if len(rank_values) > 0:
            best_rank = rank_values.min()
    
    # For completion rate, exclude dummy agent when finding best
    non_dummy_data = stats_df[stats_df['agent'].str.lower() != 'dummy']
    best_completion = non_dummy_data['completion_rate'].max() if len(non_dummy_data) > 0 and 'completion_rate' in non_dummy_data.columns else None
    
    for _, row in stats_df.iterrows():
        # Format values
        lb_str = format_value_with_sem(row['leaderboard_mean'], row['leaderboard_sem'])
        above_median_str = format_value_with_sem(row['above_median_mean'], row['above_median_sem'], is_percentage=True)
        any_medal_str = format_value_with_sem(row['any_medal_mean'], row['any_medal_sem'], is_percentage=True)
        
        # Bold if best
        lb_str = bold_if_best(lb_str, abs(row['leaderboard_mean'] - best_lb) < 0.01)
        above_median_str = bold_if_best(above_median_str, abs(row['above_median_mean'] - best_above_median) < 0.01)
        any_medal_str = bold_if_best(any_medal_str, abs(row['any_medal_mean'] - best_medal) < 0.01)
        
        # Completion rate - special handling for dummy
        if row['agent'].lower() == 'dummy':
            comp_str = "NA"
        elif row['completion_rate'] is not None:
            comp_str = f"{row['completion_rate']:.1f}"
            if best_completion is not None:
                comp_str = bold_if_best(comp_str, abs(row['completion_rate'] - best_completion) < 0.01)
        else:
            comp_str = "--"
        
        # Format agent name
        agent_display = format_agent_name(row['agent'])
        
        # Format with rank if available
        if has_rank and row.get('rank_mean') is not None and row.get('rank_sem') is not None:
            rank_str = format_value_with_sem(row['rank_mean'], row['rank_sem'])
            if best_rank is not None:
                rank_str = bold_if_best(rank_str, abs(row['rank_mean'] - best_rank) < 0.01)
            latex.append(f"{agent_display} & {lb_str} & {rank_str} & {above_median_str} & {any_medal_str} & {comp_str} \\\\")
        else:
            latex.append(f"{agent_display} & {lb_str} & {above_median_str} & {any_medal_str} & {comp_str} \\\\")
    
    latex.append("\\bottomrule")
    latex.append("\\end{tabular}%")
    latex.append("}")
    latex.append("\\end{table}")
    
    return '\n'.join(latex)

def main():
    parser = argparse.ArgumentParser(
        description="Generate LaTeX summary tables for BioML-bench results"
    )
    parser.add_argument(
        "--data-dir",
        required=True,
        help="Directory containing all_replicates.csv and completion_rates.csv"
    )
    parser.add_argument(
        "--output-dir",
        default="latex_tables",
        help="Directory to save LaTeX tables (default: latex_tables)"
    )
    
    args = parser.parse_args()
    
    # Construct file paths from data directory
    data_dir = Path(args.data_dir)
    replicates_file = data_dir / 'all_replicates.csv'
    completion_file = data_dir / 'completion_rates.csv'
    
    # Check files exist
    if not replicates_file.exists():
        raise FileNotFoundError(f"all_replicates.csv not found in {data_dir}")
    if not completion_file.exists():
        log(f"⚠️  completion_rates.csv not found in {data_dir} - will skip completion rate analysis")
        completion_df = None
    else:
        log(f"Loading completion rates from {completion_file}")
        completion_df = pd.read_csv(completion_file)
    
    # Load data
    log(f"Loading replicate data from {replicates_file}")
    df = pd.read_csv(replicates_file)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Calculate statistics
    log("Calculating subdomain statistics...")
    subdomain_stats = calculate_subdomain_statistics(df, completion_df)
    
    log("Calculating overall statistics...")
    overall_stats = calculate_overall_statistics(df, completion_df)
    
    # Generate LaTeX tables
    log("Generating LaTeX tables...")
    
    subdomain_table = generate_subdomain_latex_table(subdomain_stats)
    overall_table = generate_overall_latex_table(overall_stats)
    
    # Save tables
    with open(output_dir / "table_by_domain.tex", 'w') as f:
        f.write(subdomain_table)
    
    with open(output_dir / "table_overall.tex", 'w') as f:
        f.write(overall_table)
    
    # Also save as standalone LaTeX document for easy compilation
    standalone_doc = f"""\\documentclass{{article}}
\\usepackage{{booktabs}}
\\usepackage{{graphicx}}
\\begin{{document}}

{subdomain_table}

\\vspace{{1cm}}

{overall_table}

\\end{{document}}
"""
    
    with open(output_dir / "tables_standalone.tex", 'w') as f:
        f.write(standalone_doc)
    
    # Save CSV versions for inspection
    subdomain_stats.to_csv(output_dir / "subdomain_statistics.csv", index=False)
    overall_stats.to_csv(output_dir / "overall_statistics.csv", index=False)
    
    log("LaTeX tables generated successfully!")
    log(f"📊 Tables saved to: {output_dir}")
    log("  - table_by_domain.tex: Results broken down by domain")
    log("  - table_overall.tex: Overall results across all tasks")
    log("  - tables_standalone.tex: Complete LaTeX document with both tables")
    log("  - subdomain_statistics.csv: Raw statistics by domain")
    log("  - overall_statistics.csv: Raw overall statistics")

if __name__ == "__main__":
    main() 