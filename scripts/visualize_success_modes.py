#!/usr/bin/env python3
"""
BioML-bench Success Modes Visualization

Creates visualizations showing the strategies agents used to achieve successful runs
(above median human performance), including:
1. Agent-wise strategy usage heatmap 
2. Model architecture distribution stacked bar chart
3. Summary statistics tables in LaTeX format
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
from collections import Counter

# Set publication-ready style with better color contrast
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("tab10")  # More distinguishable colors

def load_and_process_success_data(csv_path):
    """Load success modes data and process strategy columns."""
    
    # Read the CSV
    df = pd.read_csv(csv_path)
    
    # Define strategy columns (everything between Task and Model columns)
    strategy_cols = ['Deep learning', 'Pretraing model', 'Exotic', 'Tuning', 'CV', 
                    'Adversarial training', 'Data augmentation', 'Stacked', 
                    'FeatureSelection', 'FeatureEngineering', 'Model select']
    
    print(f"Strategy columns found: {strategy_cols}")
    print(f"Total runs in dataset: {len(df)}")
    
    # Fill NaN values with 0 and convert to numeric
    for col in strategy_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
    
    # Calculate strategy usage by agent
    agent_strategy_usage = df.groupby('Agent')[strategy_cols].agg(['sum', 'count'])
    
    # Calculate proportions (what % of each agent's runs used each strategy)
    proportions_df = pd.DataFrame()
    for col in strategy_cols:
        proportions_df[col] = agent_strategy_usage[(col, 'sum')] / agent_strategy_usage[(col, 'count')] * 100
    
    # Get total runs per agent
    runs_per_agent = df.groupby('Agent').size()
    
    print(f"\nRuns per agent:")
    for agent, count in runs_per_agent.items():
        print(f"  {agent}: {count} successful runs")
    
    print(f"\nStrategy usage proportions by agent (%):")
    print(proportions_df.round(1))
    
    return proportions_df, runs_per_agent, strategy_cols, df

def clean_model_data(df):
    """Clean and standardize the Model column data."""
    
    # Create a copy to work with
    df_clean = df.copy()
    
    # Fill empty model entries
    df_clean['Model'] = df_clean['Model'].fillna('Unknown').astype(str)
    
    # Standardize model names
    def standardize_model(model_str):
        if pd.isna(model_str) or model_str.strip() == '' or model_str == 'Unknown':
            return 'Unknown'
        
        model_str = model_str.strip()
        
        # Check for ensemble/stacked models (multiple models separated by commas)
        if ',' in model_str:
            return 'Stacked Models'
        
        # Standardize common model names
        model_lower = model_str.lower()
        
        if 'random forest' in model_lower or model_lower in ['rf', 'randomforestregressor', 'randomforestclassifier']:
            return 'Random Forest'
        elif 'xgb' in model_lower or 'xgboost' in model_lower:
            return 'XGBoost'
        elif 'lgbm' in model_lower or 'lightgbm' in model_lower:
            return 'LightGBM'
        elif 'gradient boost' in model_lower or 'gbm' in model_lower or 'gradientboostingregressor' in model_lower:
            return 'GBM'
        elif 'catboost' in model_lower:
            return 'CatBoost'
        elif 'histgradient' in model_lower:
            return 'HistGradientBoosting'
        elif model_lower in ['lr', 'logistic'] or 'logisticregression' in model_lower:
            return 'Logistic Regression'
        elif 'ridge' in model_lower:
            return 'Ridge'
        elif 'lasso' in model_lower:
            return 'Lasso'
        elif 'elastic' in model_lower:
            return 'Elastic Net'
        elif 'svm' in model_lower or 'svr' in model_lower:
            return 'SVM'
        elif 'knn' in model_lower:
            return 'KNN'
        elif 'resnet' in model_lower:
            return 'ResNet'
        elif 'neural' in model_lower or 'nn' in model_lower or 'mlp' in model_lower:
            return 'Neural Network'
        elif 'baseline' in model_lower:
            return 'Baseline'
        elif 'kernelridge' in model_lower:
            return 'Kernel Ridge'
        else:
            return 'Other'
    
    df_clean['Model_Clean'] = df_clean['Model'].apply(standardize_model)
    
    print(f"\nModel standardization:")
    model_counts = df_clean['Model_Clean'].value_counts()
    for model, count in model_counts.items():
        print(f"  {model}: {count} runs")
    
    return df_clean

def create_strategy_heatmap(proportions_df, runs_per_agent, output_path="strategy_heatmap.png"):
    """Create a heatmap showing strategy usage proportions by agent."""
    
    # Create shorter labels for better readability
    short_labels = {
        'Deep learning': 'Deep Learning',
        'Pretraing model': 'Pretrained Model',
        'Exotic': 'Exotic Methods',
        'Tuning': 'Hyperparameter Tuning',
        'CV': 'Cross Validation',
        'Adversarial training': 'Adversarial Training',
        'Data augmentation': 'Data Augmentation',
        'Stacked': 'Model Stacking',
        'FeatureSelection': 'Feature Selection',
        'FeatureEngineering': 'Feature Engineering',
        'Model select': 'Model Selection'
    }
    
    # Rename columns
    plot_df = proportions_df.rename(columns=short_labels)
    
    # Remove columns with all zeros (unused strategies)
    plot_df = plot_df.loc[:, (plot_df != 0).any(axis=0)]
    
    # Sort by sum values for better visualization
    # Sort columns (strategies) by total usage across all agents
    col_sums = plot_df.sum(axis=0).sort_values(ascending=False)
    plot_df = plot_df[col_sums.index]
    
    # Sort rows (agents) by total strategy usage
    row_sums = plot_df.sum(axis=1).sort_values(ascending=False)
    plot_df = plot_df.loc[row_sums.index]
    
    # Create the figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Create heatmap with better color contrast
    sns.heatmap(plot_df, annot=True, fmt='.1f', cmap='viridis', 
                cbar_kws={'label': 'Usage Percentage (%)'}, 
                annot_kws={'fontsize': 14, 'fontweight': 'bold'},
                ax=ax, linewidths=0.5, linecolor='white')
    
    # Customize the plot
    ax.set_title('Agent Strategies Usage in Above-Median Runs', 
                 fontsize=24, fontweight='bold', pad=25)
    
    ax.set_xlabel('Strategy', fontsize=18, fontweight='bold')
    ax.set_ylabel('Agent', fontsize=18, fontweight='bold')
    
    # Format agent names with run counts
    y_labels = [f"{agent.upper()}\n(n={runs_per_agent[agent]})" for agent in plot_df.index]
    ax.set_yticklabels(y_labels, rotation=0, fontsize=16, fontweight='bold')
    ax.set_xticklabels(plot_df.columns, rotation=45, ha='right', fontsize=17)
    
    # Color bar formatting
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('Usage Percentage (%)', fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\nStrategy heatmap saved to: {output_path}")

def create_model_stacked_bar_chart(df_clean, output_path="model_stacked_bar.png"):
    """Create a stacked bar chart showing model usage proportions by agent."""
    
    # Calculate model usage proportions by agent
    model_counts = df_clean.groupby(['Agent', 'Model_Clean']).size().unstack(fill_value=0)
    model_proportions = model_counts.div(model_counts.sum(axis=1), axis=0) * 100
    
    # Get only the most common models to avoid clutter
    overall_model_usage = df_clean['Model_Clean'].value_counts()
    top_models = overall_model_usage.head(8).index.tolist()
    
    # Keep only top models, group rest as "Other"
    plot_df = model_proportions[top_models].copy()
    if len(model_proportions.columns) > len(top_models):
        other_cols = [col for col in model_proportions.columns if col not in top_models]
        plot_df['Other'] = model_proportions[other_cols].sum(axis=1)
    
    # Create the figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Define more distinguishable colors for different model types
    distinct_colors = [
        '#1f77b4',  # Blue
        '#ff7f0e',  # Orange  
        '#2ca02c',  # Green
        '#d62728',  # Red
        '#9467bd',  # Purple
        '#8c564b',  # Brown
        '#e377c2',  # Pink
        '#7f7f7f',  # Gray
        '#bcbd22',  # Olive
        '#17becf'   # Cyan
    ]
    colors = distinct_colors[:len(plot_df.columns)]
    
    # Create stacked bar chart
    plot_df.plot(kind='bar', stacked=True, ax=ax, color=colors, 
                 width=0.7, edgecolor='black', linewidth=0.8)
    
    # Customize the plot
    ax.set_title('Model Architecture Usage by Agent', 
                 fontsize=24, fontweight='bold', pad=25)
    
    ax.set_xlabel('Agent', fontsize=18, fontweight='bold')
    ax.set_ylabel('Proportion of Above-Median Runs (%)', fontsize=18, fontweight='bold')
    
    # Format x-axis with agent names
    x_labels = [agent.upper() for agent in plot_df.index]
    ax.set_xticklabels(x_labels, rotation=0, fontsize=16, fontweight='bold')
    
    # Format y-axis
    ax.set_ylim(0, 100)
    ax.tick_params(axis='y', labelsize=14)
    
    # Add grid
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_axisbelow(True)
    
    # Customize legend
    ax.legend(title='Model Type', bbox_to_anchor=(1.05, 1), loc='upper left',
              fontsize=12, title_fontsize=14, frameon=False)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\nModel architecture chart saved to: {output_path}")
    
    return model_proportions

def create_latex_tables(proportions_df, runs_per_agent, model_proportions, df, output_dir="analysis_outputs"):
    """Create LaTeX tables following the paper style."""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Table 1: Strategy Usage Summary
    strategy_table = """\\begin{table}[htbp]
\\centering
\\caption{Strategy Usage in Successful Runs by Agent}
\\label{tab:strategy_usage}
\\resizebox{\\textwidth}{!}{%
\\begin{tabular}{l""" + "c" * len(proportions_df.columns) + """c}
\\toprule
Agent & """ + " & ".join([col.replace('_', ' ').replace('Feature', 'Feat.') for col in proportions_df.columns]) + """ & Total Runs \\\\
\\midrule
"""
    
    for agent in proportions_df.index:
        values = [f"{proportions_df.loc[agent, col]:.0f}" for col in proportions_df.columns]
        total_runs = runs_per_agent[agent]
        strategy_table += f"{agent.upper()} & " + " & ".join(values) + f" & {total_runs} \\\\\n"
    
    strategy_table += """\\bottomrule
\\end{tabular}%
}
\\end{table}"""
    
    # Save strategy table
    with open(output_dir / "table_strategy_usage.tex", "w") as f:
        f.write(strategy_table)
    
    # Table 2: Model Usage Summary (top models only)
    top_models = model_proportions.sum().nlargest(6).index.tolist()
    model_table = """\\begin{table}[htbp]
\\centering
\\caption{Model Architecture Usage in Successful Runs by Agent (\\%)}
\\label{tab:model_usage}
\\resizebox{\\textwidth}{!}{%
\\begin{tabular}{l""" + "c" * len(top_models) + """}
\\toprule
Agent & """ + " & ".join(top_models) + """ \\\\
\\midrule
"""
    
    for agent in model_proportions.index:
        values = [f"{model_proportions.loc[agent, model]:.0f}" if model in model_proportions.columns else "0" 
                 for model in top_models]
        model_table += f"{agent.upper()} & " + " & ".join(values) + " \\\\\n"
    
    model_table += """\\bottomrule
\\end{tabular}%
}
\\end{table}"""
    
    # Save model table
    with open(output_dir / "table_model_usage.tex", "w") as f:
        f.write(model_table)
    
    print(f"\nLaTeX tables saved to:")
    print(f"  {output_dir}/table_strategy_usage.tex")
    print(f"  {output_dir}/table_model_usage.tex")

def create_summary_stats(proportions_df, runs_per_agent, df, strategy_cols):
    """Create summary statistics for success modes."""
    
    print("\n" + "="*80)
    print("SUCCESS MODE ANALYSIS SUMMARY")
    print("="*80)
    
    # Overall statistics
    total_successful_runs = len(df)
    print(f"\nTotal successful runs across all agents: {total_successful_runs}")
    
    print(f"\nSuccessful runs per agent:")
    for agent in runs_per_agent.index:
        count = runs_per_agent[agent]
        pct = (count / total_successful_runs) * 100
        print(f"  {agent.upper()}: {count:3d} runs ({pct:.1f}%)")
    
    # Most commonly used strategies overall
    strategy_totals = df[strategy_cols].sum().sort_values(ascending=False)
    print(f"\nMost commonly used strategies (absolute counts):")
    for strategy, count in strategy_totals.items():
        pct = (count / total_successful_runs) * 100
        print(f"  {strategy}: {int(count):3d} uses ({pct:.1f}%)")
    
    # Agent strategy preferences (highest usage strategies)
    print(f"\nAgent strategy preferences (highest usage %):")
    for agent in proportions_df.index:
        max_strategy = proportions_df.loc[agent].idxmax()
        max_pct = proportions_df.loc[agent, max_strategy]
        print(f"  {agent.upper()}: {max_pct:.1f}% usage of '{max_strategy}'")
    
    # Performance correlation (basic)
    print(f"\nAverage performance by agent:")
    perf_by_agent = df.groupby('Agent')['Leaderboard Score'].mean()
    for agent, score in perf_by_agent.items():
        print(f"  {agent.upper()}: {score:.1f} average score")

def main():
    """Main function to run the success mode analysis."""
    
    # Input and output paths
    
    csv_path = Path("scripts/analysis_outputs/Success modes BioML-Bench - Sheet1.csv")
    output_dir = Path("scripts/analysis_outputs")
    output_dir.mkdir(exist_ok=True)
    
    if not csv_path.exists():
        print(f"Error: Could not find input file: {csv_path}")
        print("Please make sure the CSV file exists.")
        return
    
    # Load and process the data
    print("Loading success modes data...")
    proportions_df, runs_per_agent, strategy_cols, df = load_and_process_success_data(csv_path)
    
    # Clean model data
    print("\nCleaning model data...")
    df_clean = clean_model_data(df)
    
    # Create visualizations
    print("\nCreating strategy heatmap...")
    create_strategy_heatmap(proportions_df, runs_per_agent, output_dir / "strategy_heatmap.png")
    
    print("\nCreating model architecture chart...")
    model_proportions = create_model_stacked_bar_chart(df_clean, output_dir / "model_stacked_bar.png")
    
    # Create LaTeX tables
    print("\nGenerating LaTeX tables...")
    create_latex_tables(proportions_df, runs_per_agent, model_proportions, df_clean)
    
    # Print summary statistics
    create_summary_stats(proportions_df, runs_per_agent, df, strategy_cols)
    
    print("\n✅ Success mode analysis complete!")

if __name__ == "__main__":
    main() 