#!/usr/bin/env python3
"""
BioML-bench Failure Modes Visualization

Creates a stacked bar chart showing the proportional failure modes for each agent.
Each bar represents an agent, and each segment shows what percentage of that 
agent's total failures are of each type.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set publication-ready style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("Set2")

def load_and_process_failure_data(csv_path):
    """Load failure modes data and calculate proportions per agent."""
    
    # Read the CSV
    df = pd.read_csv(csv_path)
    
    # Get failure mode columns (everything after 'Agent')
    failure_mode_cols = df.columns[2:].tolist()
    
    print(f"Failure modes found: {failure_mode_cols}")
    
    # Fill NaN values with 0 and convert to numeric
    for col in failure_mode_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
    
    # Group by agent and sum up failure counts
    agent_totals = df.groupby('Agent')[failure_mode_cols].sum()
    
    # Calculate total failures per agent
    agent_totals['total_failures'] = agent_totals.sum(axis=1)
    
    # Calculate proportions (what % of each agent's failures are each type)
    proportions_df = agent_totals[failure_mode_cols].div(agent_totals['total_failures'], axis=0) * 100
    
    # Replace NaN (division by 0) with 0
    proportions_df = proportions_df.fillna(0)
    
    print("\nFailure counts by agent:")
    print(agent_totals)
    print("\nFailure proportions by agent (%):")
    print(proportions_df.round(1))
    
    return proportions_df, agent_totals, failure_mode_cols

def create_stacked_bar_chart(proportions_df, agent_totals, failure_mode_cols, output_path="failure_modes_stacked_bar.png"):
    """Create a stacked bar chart showing failure mode proportions for each agent."""
    
    # Create shorter labels for better readability
    short_labels = []
    for label in failure_mode_cols:
        if 'exception' in label.lower():
            short_labels.append('Code Exceptions')
        elif 'finish running' in label.lower():
            short_labels.append('Failed to Finish Running')
        elif 'submit answers' in label.lower() or 'format' in label.lower():
            short_labels.append('Submission Error')
        elif 'time limit' in label.lower():
            short_labels.append('Time Limit')
        elif 'resources' in label.lower():
            short_labels.append('Resource Exhaustion')
        else:
            # Fallback: use first 3 words
            short_labels.append(' '.join(label.split()[:3]))
    
    # Create a copy with shorter column names for plotting
    plot_df = proportions_df.copy()
    plot_df.columns = short_labels
    
    # Create the figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Define medium-dark pastel colors for each failure mode
    failure_colors = ['#ff6666', '#ffb366', '#ffff66', '#66b3ff', '#b366ff']
    
    # Create stacked bar chart with visible borders
    plot_df.plot(kind='bar', stacked=True, ax=ax, color=failure_colors, 
                 width=0.7, edgecolor='black', linewidth=0.8)
    
    # Customize the plot with larger fonts
    ax.set_title('Agent Failure Modes', 
                 fontsize=28, fontweight='bold', pad=25)
    
    ax.set_xlabel('Agent', fontsize=22, fontweight='bold')
    ax.set_ylabel('Proportion of Failures (%)', fontsize=22, fontweight='bold')
    
    # Format x-axis with agent names and failure counts
    x_labels = []
    for agent in plot_df.index:
        total_failures = int(agent_totals.loc[agent, 'total_failures'])
        x_labels.append(f"{agent.upper()}")
    
    ax.set_xticklabels(x_labels, rotation=0, fontsize=18, fontweight='bold')
    
    # Format y-axis
    ax.set_ylim(0, 100)
    ax.tick_params(axis='y', labelsize=16)
    
    # Add grid
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_axisbelow(True)
    
    # Customize legend with larger fonts
    ax.legend(title='Failure Mode', bbox_to_anchor=(1.05, 1), loc='upper left',
              fontsize=16, title_fontsize=18, frameon=False)
    
    # Add percentage labels on bars (only if segment is large enough)
    for container in ax.containers:
        # Add labels only for segments > 5%
        labels = [f'{v:.0f}%' if v > 5 else '' for v in container.datavalues]
        ax.bar_label(container, labels=labels, label_type='center', 
                    fontsize=14, color='black')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\nStacked bar chart saved to: {output_path}")

def create_summary_stats(proportions_df, agent_totals, failure_mode_cols):
    """Create summary statistics table."""
    
    print("\n" + "="*80)
    print("FAILURE MODE ANALYSIS SUMMARY")
    print("="*80)
    
    # Overall statistics
    total_failures = agent_totals['total_failures'].sum()
    print(f"\nTotal failures across all agents: {total_failures}")
    
    print(f"\nFailures per agent:")
    for agent in agent_totals.index:
        count = agent_totals.loc[agent, 'total_failures']
        pct = (count / total_failures) * 100
        print(f"  {agent.upper()}: {count:3d} failures ({pct:.1f}%)")
    
    # Most common failure modes overall
    print(f"\nMost common failure modes (absolute counts):")
    mode_totals = agent_totals[failure_mode_cols].sum().sort_values(ascending=False)
    for mode, count in mode_totals.items():
        pct = (count / total_failures) * 100
        print(f"  {mode}: {count:3d} ({pct:.1f}%)")
    
    # Agent specializations (highest proportion failure modes)
    print(f"\nAgent failure specializations (highest proportion):")
    for agent in proportions_df.index:
        max_mode = proportions_df.loc[agent].idxmax()
        max_pct = proportions_df.loc[agent, max_mode]
        print(f"  {agent.upper()}: {max_pct:.1f}% of failures are '{max_mode}'")

def main():
    """Main function to run the failure mode analysis."""
    
    # Input and output paths
    csv_path = Path("scripts/analysis_outputs/Failure Modes for BioMLBench - Sheet3 (1).csv")
    output_path = Path("scripts/analysis_outputs/failure_modes_stacked_bar.png")
    
    if not csv_path.exists():
        print(f"Error: Could not find input file: {csv_path}")
        print("Please make sure the CSV file exists in the analysis_outputs directory.")
        return
    
    # Load and process the data
    print("Loading failure modes data...")
    proportions_df, agent_totals, failure_mode_cols = load_and_process_failure_data(csv_path)
    
    # Create the stacked bar chart
    print("\nCreating stacked bar chart...")
    create_stacked_bar_chart(proportions_df, agent_totals, failure_mode_cols, output_path)
    
    # Print summary statistics
    create_summary_stats(proportions_df, agent_totals, failure_mode_cols)
    
    print("\n✅ Failure mode analysis complete!")

if __name__ == "__main__":
    main() 