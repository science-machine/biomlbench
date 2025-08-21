#!/usr/bin/env python3
"""
Generate production jobs file for BioML-bench.

This script creates a complete jobs file with all agents (aide, biomni, mlagentbench, stella)
for all production tasks across all biomedical domains.
"""

from pathlib import Path
import os

def get_biomlbench_root():
    """Get the biomlbench project root directory."""
    return Path(__file__).parent.parent

def verify_task_exists(task_id):
    """Verify that a task exists in the tasks directory."""
    root = get_biomlbench_root()
    task_path = root / "biomlbench" / "tasks" / task_id
    config_path = task_path / "config.yaml"
    return config_path.exists()

def main():
    # All agents to test
    agents = ["aide", "biomni", "mlagentbench", "stella"]
    
    # Production tasks organized by domain
    production_tasks = {
        "Drug Discovery and Development": [
            "polarishub/polaris-pkis2-egfr-wt-c-1",
            "polarishub/polaris-adme-fang-hclint-1", 
            "polarishub/polaris-adme-fang-hppb-1",
            "polarishub/polaris-adme-fang-solu-1",
            "polarishub/polaris-adme-fang-r-1",
            "polarishub/tdcommons-cyp2d6-substrate-carbonmangels",
            "polarishub/tdcommons-lipophilicity-astrazeneca",
            "polarishub/tdcommons-herg",
            "polarishub/tdcommons-bbb-martins",
        ],
        "Biomedical Image Analysis": [
            "kaggle/osic-pulmonary-fibrosis-progression",
            "kaggle/histopathologic-cancer-detection", 
            "kaggle/rsna-miccai-brain-tumor-radiogenomic-classification",
            "kaggle/uw-madison-gi-tract-image-segmentation",
        ],
        "Protein Engineering": [
            "proteingym-dms/SPIKE_SARS2_Starr_2020_binding",  # single-substitution
            "proteingym-dms/SBI_STAAM_Tsuboyama_2023_2JVG",   # single-substitution  
            "proteingym-dms/PSAE_PICP2_Tsuboyama_2023_1PSE",  # multi-substitution
            "proteingym-dms/CBX4_HUMAN_Tsuboyama_2023_2K28",  # multi-substitution
            "proteingym-dms/Q8EG35_SHEON_Campbell_2022_indels",  # indels
            "proteingym-dms/CSN4_MOUSE_Tsuboyama_2023_1UFM_indels",  # indels
        ],
        "Single Cell Omics": [
            "manual/open-problems-predict-modality",          # multimodal prediction
            "manual/open-problems-single-cell-perturbations", # single cell perturbations
            "manual/open-problems-cell-cell-communication-ligand-target",  # cell-cell communication ligand::target  
            "manual/open-problems-spatially-variable-genes",  # spatially variable genes
        ]
    }
    
    # Verify all tasks exist
    print("🔍 Verifying task IDs...")
    missing_tasks = []
    total_tasks = 0
    
    for domain, tasks in production_tasks.items():
        print(f"\n📂 {domain}:")
        for task_id in tasks:
            total_tasks += 1
            if verify_task_exists(task_id):
                print(f"  ✅ {task_id}")
            else:
                print(f"  ❌ {task_id} - NOT FOUND")
                missing_tasks.append(task_id)
    
    if missing_tasks:
        print(f"\n❌ Found {len(missing_tasks)} missing tasks:")
        for task in missing_tasks:
            print(f"  - {task}")
        print("\nPlease fix task IDs before generating jobs file.")
        return 1
    
    # Generate jobs file
    output_file = get_biomlbench_root() / "production-jobs.txt"
    total_jobs = 0
    
    print(f"\n📝 Generating production jobs file: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write("# BioML-bench Production Jobs - Complete Task List\n")
        f.write("# Format: agent,task_id\n")
        f.write("# Run from project root: python scripts/gcp-deploy/deploy.py --jobs production-jobs.txt --concurrent 15\n")
        f.write(f"# Generated with {len(agents)} agents × {total_tasks} tasks = {len(agents) * total_tasks} total jobs\n\n")
        
        for domain, tasks in production_tasks.items():
            f.write(f"# {'=' * 50}\n")
            f.write(f"# {domain}\n") 
            f.write(f"# {'=' * 50}\n\n")
            
            for task_id in tasks:
                for agent in agents:
                    f.write(f"{agent},{task_id}\n")
                    total_jobs += 1
                f.write("\n")  # Blank line after each task
    
    print(f"✅ Generated {total_jobs} jobs ({len(agents)} agents × {total_tasks} tasks)")
    print(f"📊 Breakdown:")
    for domain, tasks in production_tasks.items():
        jobs_in_domain = len(tasks) * len(agents)
        print(f"  - {domain}: {jobs_in_domain} jobs ({len(tasks)} tasks × {len(agents)} agents)")
    
    # Estimate runtime and cost
    avg_runtime_minutes = 25  # 15-30 minutes per job
    concurrent_vms = 15
    
    estimated_hours = (total_jobs * avg_runtime_minutes) / (concurrent_vms * 60)
    estimated_cost = total_jobs * 0.60  # ~$0.60 per job
    
    print(f"\n💰 Estimates with {concurrent_vms} concurrent VMs:")
    print(f"  - Runtime: {estimated_hours:.1f} hours")
    print(f"  - Cost: ${estimated_cost:.0f}")
    
    print(f"\n🚀 To run: python scripts/gcp-deploy/deploy.py --jobs production-jobs.txt --concurrent 15")
    
    return 0

if __name__ == "__main__":
    exit(main()) 