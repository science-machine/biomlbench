#!/usr/bin/env python3
"""
BioML-bench System Health Monitor
Continuous monitoring of system health and early warning for potential issues
"""

import os
import sys
import time
import json
import psutil
import logging
import subprocess
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class SystemHealthMonitor:
    def __init__(self, config_file: Optional[str] = None):
        self.config = self.load_config(config_file)
        self.last_check = {}
    
    def load_config(self, config_file: Optional[str]) -> Dict:
        """Load monitoring configuration"""
        default_config = {
            "check_interval": 300,  # 5 minutes
            "disk_threshold": 85,   # 85% disk usage warning
            "memory_threshold": 90, # 90% memory usage warning
            "log_retention_days": 30,
            "max_container_runtime": 7200,  # 2 hours max for agent runs
        }
        
        if config_file and Path(config_file).exists():
            try:
                with open(config_file, 'r') as f:
                    user_config = json.load(f)
                    default_config.update(user_config)
            except Exception as e:
                logger.warning(f"Failed to load config {config_file}: {e}")
        
        return default_config
    
    def check_system_resources(self) -> Dict:
        """Check system resource usage"""
        try:
            # CPU usage
            cpu_percent = psutil.cpu_percent(interval=1)
            
            # Memory usage
            memory = psutil.virtual_memory()
            memory_percent = memory.percent
            
            # Disk usage
            disk = psutil.disk_usage('/')
            disk_percent = (disk.used / disk.total) * 100
            
            return {
                'timestamp': datetime.now().isoformat(),
                'cpu_percent': cpu_percent,
                'memory_percent': memory_percent,
                'memory_available_gb': memory.available / (1024**3),
                'disk_percent': disk_percent,
                'disk_free_gb': disk.free / (1024**3),
                'status': 'healthy'
            }
            
        except Exception as e:
            logger.error(f"Failed to check system resources: {e}")
            return {'status': 'error', 'error': str(e)}
    
    def check_environment_integrity(self) -> Dict:
        """Check if the environment can still run basic operations"""
        try:
            # Test basic Python environment
            test_script = '''
import pandas as pd
import numpy as np
import sklearn
print("Basic environment: OK")
'''
            
            result = subprocess.run([
                sys.executable, '-c', test_script
            ], capture_output=True, text=True, timeout=30)
            
            if result.returncode != 0:
                return {
                    'status': 'degraded',
                    'basic_env_test': 'failed',
                    'error': result.stderr
                }
            
            return {
                'status': 'healthy',
                'basic_env_test': 'passed',
                'last_check': datetime.now().isoformat()
            }
            
        except Exception as e:
            return {'status': 'error', 'error': str(e)}
    
    def generate_health_report(self) -> Dict:
        """Generate comprehensive health report"""
        logger.info("Generating health report...")
        
        report = {
            'timestamp': datetime.now().isoformat(),
            'system_resources': self.check_system_resources(),
            'environment_integrity': self.check_environment_integrity()
        }
        
        # Determine overall health status
        component_statuses = [
            report['system_resources'].get('status'),
            report['environment_integrity'].get('status')
        ]
        
        if 'error' in component_statuses:
            overall_status = 'critical'
        elif 'degraded' in component_statuses:
            overall_status = 'degraded'
        else:
            overall_status = 'healthy'
        
        report['overall_status'] = overall_status
        
        # Generate warnings and recommendations
        warnings = []
        recommendations = []
        
        # Resource warnings
        if report['system_resources'].get('memory_percent', 0) > self.config['memory_threshold']:
            warnings.append(f"High memory usage: {report['system_resources']['memory_percent']:.1f}%")
            recommendations.append("Consider freeing memory or increasing available RAM")
        
        if report['system_resources'].get('disk_percent', 0) > self.config['disk_threshold']:
            warnings.append(f"High disk usage: {report['system_resources']['disk_percent']:.1f}%")
            recommendations.append("Clean up old runs and logs")
        
        report['warnings'] = warnings
        report['recommendations'] = recommendations
        
        return report

def main():
    """Main monitoring function"""
    monitor = SystemHealthMonitor()
    
    # Single check mode
    report = monitor.generate_health_report()
    
    print(f"Overall Status: {report['overall_status'].upper()}")
    
    if report['warnings']:
        print("\n⚠️ Warnings:")
        for warning in report['warnings']:
            print(f"  - {warning}")
    
    if report['recommendations']:
        print("\n💡 Recommendations:")
        for rec in report['recommendations']:
            print(f"  - {rec}")
    
    # Exit with appropriate code
    sys.exit(0 if report['overall_status'] == 'healthy' else 1)

if __name__ == "__main__":
    main() 