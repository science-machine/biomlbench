"""
S3 upload utilities for biomlbench artifacts.

This module provides functionality to automatically upload run artifacts and grading results
to S3 buckets for backup and persistence during long-running benchmarks.
"""

import gzip
import json
import os
import tarfile
import tempfile
import time
from pathlib import Path
from typing import Optional, Dict, Any
import logging

from biomlbench.utils import get_logger

logger = get_logger(__name__)

import boto3
from botocore.exceptions import BotoCoreError, ClientError


class S3Config:
    """Configuration for S3 uploads."""
    
    def __init__(
        self,
        bucket_name: Optional[str] = None,
        prefix: Optional[str] = None,
        region: Optional[str] = None,
        aws_access_key_id: Optional[str] = None,
        aws_secret_access_key: Optional[str] = None,
        compress: bool = True,
        enabled: bool = False,
    ):
        self.bucket_name = bucket_name or os.environ.get("BIOMLBENCH_S3_BUCKET")
        self.prefix = prefix or ""
            
        self.region = region or os.environ.get("AWS_DEFAULT_REGION", "us-east-1")
        self.aws_access_key_id = aws_access_key_id or os.environ.get("AWS_ACCESS_KEY_ID")
        self.aws_secret_access_key = aws_secret_access_key or os.environ.get("AWS_SECRET_ACCESS_KEY")
        self.compress = compress
        self.enabled = enabled and self.bucket_name is not None
        
        # Auto-enable if S3 bucket is configured
        if self.bucket_name and not enabled:
            self.enabled = True
    
    @classmethod
    def from_env(cls) -> "S3Config":
        """Create S3Config from environment variables."""
        bucket_name = os.environ.get("BIOMLBENCH_S3_BUCKET")
        default_prefix = os.environ.get("BIOMLBENCH_S3_PREFIX")
        
        return cls(
            bucket_name=bucket_name,
            prefix=default_prefix or "",  # Use empty string if not specified
            region=os.environ.get("AWS_DEFAULT_REGION", "us-east-1"),
            compress=os.environ.get("BIOMLBENCH_S3_COMPRESS", "true").lower() in ("true", "1"),
            enabled=os.environ.get("BIOMLBENCH_S3_ENABLED", "auto") != "false",
        )
    
    def is_valid(self) -> bool:
        """Check if S3 configuration is valid."""
        return (
            True
            and self.enabled
            and self.bucket_name is not None
            and len(self.bucket_name.strip()) > 0
        )


class S3Uploader:
    """Handles S3 uploads for biomlbench artifacts."""
    
    def __init__(self, config: S3Config):
        self.config = config
        self.client = None
        
        if config.is_valid():
            try:
                self.client = boto3.client(
                    's3',
                    region_name=config.region,
                    aws_access_key_id=config.aws_access_key_id,
                    aws_secret_access_key=config.aws_secret_access_key,
                )
                # Test connection
                self.client.head_bucket(Bucket=config.bucket_name)
                logger.info(f"S3 uploader initialized for bucket: {config.bucket_name}")
            except Exception as e:
                logger.error(f"Failed to initialize S3 client: {e}")
                self.client = None
                self.config.enabled = False
    
    def is_enabled(self) -> bool:
        """Check if S3 uploading is enabled and configured properly."""
        return self.config.is_valid() and self.client is not None
    
    def _compress_directory(self, directory: Path) -> Path:
        """Compress a directory to a tar.gz file."""
        with tempfile.NamedTemporaryFile(suffix='.tar.gz', delete=False) as tmp_file:
            tmp_path = Path(tmp_file.name)
        
        logger.info(f"Compressing {directory.name} to {tmp_path}")
        with tarfile.open(tmp_path, 'w:gz') as tar:
            tar.add(directory, arcname=directory.name)
        
        return tmp_path
    
    def _compress_file(self, file_path: Path) -> Path:
        """Compress a single file using gzip."""
        with tempfile.NamedTemporaryFile(suffix='.gz', delete=False) as tmp_file:
            tmp_path = Path(tmp_file.name)
        
        logger.info(f"Compressing {file_path.name} to {tmp_path}")
        with open(file_path, 'rb') as f_in:
            with gzip.open(tmp_path, 'wb') as f_out:
                f_out.writelines(f_in)
        
        return tmp_path
    
    def _upload_file(self, local_path: Path, s3_key: str, metadata: Optional[Dict[str, Any]] = None) -> bool:
        """Upload a single file to S3."""
        try:
            extra_args = {}
            if metadata:
                extra_args['Metadata'] = {k: str(v) for k, v in metadata.items()}
            
            file_size = local_path.stat().st_size
            logger.info(f"Uploading {local_path.name} ({file_size / (1024*1024):.1f} MB) to s3://{self.config.bucket_name}/{s3_key}")
            
            self.client.upload_file(
                str(local_path),
                self.config.bucket_name,
                s3_key,
                ExtraArgs=extra_args
            )
            logger.info(f"✅ Successfully uploaded {s3_key}")
            return True
            
        except (BotoCoreError, ClientError) as e:
            logger.error(f"❌ Failed to upload {s3_key}: {e}")
            return False
        except Exception as e:
            logger.error(f"❌ Unexpected error uploading {s3_key}: {e}")
            return False
    
    def upload_run_group(self, run_group_dir: Path, run_group_id: str, agent_id: str = None, task_ids: list = None) -> bool:
        """
        Upload a complete run group directory to S3.
        
        Args:
            run_group_dir: Path to the run group directory
            run_group_id: Identifier for the run group
            agent_id: Agent identifier for organizing uploads
            task_ids: List of task IDs in this run group
            
        Returns:
            True if upload succeeded, False otherwise
        """
        if not self.is_enabled():
            logger.debug("S3 upload disabled, skipping run group upload")
            return False
        
        if not run_group_dir.exists():
            logger.error(f"Run group directory does not exist: {run_group_dir}")
            return False
        
        try:
            # Determine S3 path structure based on agent and task information
            if agent_id and task_ids:
                if len(task_ids) == 1:
                    # Single task: organize by agent/task_id/
                    task_id_safe = task_ids[0].replace('/', '-').replace('_', '-')
                    if self.config.prefix:
                        s3_key_base = f"{self.config.prefix}/runs/{agent_id}/{task_id_safe}/{run_group_id}"
                    else:
                        s3_key_base = f"runs/{agent_id}/{task_id_safe}/{run_group_id}"
                else:
                    # Multiple tasks: organize by agent/multi-task/
                    if self.config.prefix:
                        s3_key_base = f"{self.config.prefix}/runs/{agent_id}/multi-task/{run_group_id}"
                    else:
                        s3_key_base = f"runs/{agent_id}/multi-task/{run_group_id}"
            else:
                # Fallback to original structure if agent/task info not available
                if self.config.prefix:
                    s3_key_base = f"{self.config.prefix}/runs/{run_group_id}"
                else:
                    s3_key_base = f"runs/{run_group_id}"
            
            if self.config.compress:
                # Compress the entire directory
                compressed_path = self._compress_directory(run_group_dir)
                s3_key = f"{s3_key_base}.tar.gz"
                
                metadata = {
                    'biomlbench_type': 'run_group',
                    'run_group_id': run_group_id,
                    'agent_id': agent_id or 'unknown',
                    'task_ids': ','.join(task_ids) if task_ids else 'unknown',
                    'compressed': 'true',
                    'upload_timestamp': str(int(time.time()))
                }
                
                success = self._upload_file(compressed_path, s3_key, metadata)
                
                # Clean up temporary file
                compressed_path.unlink()
                return success
            else:
                # Upload individual files (not recommended for large run groups)
                logger.warning("Uploading uncompressed run groups is not recommended for large benchmarks")
                # Implementation for individual file upload would go here
                return False
                
        except Exception as e:
            logger.error(f"Failed to upload run group {run_group_id}: {e}")
            return False
    
    def upload_grading_results(self, summary_report_path: Path, individual_reports_dir: Path, timestamp: str, agent_id: str = None, task_id: str = None) -> bool:
        """
        Upload grading results to S3.
        
        Args:
            summary_report_path: Path to the summary grading report JSON
            individual_reports_dir: Path to the directory with individual task reports
            timestamp: Timestamp identifier for the grading run
            agent_id: Agent identifier for organizing uploads
            task_id: Task identifier for organizing uploads
            
        Returns:
            True if upload succeeded, False otherwise
        """
        if not self.is_enabled():
            logger.debug("S3 upload disabled, skipping grading results upload")
            return False
        
        try:
            # Determine S3 path structure based on agent and task information
            if agent_id and task_id:
                task_id_safe = task_id.replace('/', '-').replace('_', '-')
                if self.config.prefix:
                    s3_key_base = f"{self.config.prefix}/grades/{agent_id}/{task_id_safe}/{timestamp}"
                else:
                    s3_key_base = f"grades/{agent_id}/{task_id_safe}/{timestamp}"
            else:
                # Fallback to original structure if agent/task info not available
                if self.config.prefix:
                    s3_key_base = f"{self.config.prefix}/grades/{timestamp}"
                else:
                    s3_key_base = f"grades/{timestamp}"
            success = True
            
            # Upload summary report
            if summary_report_path.exists():
                if self.config.compress:
                    compressed_summary = self._compress_file(summary_report_path)
                    s3_key = f"{s3_key_base}_grading_report.json.gz"
                    metadata = {
                        'biomlbench_type': 'grading_summary',
                        'timestamp': timestamp,
                        'agent_id': agent_id or 'unknown',
                        'task_id': task_id or 'unknown',
                        'compressed': 'true'
                    }
                    success &= self._upload_file(compressed_summary, s3_key, metadata)
                    compressed_summary.unlink()
                else:
                    s3_key = f"{s3_key_base}_grading_report.json"
                    metadata = {
                        'biomlbench_type': 'grading_summary',
                        'timestamp': timestamp,
                        'agent_id': agent_id or 'unknown',
                        'task_id': task_id or 'unknown',
                        'compressed': 'false'
                    }
                    success &= self._upload_file(summary_report_path, s3_key, metadata)
            
            # Upload individual reports directory
            if individual_reports_dir.exists():
                if self.config.compress:
                    compressed_individual = self._compress_directory(individual_reports_dir)
                    s3_key = f"{s3_key_base}_individual_reports.tar.gz"
                    metadata = {
                        'biomlbench_type': 'individual_reports',
                        'timestamp': timestamp,
                        'agent_id': agent_id or 'unknown',
                        'task_id': task_id or 'unknown',
                        'compressed': 'true'
                    }
                    success &= self._upload_file(compressed_individual, s3_key, metadata)
                    compressed_individual.unlink()
                else:
                    # Upload individual files in the directory
                    for report_file in individual_reports_dir.glob("*.json"):
                        s3_key = f"{s3_key_base}_individual_reports/{report_file.name}"
                        metadata = {
                            'biomlbench_type': 'individual_report',
                            'timestamp': timestamp,
                            'agent_id': agent_id or 'unknown',
                            'task_id': report_file.stem.replace('_', '/'),
                            'compressed': 'false'
                        }
                        success &= self._upload_file(report_file, s3_key, metadata)
            
            return success
            
        except Exception as e:
            logger.error(f"Failed to upload grading results for {timestamp}: {e}")
            return False
    
    def upload_incremental_run(self, run_dir: Path, run_group_id: str, task_id: str, run_id: str) -> bool:
        """
        Upload a single task run directory incrementally during execution.
        
        This is useful for long-running benchmarks to avoid losing progress.
        
        Args:
            run_dir: Path to individual run directory
            run_group_id: Run group identifier
            task_id: Task identifier
            run_id: Individual run identifier
            
        Returns:
            True if upload succeeded, False otherwise
        """
        if not self.is_enabled():
            return False
        
        if not run_dir.exists():
            logger.error(f"Run directory does not exist: {run_dir}")
            return False
        
        try:
            # Sanitize task_id for S3 key
            safe_task_id = task_id.replace('/', '_').replace('\\', '_')
            if self.config.prefix:
                s3_key_base = f"{self.config.prefix}/runs/{run_group_id}/incremental/{safe_task_id}_{run_id}"
            else:
                s3_key_base = f"runs/{run_group_id}/incremental/{safe_task_id}_{run_id}"
            
            if self.config.compress:
                compressed_path = self._compress_directory(run_dir)
                s3_key = f"{s3_key_base}.tar.gz"
                
                metadata = {
                    'biomlbench_type': 'incremental_run',
                    'run_group_id': run_group_id,
                    'task_id': task_id,
                    'run_id': run_id,
                    'compressed': 'true',
                    'upload_timestamp': str(int(time.time()))
                }
                
                success = self._upload_file(compressed_path, s3_key, metadata)
                compressed_path.unlink()
                return success
            else:
                # For uncompressed, we could upload key files individually
                # But this is not recommended for frequent incremental uploads
                logger.warning("Incremental uploads work best with compression enabled")
                return False
                
        except Exception as e:
            logger.error(f"Failed to upload incremental run {run_id}: {e}")
            return False

    def upload_failed_run_group(self, run_group_dir: Path, run_group_id: str) -> bool:
        """Upload failed run group to failed_runs folder."""
        if not self.is_enabled():
            return False
        
        if self.config.prefix:
            s3_key_base = f"{self.config.prefix}/failed_runs/{run_group_id}"
        else:
            s3_key_base = f"failed_runs/{run_group_id}"
        
        if self.config.compress:
            compressed_path = self._compress_directory(run_group_dir)
            s3_key = f"{s3_key_base}.tar.gz"
            metadata = {'biomlbench_type': 'failed_run_group', 'run_group_id': run_group_id}
            success = self._upload_file(compressed_path, s3_key, metadata)
            compressed_path.unlink()
            return success
        return False

    def upload_failed_grading_results(self, summary_report_path: Path, individual_reports_dir: Path, timestamp: str, agent_id: str = None, task_id: str = None) -> bool:
        """Upload failed grading results to failed_grades folder."""
        if not self.is_enabled():
            return False
        
        # Determine S3 path structure based on agent and task information
        if agent_id and task_id:
            task_id_safe = task_id.replace('/', '-').replace('_', '-')
            if self.config.prefix:
                s3_key_base = f"{self.config.prefix}/failed_grades/{agent_id}/{task_id_safe}/{timestamp}"
            else:
                s3_key_base = f"failed_grades/{agent_id}/{task_id_safe}/{timestamp}"
        else:
            # Fallback to original structure if agent/task info not available
            if self.config.prefix:
                s3_key_base = f"{self.config.prefix}/failed_grades/{timestamp}"
            else:
                s3_key_base = f"failed_grades/{timestamp}"
        
        success = True
        if summary_report_path and summary_report_path.exists() and self.config.compress:
            compressed_summary = self._compress_file(summary_report_path)
            s3_key = f"{s3_key_base}_grading_report.json.gz"
            success &= self._upload_file(compressed_summary, s3_key, {'biomlbench_type': 'failed_grading_summary'})
            compressed_summary.unlink()
        
        if individual_reports_dir and individual_reports_dir.exists() and self.config.compress:
            compressed_individual = self._compress_directory(individual_reports_dir)
            s3_key = f"{s3_key_base}_individual_reports.tar.gz"
            success &= self._upload_file(compressed_individual, s3_key, {'biomlbench_type': 'failed_individual_reports'})
            compressed_individual.unlink()
        
        return success


# Global S3 uploader instance
_s3_uploader: Optional[S3Uploader] = None


def get_s3_uploader() -> Optional[S3Uploader]:
    """Get or create the global S3 uploader instance."""
    global _s3_uploader
    if _s3_uploader is None:
        config = S3Config.from_env()
        _s3_uploader = S3Uploader(config)
    return _s3_uploader


def upload_run_group_artifacts(run_group_dir: Path, run_group_id: str, agent_id: str = None, task_ids: list = None) -> bool:
    """
    Convenience function to upload run group artifacts.
    
    Args:
        run_group_dir: Path to the run group directory  
        run_group_id: Run group identifier
        agent_id: Agent identifier for organizing uploads
        task_ids: List of task IDs in this run group
        
    Returns:
        True if upload succeeded or S3 is disabled, False on error
    """
    uploader = get_s3_uploader()
    if uploader and uploader.is_enabled():
        return uploader.upload_run_group(run_group_dir, run_group_id, agent_id, task_ids)
    return True  # Return True if S3 is disabled (not an error condition)


def upload_grading_artifacts(summary_report_path: Path, individual_reports_dir: Path, timestamp: str, agent_id: str = None, task_id: str = None) -> bool:
    """
    Convenience function to upload grading artifacts.
    
    Args:
        summary_report_path: Path to the summary report
        individual_reports_dir: Path to individual reports directory
        timestamp: Grading timestamp
        agent_id: Agent identifier for organizing uploads
        task_id: Task identifier for organizing uploads
        
    Returns:
        True if upload succeeded or S3 is disabled, False on error
    """
    uploader = get_s3_uploader()
    if uploader and uploader.is_enabled():
        return uploader.upload_grading_results(summary_report_path, individual_reports_dir, timestamp, agent_id, task_id)
    return True  # Return True if S3 is disabled (not an error condition)


def upload_incremental_artifacts(run_dir: Path, run_group_id: str, task_id: str, run_id: str) -> bool:
    """
    Convenience function to upload individual run artifacts incrementally.
    
    Args:
        run_dir: Path to the run directory
        run_group_id: Run group identifier  
        task_id: Task identifier
        run_id: Run identifier
        
    Returns:
        True if upload succeeded or S3 is disabled, False on error
    """
    uploader = get_s3_uploader()
    if uploader and uploader.is_enabled():
        return uploader.upload_incremental_run(run_dir, run_group_id, task_id, run_id)
    return True  # Return True if S3 is disabled (not an error condition)


def upload_failed_run_group_artifacts(run_group_dir: Path, run_group_id: str) -> bool:
    """Upload failed run group artifacts to failed_runs folder."""
    uploader = get_s3_uploader()
    if uploader and uploader.is_enabled():
        return uploader.upload_failed_run_group(run_group_dir, run_group_id)
    return True


def upload_failed_grading_artifacts(summary_report_path: Path, individual_reports_dir: Path, timestamp: str, agent_id: str = None, task_id: str = None) -> bool:
    """Upload failed grading artifacts to failed_grades folder."""
    uploader = get_s3_uploader()
    if uploader and uploader.is_enabled():
        return uploader.upload_failed_grading_results(summary_report_path, individual_reports_dir, timestamp, agent_id, task_id)
    return True 