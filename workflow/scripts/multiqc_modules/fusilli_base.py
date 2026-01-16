"""
Base module template for FUSILLI MultiQC modules.

This module provides utility functions for FUSILLI custom MultiQC modules
to use for data parsing and common operations.
"""

import logging
from pathlib import Path
from typing import Optional, Any
import pandas as pd


logger = logging.getLogger(__name__)


def parse_csv_file(file_path: str) -> Optional[pd.DataFrame]:
    """
    Parse a CSV file into a pandas DataFrame.
    
    Args:
        file_path: Path to CSV file
        
    Returns:
        pandas DataFrame or None if file doesn't exist or parsing fails
    """
    try:
        if not Path(file_path).exists():
            logger.warning(f"File not found: {file_path}")
            return None
        df = pd.read_csv(file_path)
        logger.info(f"Parsed {len(df)} rows from {file_path}")
        return df
    except Exception as e:
        logger.error(f"Error parsing CSV file {file_path}: {e}")
        return None


def parse_json_file(file_path: str) -> Optional[dict]:
    """
    Parse a JSON file into a dictionary.
    
    Args:
        file_path: Path to JSON file
        
    Returns:
        Dictionary or None if file doesn't exist or parsing fails
    """
    try:
        import json
        if not Path(file_path).exists():
            logger.warning(f"File not found: {file_path}")
            return None
        with open(file_path, 'r') as f:
            data = json.load(f)
        return data
    except Exception as e:
        logger.error(f"Error parsing JSON file {file_path}: {e}")
        return None
