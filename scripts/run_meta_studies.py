#!/usr/bin/env python3
"""
CLI wrapper for ``meta_studies.py``, a module in analysis_tools (lives in Privitera_335 directory)

Usage (same as meta_studies.py):
    python run_meta_studies.py -j config/VR_study.json
"""

import sys
from pathlib import Path
from analysis_tools.meta_studies import main

if __name__ == "__main__":
    main()
