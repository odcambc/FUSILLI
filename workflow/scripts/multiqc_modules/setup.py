"""
Setup script for FUSILLI MultiQC custom modules.

This package registers FUSILLI-specific MultiQC modules as plugins
so they can be discovered and loaded by MultiQC.
"""

from setuptools import setup

# Note: This setup.py is in workflow/scripts/multiqc_modules/
# When installed with pip install -e ., the package name will be 'multiqc_modules'
# To avoid conflicts, we'll use a custom package_dir mapping
setup(
    name="fusilli-multiqc-modules",
    version="1.0.0",
    description="Custom MultiQC modules for FUSILLI fusion detection pipeline",
    py_modules=[
        "fusilli_detection",
        "fusilli_diversity", 
        "fusilli_preprocessing",
        "fusilli_partners",
        "fusilli_base",
    ],
    package_dir={"": "."},
    install_requires=[
        "multiqc>=1.0",
        "pandas>=1.0",
        "numpy>=1.0",
    ],
    entry_points={
        "multiqc.modules.v1": [
            "fusilli_detection = fusilli_detection:MultiqcModule",
            "fusilli_diversity = fusilli_diversity:MultiqcModule",
            "fusilli_preprocessing = fusilli_preprocessing:MultiqcModule",
            "fusilli_partners = fusilli_partners:MultiqcModule",
        ],
    },
    python_requires=">=3.8",
)
