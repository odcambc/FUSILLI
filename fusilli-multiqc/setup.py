"""
Setup script for FUSILLI MultiQC custom modules.

This package registers FUSILLI-specific MultiQC modules as plugins
so they can be discovered and loaded by MultiQC.
"""

from setuptools import setup

setup(
    name="fusilli-multiqc",
    version="1.0.0",
    description="Custom MultiQC modules for FUSILLI fusion detection pipeline",
    author="Chris Macdonald",
    license="MIT",
    python_requires=">=3.8",
    packages=["fusilli_multiqc", "fusilli_multiqc.modules"],
    install_requires=[
        "multiqc>=1.0",
        "pandas>=1.0",
        "numpy>=1.0",
    ],
    entry_points={
        "multiqc.modules.v1": [
            "fusilli_detection = fusilli_multiqc.modules.detection:MultiqcModule",
            "fusilli_diversity = fusilli_multiqc.modules.diversity:MultiqcModule",
            "fusilli_preprocessing = fusilli_multiqc.modules.preprocessing:MultiqcModule",
            "fusilli_partners = fusilli_multiqc.modules.partners:MultiqcModule",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
