#!/usr/bin/env python
"""
Setup configuration for TETyper 2.0
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

setup(
    name="tetyper2",
    version="2.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Enhanced Transposable Element Typing Platform",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/TETyper2.0",
    project_urls={
        "Bug Tracker": "https://github.com/yourusername/TETyper2.0/issues",
        "Documentation": "https://tetyper2-docs.readthedocs.io",
        "Source Code": "https://github.com/yourusername/TETyper2.0",
    },
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.78",
        "numpy>=1.19.0",
        "pandas>=1.1.0",
        "pysam>=0.16.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.10.0",
            "black>=21.0",
            "flake8>=3.9.0",
            "isort>=5.9.0",
        ],
        "docs": [
            "sphinx>=3.0",
            "sphinx-rtd-theme>=0.5.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "tetyper2=cli:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
