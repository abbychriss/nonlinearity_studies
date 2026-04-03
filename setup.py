from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="nonlinearity_studies",
    version="0.1.0",
    author="Abby Chriss",
    description="Analysis tools for detector nonlinearity studies using FITS data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/nonlinearity_studies",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.26.0",
        "matplotlib>=3.8.0",
        "scipy>=1.11.0",
        "astropy>=5.3",
    ],
    entry_points={
        "console_scripts": [
            "run-nonlinearity-studies=nonlinearity_studies.run_nonlinearity_studies:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
