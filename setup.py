#!/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as input:
    long_description = input.read()

setup(
    name="CLCR",
    version="1.0.0",
    python_requires='>=3.7.0',
    description="Tool fo low coverage region correction",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Johannes Zieres",
    author_email="johannes.zieres@gmail.com",
    url="https://github.com/Johannes-Zi/CLCR",
    packages=find_packages(),
    package_data={'': ['*']},
    install_requires=[
        'pandas',
    ],
    entry_points={
        'console_scripts': ["clcr.query_creation = CLCR_func.query_creation:main",
                            "clcr.cluster_run = CLCR_func.slurmarray_creation:main",
                            "clcr.assembly_healing = CLCR_func.assembly_healing:main"],
    },
    license="GPL-3.0",
    classifiers=[
        "Operating System :: POSIX :: Linux",
        "Environment :: Console",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: End Users/Desktop",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
)