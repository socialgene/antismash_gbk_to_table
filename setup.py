#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

requirements = []

test_requirements = []

setup(
    author="Chase M Clark",
    author_email="chasingmicrobes@gmail.com",
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="Python Boilerplate contains all the boilerplate you need to create a Python package.",
    entry_points={
        "console_scripts": [
            "antismash_gbk_to_table=antismash_gbk_to_table.cli:main",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords="antismash_gbk_to_table",
    name="antismash_gbk_to_table",
    packages=find_packages(
        include=["antismash_gbk_to_table", "antismash_gbk_to_table.*"]
    ),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/chasemc/antismash_gbk_to_table",
    version="0.1.0",
    zip_safe=False,
)
