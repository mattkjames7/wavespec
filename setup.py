import setuptools
from setuptools.command.install import install
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="wavespec",
    version="0.0.3",
    author="Matthew Knight James",
    author_email="mattkjames7@gmail.com",
    description="Some spectral analysis tools for analyzing waves in data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mattkjames7/wavespec",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: POSIX",
    ],
    install_requires=[
		'numpy',
		'matplotlib',
		'scipy',
		'DateTimeTools',
	],
	include_package_data=True,
)



