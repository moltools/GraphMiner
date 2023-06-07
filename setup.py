import setuptools


setuptools.setup(
    name="GraphMiner",
    version="0.0.1",
    author="",
    author_email="",
    description="",
    url="",
    install_requires=["rdkit>=2022.03.5", "pandas>=1.1.3", "timer", "statsmodels>=0.13.5", "scipy>=1.10.1", "numpy>=1.24.2", "matplotlib", "IPython"],
    packages=setuptools.find_packages(),
    python_requires=">=3.8",
    entry_points={"console_scripts": ["GraphMiner = GraphMiner.main:main"]}
)