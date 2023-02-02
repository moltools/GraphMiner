import setuptools


setuptools.setup(
    name="GraphMiner",
    version="0.0.1",
    author="",
    author_email="",
    description="",
    url="",
    install_requires=["rdkit>=2022.03.5", "pandas"],
    packages=setuptools.find_packages(),
    python_requires=">=3.9",
    entry_points={"console_scripts": ["GraphMiner = GraphMiner.main:main"]}
)