from setuptools import setup, find_packages

setup(
    name="hsp",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'requests>=2.25.0',
        'numpy',
        'flask',
        'scipy',
        'websockets',
        'rasterio',
        'pyproj',
    ],
)