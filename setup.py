from os.path import join
from glob import glob

import setuptools


DATA_EXTS = {'*.h5'}


def getDataFiles():
    """Return all data files from ``pyIsoDep/data``"""

    files = []
    for ext in DATA_EXTS:
        files.extend(glob(join('pyIsoDep', 'data', ext)))
    return files


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyIsoDep",
    version="0.0.1",
    author="Dan Kotlyar",
    author_email="dan.kotlyar@me.gatech.edu",
    description="Python-based Isotopic Depletion tool for transmutation and "
                "decay analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DanKotlyar/PYTHON-ISOTOPIC-DEPLETION-PACKAGE",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    packages=['pyIsoDep', 'pyIsoDep.functions',
              'pyIsoDep.tests', 'pyIsoDep.data'],
    package_data={
        'pyIsoDep.data': ['data/{}'.format(ext) for ext in DATA_EXTS],
    },
    include_package_data=True,
    data_files=[('pyIsoDep/data', getDataFiles()), ],
)
