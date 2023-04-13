try:
    from setuptools import setup, find_packages
except ImportError:
    exit("Please install setuptools.")

import os
import urllib

try:
    from urllib.request import urlretrieve

    classifiers = [
        "Programming Language :: Python",
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
except ImportError:
    from urllib.request import urlretrieve

    classifiers = [
        "Programming Language :: Python",
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]

VERSION = "1.0.1"
AUTHOR = "Mahdi Baghbanzadeh"
AUTHOR_EMAIL = "mbagh@gwu.edu"
MAINTAINER = "Mahdi Baghbanzadeh"
MAINTAINER_EMAIL = "mbagh@gwu.edu"

# try to download the bitbucket counter file to count downloads
# this has been added since PyPI has turned off the download stats
# this will be removed when PyPI Warehouse is production as it
# will have download stats
COUNTER_URL = "https://github.com/omicsEye/supp1/blob/master/README.md"
counter_file = "README.md"
if not os.path.isfile(counter_file):
    print("Downloading counter file to track supp1 downloads" +
          " since the global PyPI download stats are currently turned off.")
    try:
        pass  # file, headers = urlretrieve(COUNTER_URL,counter_file)
    except EnvironmentError:
        print("Unable to download counter")

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name="supp1",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    version=VERSION,
    license="MIT",
    description="supp1: creating a summary figure of PUBMED publications directions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/omicsEye/supp1",
    keywords=['data visualization', 'pubmed'],
    platforms=['Linux', 'MacOS', "Windows"],
    classifiers=classifiers,
    # long_description=open('readme.md').read(),

    install_requires=required,
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'supp1 = supp1.supp1:main'
        ]},
    test_suite='supp1.tests.supp1_test',
    zip_safe=False
)