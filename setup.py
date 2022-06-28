"""
Setup rewalt package
"""

from re import search, M
from setuptools import find_packages, setup

with open('rewalt/__init__.py', 'r') as file:
    MATCH = search(r"^__version__ = ['\"]([^'\"]*)['\"]", file.read(), M)
    if MATCH:
        VERSION = MATCH.group(1)
    else:
        raise RuntimeError('Unable to find version string.')

setup(name='rewalt',
      package_dir={'rewalt': 'rewalt'},
      packages=find_packages(),
      description='Rewriting Algebraic Topology',
      version=VERSION,
      tests_suite='tests'
      )
