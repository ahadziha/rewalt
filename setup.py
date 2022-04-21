"""
Setup rewal package
"""

from re import search, M
from setuptools import find_packages, setup

with open('discopy/__init__.py', 'r') as file:
    MATCH = search(r"^__version__ = ['\"]([^'\"]*)['\"]", file.read(), M)
    if MATCH:
        VERSION = MATCH.group(1)
    else:
        raise RuntimeError('Unable to find version string.')

setup(name='rewal',
      package_dir={'rewal': 'rewal'},
      packages=find_packages(),
      description='Rewriting and Algebraic Topology',
      version=VERSION,
      tests_suite='tests'
      )
