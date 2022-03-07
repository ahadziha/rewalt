"""
Setup rewal package
"""

from setuptools import find_packages, setup


setup(name='rewal',
      package_dir={'rewal': 'rewal'},
      packages=find_packages(),
      description='Rewriting and Algebraic Topology',
      version='0.0.1',
      tests_suite='tests',
      )
