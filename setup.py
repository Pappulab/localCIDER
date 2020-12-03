from setuptools import setup
import sys

setup(
    name='localcider',
    version='0.1.18',
    author='Alex Holehouse',
    author_email='alex.holehouse@wustl.edu',
    packages=['localcider', 'localcider.tests', 'localcider.backend', 'localcider.backend.data'],
    scripts=[],
    url='http://pappulab.github.io/localCIDER/',
    license='LICENSE.txt',
    description='Tools for calculating sequence properties of disordered proteins [from the Pappu Lab at Washington University in St. Louis]',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy","matplotlib","scipy"],
    test_suite='localcider.tests.suite')
#    **extras)
