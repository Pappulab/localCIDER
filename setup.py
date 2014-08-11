from distutils.core import setup

setup(
    name='localcider',
    version='0.1.0',
    author='Alex Holehouse',
    author_email='alex.holehouse@wustl.edu',
    packages=['localcider', 'localcider.tests', 'localcider.backend', 'localcider.backend.data'],
    scripts=[],
    url='http://pappulab.github.io/localCIDER/',
    license='LICENSE.txt',
    description='Tools for calculating sequence properties of disordered proteins [from the Pappu Lab in at Washington University in St. Louis]',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy","matplotlib"],
)
