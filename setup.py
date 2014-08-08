from distutils.core import setup

setup(
    name='localCIDER',
    version='0.1.0',
    author='Alex Holehouse',
    author_email='alex.holehouse@wustl.edu',
    packages=['localCIDER', 'localCIDER.test, localCIDER.kappa, localCIDER.wl, localCIDER.backend'],
    scripts=[''],
    url='http://pypi.python.org/pypi/TowelStuff/',
    license='LICENSE.txt',
    description='Tools for calculating sequence properties of disodered proteins',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy","matplotlib"],
    test_suite="tests",
)
