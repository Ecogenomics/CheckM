from distutils.core import setup

setup(
    name='checkm-ACE',
    version='1.0.0',
    author='Donovan Parks, Michael Imelfort, Connor Skennerton',
    author_email='donovan.parks@gmail.com',
    packages=['checkm', 'checkm.plot', 'checkm.test'],
    scripts=['bin/checkm'],
    url='http://pypi.python.org/pypi/checkm/',
    license='GPL3',
    description='Estimate completeness and contamination of putative genomes.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.6.1",
        "scipy >= 0.10.1",
        "matplotlib >= 1.1.0",
        "pysam >= 0.7.4",
        "sqlite3 >= 2.6.0"],
)
