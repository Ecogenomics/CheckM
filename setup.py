from distutils.core import setup

setup(
    name='checkm-ACE',
    version='1.0.0',
    author='Donovan Parks, Connor Skennerton, Michael Imelfort',
    author_email='donovan.parks@gmail.com',
    packages=['checkm', 'checkm.test'],
    scripts=['bin/checkm'],
    url='http://pypi.python.org/pypi/checkm/',
    license='GPL3',
    description='Estimate completeness and contamination of putative genomes.',
    long_description=open('README.md').read(),
    install_requires=[],
)
