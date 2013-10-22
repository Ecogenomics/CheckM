from distutils.core import setup

setup(
    name='checkm-ACE',
    version='1.0.0',
    author='Donovan Parks, Connor Skennerton, Michael Imelfort',
    author_email='donovan.parks@gmail.com',
    packages=['checkm', 'checkm.test'],
    scripts=['bin/checkm', 'bin/uniqueMarkers.py'],
    url='http://pypi.python.org/pypi/checkm/',
    license='GPL3',
    description='Estimate completeness and contamination of putative genomes.',
    long_description=open('README.md').read(),
    install_requires=[
                      'simplehmmer >= 0.2.4',
                      'biopython >= 1.58'
                      ],
)
