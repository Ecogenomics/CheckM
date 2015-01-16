from distutils.core import setup

setup(
    name='checkm-genome',
    version='0.9.7',
    author='Donovan Parks, Michael Imelfort, Connor Skennerton',
    author_email='donovan.parks@gmail.com',
    packages=['checkm', 'checkm.plot', 'checkm.test', 'checkm.util'],
    scripts=['bin/checkm'],
    url='http://pypi.python.org/pypi/checkm/',
    license='GPL3',
    description='Assess the quality of putative genome bins.',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy >= 1.8.0",
        "scipy >= 0.9.0",
        "matplotlib >= 1.3.1",
        "pysam >= 0.7.4",
        "dendropy >= 3.12.0",
        "ScreamingBackpack >= 0.2.333"],
    package_data={'checkm' : ['DATA_CONFIG']}
)
