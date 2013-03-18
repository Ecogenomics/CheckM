from distutils.core import setup

setup(
    name='checkM',
    version='0.3.1',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['metachecka2000', 'metachecka2000.test'],
    scripts=['bin/metachecka2000', 'bin/checkM', 'bin/uniqueMarkers.py'],
    url='http://pypi.python.org/pypi/MetaChecka2000/',
    license='GPL3',
    description='Beta check your meta before you wreck your meta',
    long_description=open('README.md').read(),
    install_requires=[
                      'simplehmmer >= 0.2.2',
                      'biopython >= 1.58'
                      ],
)
