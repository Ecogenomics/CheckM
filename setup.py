from distutils.core import setup

setup(
    name='MetaChecka2000',
    version='0.0.1',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['metachecka2000', 'metachecka2000.test'],
    scripts=['bin/metachecka2000.py'],
    url='http://pypi.python.org/pypi/MetaChecka2000/',
    license='LICENSE.txt',
    description='MetaChecka2000',
    long_description=open('README.txt').read(),
    install_requires=[],
)
