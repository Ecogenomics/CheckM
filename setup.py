from distutils.core import setup

setup(
    name='MetaChecka2000',
    version='0.1.0',
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['metachecka2000', 'metachecka2000.test'],
    scripts=['bin/metachecka2000'],
    url='http://pypi.python.org/pypi/MetaChecka2000/',
    license='GPL3',
    description='MetaChecka2000',
    long_description=open('README.txt').read(),
    install_requires=[
                      'simplehmmer >= 0.2.1'
                      ],
)
