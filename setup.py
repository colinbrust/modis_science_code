from setuptools import setup, find_packages

setup(
    name='msc',
    version='0.0.1',
    description='Code to run and calibrate MOD16 and MOD17 algorithms.',
    author='Colin Brust',
    author_email='colin.brust@umontana.edu',
    license='MIT',
    url='https://github.com/colinbrust/modis_science_code',
    packages=find_packages(exclude=["test*", "data*", "examples*"])
)
