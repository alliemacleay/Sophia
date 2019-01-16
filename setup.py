from setuptools import setup, find_packages

__author__ = "Allison MacLeay"
__version__ = "1.0"


def load_requirements():
    req_list = []
    with open('requirements.txt', 'r') as fileobj:
        for line in fileobj:
            req_list.append(line.strip())
    return req_list

setup(
    name='Sophia',
    version=__version__,
    py_modules=['sophia'],
    packages=find_packages(),
    include_package_data=True,
    url='https://github.com/alliemacleay/Sophia',
    author='Allison MacLeay',
    author_email='allison.macleay@gmail.com',
    package_data={
        'data': "*.*"},
    license='MIT',
    install_requires=load_requirements(),
    entry_points='''
        [console_scripts]
        sophia=main
    ''',
)