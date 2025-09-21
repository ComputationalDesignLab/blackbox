from setuptools import setup, find_packages

setup(
    name='blackbox',
    version='0.9.0',
    packages=find_packages(include=['blackbox*']),
    install_requires=["numpy>=1.21.6", 
                      "scipy>=1.7.3", 
                      "mpi4py>=3.1.4", 
                      "pygeo>=1.12.2", 
                      "prefoil>=2.0.0", 
                      "smt>=2.0.0", 
                      "psutil>=5.9.1"],
)
