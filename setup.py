from setuptools import setup, find_packages

setup(
    name="blackbox",
    version="1.0.0",
    packages=find_packages(include=["blackbox*"]),
    python_requires=">=3.11",
    install_requires=["numpy>=1.26.0", 
                      "scipy>=1.11.1",
                      "pyDOE3>=1.1.0"
                    ],
)
