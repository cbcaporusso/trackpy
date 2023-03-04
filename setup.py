from setuptools import setup, find_packages

packages = find_packages(
    where='trackpy',
    include=['trackpy.*'],
)

# print(packages)

setup(
    name='trackpy',
    version='0.0.1',
    packages=packages,
    requires=["numpy", "matplotlib", "scipy"]
)
