from setuptools import setup, find_packages


# Define your version here
version = "1.0.1"

def get_version():
    return version


setup(
    name='bempy',
    version=get_version(),
    packages=find_packages(),
    install_requires=[
        # Add any dependencies your package needs here
    ],
    # Metadata
    author='Dylan Green',
    description='Python Blade Element Momentum (BEM) solver.'
)
