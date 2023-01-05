from setuptools import setup, find_packages

__author__ = "Zheren Wang"
__version__ = "1.0"
__maintainer__ = "Zheren Wang"
__email__ = "zherenwang@berkeley.edu"
__status__ = "Development"

if __name__ == '__main__':
    setup(name='synthesis_condition_optimizer',
          version=1.0,
          author="Zheren Wang",
          author_email="zherenwang@berkeley.edu",
          license="MIT License",
          packages=find_packages(),
          include_package_data=True,
          install_requires=[
              'numpy',
          ],
          zip_safe=False)