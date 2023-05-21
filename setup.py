from setuptools import setup

setup(
  name='frastro',
  version='0.1.0',
  author='Camilo Jimenez',
  author_email='camiloejimenez@gmail.com',
  packages=['frastro', 'frastro.convolve','frastro.core','frastro.external','frastro.galaxymodel','frastro.gui','frastro.observation','frastro.telescope'],
  url='https://github.com/xcamilox/frastro',
  license='LICENSE.txt',
  description='Framework astrophysics',
  long_description=open('README.txt').read(),
  install_requires=[
     "astropy",
    "matplotlib",
    "pandas",
    "numpy"
  ],
)