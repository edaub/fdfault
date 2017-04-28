from setuptools import setup

setup(name='fdfault',
      version='1.0',
      description='Tools made to go with fdfault rupture code',
      url='https://github.com/egdaub/fdfault',
      author='Eric Daub',
      author_email='egdaub@memphis.edu',
      packages=['fdfault', 'fdfault.analysis'],
      install_requires=['numpy'])
