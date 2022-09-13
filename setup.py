
import os,sys
from setuptools import setup


# build long description
base_dir = os.path.dirname(os.path.abspath(__file__))
long_description = '\n\n'.join([open(os.path.join(base_dir,'README.md'),'r').read()])

# requires
requires = [
    'squidpy ==1.2.0',
    'gseapy ==0.10.4',
    'scrublet ==0.2.3',
    'yattag',
    'anytree',
    'mygene ==3.2.2'
]

setup(
      name = 'sctriangulate',
      version = '0.12.0',
      description= 'A Python package to mix-and-match conflicting clustering results in single cell analysis, and generate reconciled clustering solutions.',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Guangyuan(Frank) Li',
      author_email='li2g2@mail.uc.edu',
      maintainer='Guangyuan(Frank) Li',
      maintainer_email='li2g2@mail.uc.edu',
      url='https://github.com/frankligy/scTriangulate',
      project_urls={
          'Ducumentation':'https://sctriangulate.readthedocs.io',
      },
      packages=['sctriangulate'],
      package_data = {'sctriangulate':['artifact_genes.txt','viewer/*','433colorhex.txt']},
      install_requires=requires,
      python_requires='>=3.7',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Programming Language :: Python :: 3',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
      ]

)
























