from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='phase_fields_bo',
      version='1.0',
      description='Bayesian optimisation for accelerated exploration of phase fields',
      url='http://github.com/DrewNow/PhaseFieldsBO',
      author='Andrij Vasylenko',
      author_email='and.vasylenko@gmail.com',
      license='MIT',
      packages=['phase_fields_bo'],
      install_requires=['numpy', 'GPyOpt', 'pandas', 'pymatgen'],
      python_requires='>=3.7, <3.9',
      include_package_data=True,
      entry_points={"console_scripts": ["phase_fields_bo=phase_fields_bo.__main__:run"]},  
      zip_safe=False)
