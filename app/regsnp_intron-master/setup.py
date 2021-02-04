try:
    from setuptools import setup,find_packages
except ImportError:
    from distutils.core import setup
    

setup(
    name='regsnp_intron',
    version='0.1.6',
    packages=find_packages(),
    include_package_data=True,
    scripts=['bin/regsnp_intron'],
    classifiers=['Development Status :: 3 - Alpha',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 3.8',
                 'Topic :: Scientific/Engineering :: Bio-Informatics'],
    url='https://github.com/linhai86/regsnp_intron',
    license='MIT',
    author='linhai',
    author_email='linhai@iupui.edu',
    description='Predict disease-causing probability of human intronic SNVs.'
)
#print(packages)
