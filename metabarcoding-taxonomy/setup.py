from setuptools import setup, find_packages

setup(
    name='metabarcoding_taxonomy',
    version='0.1.0',
    author='Seokwoo JO',
    author_email='newgenes1031@naver.com',
    description='ASV taxonomy filtering, stats and QIIME2-style plots',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/Newgenes1031/metabarcoding_taxonomy',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'pandas>=1.3',
        'numpy>=1.21',
        'matplotlib>=3.4',
    ],
    entry_points={
        'console_scripts': [
            'metabarcoding-taxonomy=metabarcoding_taxonomy.taxonomy:main',
        ],
    },
    python_requires='>=3.7',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
)