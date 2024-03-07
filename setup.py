from setuptools import setup

setup(
    name='microbeGWAS',
    version='1.0',
    packages=['microbeGWAS'],
    entry_points={
        'console_scripts': [
            'microbeGWAS=microbeGWAS.engine:main'
        ]
    },
)
