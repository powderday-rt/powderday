from setuptools import setup, find_packages

try:
    import hyperion
except ImportError:
    raise Exception("Hyperion must be installed manually, see docs/installation.rst")

setup(
    name="powderday",
    version="0.1.0",
    packages=find_packages(),
    platforms="any",
    setup_requires=[
        'numpy==1.16; python_version=="2.7"',
        'numpy; python_version>="3.5"',
    ],
    install_requires=[
        'scipy==1.10.0; python_version=="2.7"',
        'scipy; python_version>="3.5"',
        'astropy==2.0; python_version=="2.7"',
        'astropy; python_version>="3.5"',
        'h5py>=2.9',
        'yt',
        'unyt',
        'fsps',
        'scikit-learn',
        'p_tqdm'
    ],
    scripts=["pd_front_end.py"],
    project_urls={
        'Source': 'https://github.com/dnarayanan/powderday.git',
    },
)
