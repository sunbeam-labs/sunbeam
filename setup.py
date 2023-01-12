from setuptools import setup, find_packages

setup(
    name="sunbeam",
    setup_requires=['setuptools'],
    install_requires=['pysam', 'semantic_version', 'pytest'],
    packages=find_packages(),
    include_package_data=True,
    package_data={"sunbeamlib": ["sunbeamlib/data/*.yml", "sunbeamlib/data/*.yaml"]},
    entry_points={'console_scripts': [
        'sunbeam = sunbeamlib.scripts.command:main',
        'sunbeam_init = sunbeamlib.scripts.init:main',
        'sunbeam_mod_config = sunbeamlib.scripts.config:main'
    ]},
    classifiers=[
        'Programming Language :: Python :: 3.4'
    ]
)
