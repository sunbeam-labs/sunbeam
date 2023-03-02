from setuptools import setup, find_packages

setup(
    name="sunbeam",
    setup_requires=['setuptools'],
    install_requires=['more-itertools', 'semantic_version', 'pytest'],
    packages=find_packages(),
    include_package_data=True,
    package_data={"sunbeamlib": ["sunbeamlib/data/*.yml", "sunbeamlib/data/*.yaml"]},
    entry_points={'console_scripts': [
        'sunbeam = sunbeamlib.scripts.command:main'
    ]},
    classifiers=[
        'Programming Language :: Python :: 3.9'
    ]
)
