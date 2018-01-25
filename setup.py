from distutils.core import setup

setup(
    name="sunbeam",
    use_scm_version = True,
    setup_requires=['setuptools_scm'],
    packages=["sunbeamlib"],
    include_package_data=True,
    package_data={"sunbeamlib": ["sunbeamlib/data/*.yml"]},
    entry_points={'console_scripts': [
        'sunbeam_init = sunbeamlib.initialize:main',
        'sunbeam_mod_config = sunbeamlib.mod_config:main'
    ]},
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3.4'
    ]
)
