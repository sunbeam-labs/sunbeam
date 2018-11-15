from distutils.core import setup

setup(
    name="sunbeam",
    use_scm_version = True,
    setup_requires=['setuptools_scm'],
    packages=["sunbeamlib"],
    include_package_data=True,
    package_data={"sunbeamlib": ["sunbeamlib/data/*.yml"]},
    entry_points={'console_scripts': [
        'sunbeam = sunbeamlib.scripts.command:main',
        'sunbeam_init = sunbeamlib.scripts.init:main',
        'sunbeam_get = sunbeamlib.scripts.get:main',
        'sunbeam_mod_config = sunbeamlib.scripts.config:main'
    ]},
    classifiers=[
        'Programming Language :: Python :: 3.4'
    ]
)
