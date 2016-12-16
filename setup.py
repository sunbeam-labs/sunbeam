from distutils.core import setup

setup(
    name="sunbeam",
    version="0.1",
    packages=["sunbeamlib"],
    package_data={"sunbeamlib": ["sunbeamlib/data/config_template.yml"]},
    entry_points={'console_scripts': ['sunbeam_init = sunbeamlib.initialize:main']}
)
