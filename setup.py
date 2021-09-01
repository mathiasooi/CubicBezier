from setuptools import find_packages, setup

setup(
    name="CubicBezier",
    packages=find_packages(include=["CubicBezier"]),
    version="0.1.0",
    install_requires=[
        "matplotlib",
        "numpy",],
    description="Cubic Bezier curve library",
    author="Mathias Ooi",
    author_email="mathias.ho.ooi@gmail.com",
    license="MIT",
)