import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

with open("requirements.txt", "r") as f:
    REQUIRED = f.readlines()

setuptools.setup(
    name="pymol-rosetta-utils",
    version="0.0.1",
    author="Brian D. Weitzner",
    author_email="brian.weitzner@gmail.com",
    description="PyMOL visualization shortcuts that use Rosetta.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires=">=3.7",
    install_requires=REQUIRED,
)
