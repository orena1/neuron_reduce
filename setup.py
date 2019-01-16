import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="neuron_reduce",
    version="0.0.3",
    author="Oren Amsalem, Guy Eyal, Noa Rogozinski, Idan Segev",
    author_email="oren.amsalem1@mail.huji.ac.il",
    description="Efficient analytical reduction of nonlinear detailed neuron models",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
) 
