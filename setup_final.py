#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name='pygsf',
    version='2.1',
    packages=find_packages(),

    install_requires=['numpy'],

    # metadata for upload to PyPI
    author="Mauro Alberti",
    author_email="alberti.m65@gmail.com",
    description="gsf, a library for structural geology",
    license="GPL-3",
    keywords="structural geology",
    url="https://github.com/mauroalberti/gsf",
    project_urls={
        "Bug Tracker": "https://github.com/mauroalberti/gsf/issues",
        "Documentation": "https://github.com/mauroalberti/gsf/blob/master/README.md",
        "Source Code": "https://github.com/mauroalberti/gsf/tree/master/gsf_py",
    }
)
