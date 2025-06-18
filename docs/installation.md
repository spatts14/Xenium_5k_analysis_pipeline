# Installation Guide

To install the most basic version of this project, create a python virtual environment
and run:

```bash
pip install git+https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics.git
```

However, there are some analysis modules that use the tool `MuSpAn`, which requires a
license to install. If you want to use those modules, you will need to request a
license. Do this by following the link in the [MuSpAn installation instructions]. Once
you have your username and password, run the following command and enter them when
prompted.

```bash
pip install 'recode_st[muspan] @ git+https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics.git'
```

[MuSpAn installation instructions]: https://docs.muspan.co.uk/latest/Installation.html

## Development installation

To install the development version of the software, clone this repo and run

```bash
pip install -e ".[dev,muspan]"
```

Pinned versions of the dependencies are saved in two text files:

- `requirements-muspan.txt`
- `requirements-muspan.txt`

These are managed by `pip-tools`, which is a dependency. They are built using the
`pip-compile` command at the top of each file. You can also clean your current
environment by ensuring it is in-sync with the following commands:

```bash
pip-sync requirements*.txt
pip install -e .
```
