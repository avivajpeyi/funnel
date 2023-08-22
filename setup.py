# -*- coding: utf-8 -*-
from setuptools import setup

package_dir = {"": "src"}

packages = ["funnel"]

package_data = {"": ["*"]}

install_requires = ["loguru", "matplotlib", "numpy", "pytest", "tqdm", "numba"]

entry_points = {"console_scripts": ["funnel = funnel.__main__:main"]}

setup_kwargs = {
    "name": "funnel",
    "version": "0.0.1",
    "description": "Funnel",
    "long_description": "Funnel (FoUrier iNtegratioN Evidence calcuLation)\n=================================================\n\n|PyPI| |Python Version|\n\n|Tests| |Codecov|\n\n|pre-commit| |Black|\n\n.. |PyPI| image:: https://img.shields.io/pypi/v/funnel.svg\n   :target: https://pypi.org/project/funnel/\n   :alt: PyPI\n.. |Python Version| image:: https://img.shields.io/pypi/pyversions/funnel\n   :target: https://pypi.org/project/funnel\n   :alt: Python Version\n.. |Tests| image:: https://github.com/avivajpeyi/funnel/workflows/Tests/badge.svg\n   :target: https://github.com/avivajpeyi/funnel/actions?workflow=Tests\n   :alt: Tests\n.. |Codecov| image:: https://codecov.io/gh/avivajpeyi/funnel/branch/master/graph/badge.svg\n   :target: https://codecov.io/gh/avivajpeyi/funnel\n   :alt: Codecov\n.. |pre-commit| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white\n   :target: https://github.com/pre-commit/pre-commit\n   :alt: pre-commit\n.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg\n   :target: https://github.com/psf/black\n   :alt: Black\n\nContributing\n------------\n\nHow to set up your development environment\n------------------------------------------\n\nYou need Python 3.6+ and the following tools:\n\n- Poetry_\n- Nox_\n\nInstall the package with development requirements:\n\n.. code:: console\n\n   $ poetry install\n\nYou can now run an interactive Python session,\nor the command-line interface:\n\n.. code:: console\n\n   $ poetry run python\n   $ poetry run funnel\n\n.. _Poetry: https://python-poetry.org/\n.. _Nox: https://nox.thea.codes/\n\n\nHow to test the project\n-----------------------\n\nRun the full test suite:\n\n.. code:: console\n\n   $ nox\n\nList the available Nox sessions:\n\n.. code:: console\n\n   $ nox --list-sessions\n\nYou can also run a specific Nox session.\nFor example, invoke the unit test suite like this:\n\n.. code:: console\n\n   $ nox --session=tests\n\nUnit tests are located in the ``tests`` directory,\nand are written using the pytest_ testing framework.\n\n.. _pytest: https://pytest.readthedocs.io/\n\nHow to lint the project\n-----------------------\n\nYou can ensure that your changes adhere to the code style by reformatting with Black_:\n\n.. code:: console\n\n   $ nox --session=black\n\nIt is recommended to open an issue before starting work on anything.\nThis will allow a chance to talk it over with the owners and validate your approach.\n\n.. _pull request: https://github.com/avivajpeyi/funnel/pulls\n.. _Black: https://black.readthedocs.io/\n.. github-only\n",
    "author": "Avi Vajpeyi",
    "author_email": "avi.vajpeyi@gmail.com",
    "maintainer": "None",
    "maintainer_email": "None",
    "url": "https://github.com/avivajpeyi/funnel",
    "package_dir": package_dir,
    "packages": packages,
    "package_data": package_data,
    "install_requires": install_requires,
    "entry_points": entry_points,
    "python_requires": ">=3.9.1,<4.0.0",
}


setup(**setup_kwargs)
