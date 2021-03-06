==========================================================================
NDEx Cancer Hallmark networks from CPTAC at WikiPathways Enrichment Loader
==========================================================================


.. image:: https://img.shields.io/pypi/v/ndexcptacenrichloader.svg
        :target: https://pypi.python.org/pypi/ndexcptacenrichloader

.. image:: https://img.shields.io/travis/coleslaw481/ndexcptacenrichloader.svg
        :target: https://travis-ci.org/coleslaw481/ndexcptacenrichloader

.. image:: https://readthedocs.org/projects/ndexcptacenrichloader/badge/?version=latest
        :target: https://ndexcptacenrichloader.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Copies and converts Special collection of Cancer Hallmark networks
from CPTAC at WikiPathways networks already loaded in NDEx so they
can be used in Enrichment.


* Free software: BSD license
* Documentation: https://ndexcptacenrichloader.readthedocs.io.


Tools
-----

* **ndexloadcptacenrich.py** -- Loads CPTAC into NDEx_

Dependencies
------------

* `ndex2 <https://pypi.org/project/ndex2>`_
* `ndexutil <https://pypi.org/project/ndexutil>`_

Compatibility
-------------

* Python 3.3+

Installation
------------

.. code-block::

   git clone https://github.com/coleslaw481/ndexcptacenrichloader
   cd ndexcptacenrichloader
   make dist
   pip install dist/ndexloadcptacenrich*whl


Run **make** command with no arguments to see other build/deploy options including creation of Docker image 

.. code-block::

   make

Output:

.. code-block::

   clean                remove all build, test, coverage and Python artifacts
   clean-build          remove build artifacts
   clean-pyc            remove Python file artifacts
   clean-test           remove test and coverage artifacts
   lint                 check style with flake8
   test                 run tests quickly with the default Python
   test-all             run tests on every Python version with tox
   coverage             check code coverage quickly with the default Python
   docs                 generate Sphinx HTML documentation, including API docs
   servedocs            compile the docs watching for changes
   testrelease          package and upload a TEST release
   release              package and upload a release
   dist                 builds source and wheel package
   install              install the package to the active Python's site-packages
   dockerbuild          build docker image and store in local repository
   dockerpush           push image to dockerhub


Configuration
-------------

The **ndexloadcptacenrich.py** requires a configuration file in the following format be created.
The default path for this configuration is :code:`~/.ndexutils.conf` but can be overridden with
:code:`--conf` flag.

**Format of configuration file**

.. code-block::

    [<value in --profile (default ndexcptacenrichloader)>]

    user = <NDEx username>
    password = <NDEx password>
    server = <NDEx server(omit http) ie public.ndexbio.org>
    source_user = <NDEx username for source networks>
    source_password = <NDEx password for source networks>
    source_server = <NDEx server(omit http) ie public.ndexbio.org for source networks>
    source_networkset = <id of network set containing networks>


**Example configuration file**

.. code-block::

    [ndexcptacenrichloader_dev]

    user = joe123
    password = somepassword123
    server = dev.ndexbio.org
    source_user = bob456
    source_password = anotherpassword123
    source_server = public.ndexbio.org
    source_networkset = cf04cd21-b695-44b1-80c6-32952aaba2b9


Needed files
------------

**TODO:** Add description of needed files


Usage
-----

For information invoke :code:`ndexloadcptacenrich.py -h`

**Example usage**

**TODO:** Add information about example usage

.. code-block::

   ndexloadcptacenrich.py # TODO Add other needed arguments here


Via Docker
~~~~~~~~~~~~~~~~~~~~~~

**Example usage**

**TODO:** Add information about example usage


.. code-block::

   docker run -v `pwd`:`pwd` -w `pwd` coleslawndex/ndexcptacenrichloader:0.1.0 ndexloadcptacenrich.py --conf conf # TODO Add other needed arguments here


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
