Branch ``geoframe-pomet``
========================

This branch adds the GEOframe-PoMet interpolation (Po River District, 1-km daily precipitation and temperature,
1991–2020) to the Krigings component: packages ``pokriging`` and ``krigingsCore`` in ``src/main/java``, the OMS
project and the Python pipeline in ``geoframe_pomet/`` (start with ``geoframe_pomet/README.md``), and the list of
changes to the upstream code in ``UPSTREAM_DIFF.md``.

Krigings
========

``@author`` Marialaura Bancheri

``@author`` Francesco Serafin (sidereus3), francesco.serafin.3@gmail.com

``@copyright`` GNU Public License v3 AboutHydrology (Riccardo Rigon)

Description
-----------

This is a Gradle project

.. figure:: doc/ReadMe/gradle.png
   :alt: Gradle logo

   Gradle logo

To build the project just run

::

    ~ $ gradle build

in the working directory. You will find the built ``.jar`` in
``build/libs``. But a built version of the lates release is already
available at the `GitHub release
section <https://github.com/geoframecomponents/Krigings/releases>`__

To build the ReadMe file from the markdown one in the doc/ReadMe folder

::

    pandoc doc/ReadMe/ReadMe.md -o ReadMe.rst --bibliography=doc/ReadMe/biblio.bib

Implementation
~~~~~~~~~~~~~~

Documentation
-------------

Developers' documentation
~~~~~~~~~~~~~~~~~~~~~~~~~

Linkers' documentation
~~~~~~~~~~~~~~~~~~~~~~

Users' documentation
~~~~~~~~~~~~~~~~~~~~

References
----------
