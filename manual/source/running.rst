=============
Running BDSIM
=============

Basic Operation
===============

The minimum required information to run BDSIM is an input GMAD file. ::

  bdsim --file=lattice.gmad

The following section describes other 'executable' options that may specify
various options, such as whether to run in 'batch' mode or interactively with a
visualiser and where to write output.

.. _executable-options:

Executable Options
==================

BDSIM can be executed in a terminal with extra arguments to specify various inputs.
The angular brackets here are to indicate a user-specified value and should not
be used literally.  The following command may be used to display all options::

  bdsim --help

.. tabularcolumns:: |p{6cm}|p{9cm}|

+------------------------------+------------------------------------------------+
| -\\-file=<file>              | Specifies the input gmad file                  |
+------------------------------+------------------------------------------------+
| -\\-batch                    | Batch mode - no graphics                       |
+------------------------------+------------------------------------------------+
| -\\-circular                 | Assumes circular machine - turn control. See   |
|                              | note below.                                    |
+------------------------------+------------------------------------------------+
| -\\-colours                  | Print out all colours predefined in BDSIM and  |
|                              | exit                                           |
+------------------------------+------------------------------------------------+
| -\\-distrFile=<file>         | Which file to use for the bunch                |
|                              | distribution                                   |
+------------------------------+------------------------------------------------+
| -\\-exportGeometryTo=<file>  | Exports the geometry to a file                 |
|                              | extension, which determines the format         |
|                              | where possible extensions are ("gdml")         |
+------------------------------+------------------------------------------------+
| -\\-geant4Macro=<file>       | Optional Geant4 macro file to run after the    |
|                              | visualisation has started. Only works in       |
|                              | interactive visualisation.                     |
+------------------------------+------------------------------------------------+
| -\\-generatePrimariesOnly    | Generates primary particle coordinates only    |
|                              | then exits without simulating anything         |
+------------------------------+------------------------------------------------+
| -\\-help                     | Lists all executable commands                  |
+------------------------------+------------------------------------------------+
| -\\-materials                | Lists materials included in BDSIM by default   |
+------------------------------+------------------------------------------------+
| -\\-ngenerate=N              | The number of primary events to simulate       |
|                              | overrides the ngenerate option in the input    |
|                              | file.                                          |
+------------------------------+------------------------------------------------+
| -\\-output=<fmt>             | Outputs the format "rootevent" (default) or    |
|                              | "none"                                         |
+------------------------------+------------------------------------------------+
| -\\-outfile=<file>           | Outputs file name. Will be appended with _N    |
|                              | where N = 0, 1, 2, 3...                        |
+------------------------------+------------------------------------------------+
| -\\-recreate=<file>          | The rootevent output file to recreate events   |
|                              | from.                                          |
+------------------------------+------------------------------------------------+
| -\\-seed=<N>                 | Seed for the random number generator           |
+------------------------------+------------------------------------------------+
| -\\-seedStateFileName=<file> | File containing CLHEP::Random seed state       |
|                              | NB \- this overrides other seed values         |
+------------------------------+------------------------------------------------+
| -\\-startFromEvent=N         | Event offset to start from when recreating     |
|                              | events when using `-\\-recreate`               |
+------------------------------+------------------------------------------------+
| -\\-survey=<file>            | Prints survey info to <file>                   |
+------------------------------+------------------------------------------------+
| -\\-printFractionEvents=N    | How often to print out events as a fraction    |
|                              | of the total number of events to simulate      |
|                              | (default 0.1 i.e. 10%). -1 for all.            |
+------------------------------+------------------------------------------------+
| -\\-printFractionTurns=N     | How often to print out turns as a fraction     |
|                              | of the total number of turns to simulate       |
|                              | (default 0.2 i.e. 20%). -1 for all. Will       |
|                              | only print out on an event that will print     |
|                              | out as well.                                   |
+------------------------------+------------------------------------------------+
| -\\-verbose                  | Displays general parameters before run         |
+------------------------------+------------------------------------------------+
| -\\-verbose\_event           | Displays information for every event           |
+------------------------------+------------------------------------------------+
| -\\-verbose\_step            | Displays tracking information after each       |
|                              | step                                           |
+------------------------------+------------------------------------------------+
| -\\-verbose\_event\_num=<N>  | Displays tracking information for event 'N'    |
+------------------------------+------------------------------------------------+
| -\\-vis_debug                | Displays all volumes in the visualiser         |
+------------------------------+------------------------------------------------+
| -\\-vis_mac=<file>           | File with the visualisation macro script.      |
|                              | Default provided by BDSIM: openGL (OGLSQt)).   |
+------------------------------+------------------------------------------------+
| -\\-writeSeedState           | Writes an ASCII file seed state for each       |
|                              | event                                          |
+------------------------------+------------------------------------------------+

BDSIM can be run in one of two ways: `interactively`_ or `in batch mode`_, which
are described in the following sections.

When run interactively, a Geant4 visualiser is invoked that produces a window with an image
of the BDSIM model as well as a terminal prompt to control it. No events are simulated
without user input. BDSIM provides a basic visualisation "macro" for Geant4 using the Qt
visualiser to start and add the geometry to the scene. This is found from the
BDSIM installation directory or failing that the build directory. The user may provide
their own custom visualisuation macro with the executable command :code:`--vis_mac=mymac.mac`.

Alternatively, BDSIM can be run in batch mode, where no visualiser
is used and the specified number of primary events is simulated and feedback is printed
to the terminal. Batch mode is typically much faster than the interactive mode, but
the interactive mode is very useful for understanding the model and a typical event
in the simulation, i.e. where a particle hits.

.. note:: For more information on the `-\\-circular` option and using BDSIM with circular machines,
	  see :ref:`circular-machines`.

Interactively
=============

Features:

* Default option
* Interactive visualisation of the accelerator model
* Ability to view and rotate the accelerator model
* Ability to run individual events
* Typically slower than batch mode
* No events run without user input

To execute BDSIM in interactive mode, the user must simply not use the :code:`-\\-batch` command.
The user can also specify a macro file using the :code:`-\\-vis_mac` option above; otherwise, BDSIM
will look for "vis.mac" in the current working directory. If not found, BDSIM will use its own
default visualisation settings (typically: Qt visualiser, extra convenient buttons). Once
executed a window such as this will appear (depending on the visualiser you use).

.. figure:: figures/visualisation/qtvisualiser.png
   :width: 90%
   :align: center
   :figclass: align-center

   BDSIM running interactively with OpenGL Qt visualiser from Geant4.

Note: the visualiser is part of Geant4, so if the desired visualiser isn't available, you
must recompile Geant4 with the correct visualiser (and subsequently BDSIM afterwards). Geant4
also uses the CMake configuration system. The visualiser shown is the OpenGL Qt visualiser, which
we recommend for its ease of use and high level of interactivity.

More details can be found in :ref:`Visualisation`.

.. note:: BDSIM simulates one particle at a time from the primary distribution and all of the
	  associated secondaries. Each event is independent and different particles in the input
	  bunch cannot interact with each other or their secondaries. This is an underlying feature
	  of Geant4.

In Batch Mode
=============

Features:

* No interactive visualiser
* Faster
* Number of events runs, then program quits
* No user input
* Typical use for a "job" on a farm

To execute BDSIM in batch mode, simply use the :code:`-\\-batch` execution option.

Examples
========
::

   bdsim --file=atf2.gmad --outfile=run1 --batch --seed=123

This executes BDSIM for the ATF2 example with ROOT output to a file name "run1" in batch
mode with a seed value of 123. The simulation runs the number of events specified by the
:code:`ngenerate` options parameter in the input gmad file. ::

      bdsim --file=sm.gmad --outfile=run2

This executes BDSIM for the simpleMachine example with ROOT output to a file named
"run2". The program is run interactively and the window in Figure appears. From here, the
user types into the visualiser terminal::

  /run/beamOn 1
  exit

,which runs one event and visualises it.

Recreate Mode
=============

After performing a simulation in BDSIM, it is possible to reproduce one or more events exactly
the same again. To do this, the original input gmad files (and any associated external geometry
or field maps, e.g. all the input) are required and should be the same as was originally used.
Along with this, a BDSIM ROOT output file is required. The output file is used to load the
random number generator seed states at the start of each event such that the beam and physics
processes will be the same. For example::

  bdsim --file=mymodel.gmad --outfile=run1 --batch --ngenerate=100

  # let's recreate event 87

  bdsim --file=mymodle.gmad --outfile=selectevent --batch --ngenerate=1 --recreate=run1.root --startFromEvent=87

Usually, one would reproduce something individually to investigate something unexpected and therefore
might change some physics options. If this is the case, subsequent events will still begin with the same
beam coordinates as the random number generator is loaded from the previous simulation output file
at the start of each event. If the physics proceeds differently, this can advance the random number
generator to a different state.

* If the recreation goes beyond the stored number of events, the random number generator will proceed
  as normal. e.g. starting from event 80/100 and generating 30 events, will result in 10 new events.
* Executable options override whatever options were used (and therefore stored in the output) in the
  initial run of BDSIM.
