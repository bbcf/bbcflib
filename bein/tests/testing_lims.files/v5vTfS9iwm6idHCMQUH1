=======================
Documentation: ``bein``
=======================

.. automodule:: bein

Executions
**********

.. autofunction:: execution
.. autoclass:: Execution
.. automethod:: Execution.add
.. automethod:: Execution.use
.. _minilims:

MiniLIMS
********

.. autoclass:: MiniLIMS
.. automethod:: MiniLIMS.add_alias
.. automethod:: MiniLIMS.associate_file
.. automethod:: MiniLIMS.associated_files_of
.. automethod:: MiniLIMS.copy_file
.. automethod:: MiniLIMS.delete_alias
.. automethod:: MiniLIMS.delete_execution
.. automethod:: MiniLIMS.delete_file
.. automethod:: MiniLIMS.delete_file_association
.. automethod:: MiniLIMS.export_file
.. automethod:: MiniLIMS.fetch_execution
.. automethod:: MiniLIMS.fetch_file
.. automethod:: MiniLIMS.import_file
.. automethod:: MiniLIMS.path_to_file
.. automethod:: MiniLIMS.resolve_alias
.. automethod:: MiniLIMS.search_executions
.. automethod:: MiniLIMS.browse_executions
.. automethod:: MiniLIMS.search_files

Programs
********

.. autoclass:: program

Miscellaneous
*************

.. autofunction:: unique_filename_in

.. autoclass:: bein.ProgramOutput

.. attribute:: return_code

  An integer giving the return code of the program when it exited.
  By convention this is 0 if the program succeeded, or some other
  value if it failed.

.. attribute:: pid

  An integer giving the process ID under which the program ran.

.. attribute:: arguments

  A list of strings which were the program and the exact arguments
  passed to it.

.. attribute:: stdout

  The text printed by the program to ``stdout``.  It is returned
  as a list of strings, each corresponding to one line of
  ``stdout``, and each still carrying their terminal ``\n``.

.. attribute:: stderr

  The text printed by the program to ``stderr``.  It has the same
  format as ``stdout``.

.. autoexception:: ProgramFailed

.. autofunction:: task
