.. _text:

**********************************
Text Input Files
**********************************

Text input files contain multiple sections. With one exception, the order of the sections is not important -- the file is read multiple times by various parts of the code, and each time the file is read the code starts from the beginning and scans through to find the appropriate section. Each section is designated by a text string ``[fdfault.<name>]``, where ``<name>`` refers to which part of the code will be reading this section. If the code expects to find a certain section and it reaches the end of the text file without finding it, the code will abort. At minimum, the following sections are required to fully specify a problem: ::

    [fdfault.problem]
    [fdfault.domain]
    [fdfault.fields]
    [fdfault.blockXYZ]
    [fdfault.outputlist]
    [fdfault.frontlist]

In addition to these sections, optional sections describe additional blocks and interfaces, as well as other parameter settings. These include the following sections: ::

    [fdfault.cartesian]
    [fdfault.operator]
    [fdfault.interfaceN]
    [fdfault.friction]
    [fdfault.slipweak]
    [fdfault.stz]

If the problem has more than one block or more than one interface, the sections are designated with the numeric value in place of ``XYZ`` or ``N`` included in the section header.

.. toctree::
   :maxdepth: 1

   problem
   domain
   cartesian
   fields
   operator
   block
   interface
   friction
   slipweak
   stz
   outputlist
   frontlist