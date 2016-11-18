.. _interface:

**********************************
Interface Input
**********************************

All interfaces are specified by a header with the form ``[fdfault.interfaceN]``, where ``N`` determines the interface number (zero indexed). For problems with :math:`{n}` interfaces, the input file must contain an interface header for all :math:`{n}` interfaces, or you will get an error. Interface headers are followed by the following arguments: ::

    Approximate normal direction
    Minus block indices
    Plus block indices

First, the list defines the approximate normal direction of the interface in the simulation based on the simulation geometry. For rectangular blocks, this is the true normal direction of the block, while if the block has boundaries specified through a file the normal direction may not be precisely in this direction, or the normal direction may not be uniform across the entire block. The direction is set by a string ``x``, ``y``, or ``z``.

Following the direction specification, you must set the indices of the block on the "minus" side of the interface (a list of 3 integers). This can be any block in the simulation, but must be the block with the lowest set of indices that are joined by this interface. Order is not important for setting up the interfaces -- interface 0 can join any pair of neighboring blocks in the simulation -- and the lists can appear in the input file in any order.

Next, the block on the "plus" side of the interface is given by its indices. Because the blocks in the simulation must form a regular grid, the "plus" block must differ in only one index from the "minus" block, and the index that is different must be the same as the direction specified above (this is checked by the code when initializing). For instance, if the minus block is :math:`{(0,0,0)}`, and the direction is ``x``, then the plus block must have the index :math:`{(1,0,0)}`. Line breaks are ignored when reading in the indices.
