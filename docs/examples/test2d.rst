.. _test2d:

**********************************
Example Problem in 2D
**********************************

To illustrate how to parameters in a text file, here is an example problem ``test2d.in`` (included in the ``problems`` directory): ::

    [fdfault.problem]
    test2d
    data/
    1000
    0
    0
    0.3
    50
    4

    [fdfault.domain]
    2
    2
    801 802 1
    1 2 1
    801
    401 401
    1
    1
    slipweak
    4
    elastic

    [fdfault.fields]
    0. 0. 0. 0. 0. 0.
    none
    none

    [fdfault.block000]
    2.67 32.04 32.04
    0. 0.
    40. 20.
    absorbing
    absorbing
    absorbing
    none
    none
    none
    none
    none

    [fdfault.block010]
    2.67 32.04 32.04
    0. 20.
    40. 20.
    absorbing
    absorbing
    none
    absorbing
    none
    none
    none
    none

    [fdfault.operator]
    0.

    [fdfault.interface0]
    y
    0 0 0
    0 1 0

    [fdfault.friction]
    2
    constant 0. 0. 0. 0. 0. -120. 70. 0.
    boxcar 0. 20. 1.5 0. 0. 0. 11.6 0.
    none

    [fdfault.slipweak]
    3
    constant 0. 0. 0. 0. 0. 0.4 0.677 0.525 0. 0. 0.
    boxcar 0. 2.5 2.5 0. 0. 0. 0. 0. 0. 0. 0.
    boxcar 0. 37.5 2.5 0. 0. 0. 0. 0. 0. 0. 0.
    none

    [fdfault.outputlist]
    V
    V
    0 1000 100
    0 800 1
    401 401 1
    0 0 1
    S
    S
    0 1000 100
    0 800 1
    401 401 1
    0 0 1
    U
    U
    0 1000 100
    0 800 1
    401 401 1
    0 0 1


    [fdfault.frontlist]
    0

