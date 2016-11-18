.. _slipweak:

**********************************
Slip Weakening Input
**********************************

For slip weakening laws, there are six parameters: ::

    dc mus mud c0 tc trup

``dc`` is the slip weakening distance, ``mus`` is the static friction coefficient, ``mud`` is the dynamic friction coefficient, ``c0`` is the frictional cohesion, ``tc`` is the characteristic failure time, and ``trup`` is the forced rupture time. The first four parameters are fairly standard, while the final two parameters are used to impose kinematic forcing on the rupture front: ``tc`` is the time scale over which the friction weakens from static to dynamic, and ``trup`` is the time at which this process initiates. If ``trup`` is 0, then no kinematic forcing is used.