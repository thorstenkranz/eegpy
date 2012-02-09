Tutorial
============= 

The intentions of this tutorial are to

* serve as a starting point to learn how to use eegpy
* show some of the key features of eegpy.

The

.. plot::
    :include-source:

    import matplotlib.pyplot as p
    import numpy as np

    p.plot(np.arange(0,10,0.1),
           np.sin(np.arange(0,10,0.1)),
           "rD",label="sine wave")

