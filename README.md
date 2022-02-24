# SPG: signal processing on graphs

This package contains code for constructing signal-adapted systems of spectral kernels for a given graph and a given graph signal set, based on theory presented in the following [paper] (https://doi.org/10.1109/TSP.2016.2591513):

> [1] Behjat et al., 2016. Signal-adapted tight frames on graphs. IEEE Trans. Signal Process.

The Minnesota road graph, the Alameda graph and the cerebellum graph, together with their associated datasets as described in the paper are also included in the package. The cerebellar dataset is partially included due to its large size, but the remainder can also be provided upon request.

The functions and demos also provide means to construct SGWT frame [2], Meyer-like [3] and spectrum-adapted [4] systems of spectral kernels for comparison with the proposed signal-adapted systems of spectral kernels.

> [2] Hammond, et al., 2011. Wavelets on graphs via spectral graph theory, Appl. Comput. Harmon. Anal., vol. 30, pp. 129-150.

> [3] Leonardi, et al., 2013. Tight wavelet frames on multislice graphs. IEEE Trans. Signal Process., vol. 61(13), pp. 3357-3367.  

> [4] D. I. Shuman, et al., 2015. Spectrum-adapted tight graph wavelet and vertex-frequency frames, IEEE Trans. Signal Process., vol. 63(16), pp. 4223-4235.

For installation, simply unpack the directory, and then type in Matlab command window: `spg_setpath`

You may then run the following demos which replicate the some of results presented in [1]: 

`spg_demo_uniformMeyerType`    (cf. Fig. 2 in [1])

`spg_demo_minnesota`                (cf. Fig. 4 in [1])

`spg_demo_alameda`                   (cf. Fig. 6 in [1])

`spg_demo_cerebellum`               (cf. Fig. 9 in [1])

Demos for replicating the other results can also be constructed and provided upon request.


The following demo can be used to construct signal-adapted systems of spectral kernels using your own graph and graph signal set:

`spg_demo_your_data_frame(...)`*
      

The following demo can be used to construct signal-adapted systems of spectral kernels using your own graph and graph signal set, and to decompose a set of graph signals using the constructed frame (not necessarily the same signals used to construct the frame):

`spg_demo_your_data_decompose(...)`*

*This demo requires inputs; read the demo's help function for instructions.

