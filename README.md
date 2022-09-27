This is the code for the paper ["Spontaneous emergence of fast attractor dynamics in a model of developing primary visual cortex"](https://www.nature.com/articles/ncomms13208), by T. Miconi, J. McKinstry and G. M. Edelman (published in *Nature Communications*, vol. 7, n. 13208, 2016).

This code simulates a small patch of cortex, containing 100 principal neurons
and 20 inhibitory neurons. Excitatory neurons have voltage-based spike-timing
dependent plasticity (Clopath-Gerstner algorithm) at *both* feedforward and
recurrent synapses. Connections to and from inhibitory neurons are fixed and
random. The patch is repeatedly exposed to natural images, preprocessed in a
way that emulates retinal and thelamic filtering. After many iterations, the
network develops recurrent connectivity that supports attractor dynamics
(groups of mutually-connected, similarly-tuned neurons, as observed in cortex),
as well as realistic feedforward receptive fields.



The entire code is contained in one program: stdp.cpp. It requires the Eigen v3 algebra library, which is very easy to install since it's source-only. After dowloading and unzipping the Eigen library into the `$EIGEN_DIR` directory, you can compile stdp.cpp as follows (be sure to adapt to your own compiler options and directory structure):

`g++ -I $EIGEN_DIR/ -O3 -std=c++11 stdp.cpp -o stdp`

The present code base includes pre-trained weights (`w.*` and `wff.*`). If you want to run the full training yourself, just run `./stdp learn` to train a network from scratch. The code is configured to run 500K presentations, which could take about a day. However, good results are evident by ~300K evaluations.

If you want to evaluate the trained network and generate the figures from the paper, you will need to run `./stdp test`, and also `./stdp mix` and `./stdp pulse`. Consult `makefigures.m` for more details.

We include additional Matlab programs that generate the image stimuli (makepatches.m - note that  pre-generated image patches are included in the present code base) and produce figures (showw.m, analyzewlat.m and makefigures.m)

- showw.m produces a picture of the receptive fields (subtracting feedforward weights from the ON and OFF retinal cells).
- analyzewlat.m produces various graphs and preovides information about cluster formation
- makefigures.m makes the figures in the paper. You need to read it and
manually activate the code for each figure, after having produced the necessary
data as instructed in the code for this figure. Note that figures  will look
slightly different from those in the paper, but qualitatively similar.



