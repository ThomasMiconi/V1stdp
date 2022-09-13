This is the code for the upcoming paper, "Spontaneous emergence of fast attractor dynamics in a model of developing primary visual cortex", by T. Miconi, J. McKinstry and G. M. Edelman, to appear in Nature Communications.

This code simulates a small patch of cortex, containing 100 principal neurons
and 20 inhibitory neurons. Excitatory neurons have voltage-based spike-timing
dependent plasticity (Clopath-Gerstner algorithm) at *both* feedforward and
recurrent synapses. Connections to and from inhibitory neurons are fixed and
random. The patch is repeatedly exposed to natural images, preprocessed in a
way that emulates retinal and thelamic filtering. After many iterations, the
network develops recurrent connectivity that supports attractor dynamics
(groups of mutually-connected, similarly-tuned neurons, as observed in cortex),
as well as realistic feedforward receptive fields.



The entire code is contained in one program: stdp.cpp. It requires the Eigen v3 algebra library, which is very easy to install since it's source-only. After dowloading and unzipping the Eigen library, you can compile stdp.cpp as follows (be sure to adapt to your own compiler options and directory structure):

`g++ -I $EIGEN_DIR/Eigen/ -O3 -std=c++11 stdp.cpp -o stdp`

Just run `./stdp learn latconnmult 4.0 wie .6 wpenscale .33 altpmult .75 delayparam 5.0` (as specified in the COMMAND.sh file) to build a network. The code is configured to run the full 1M presentations, which could take several days. However, good results are evident by ~300K evaluations.

We also include additional Matlab programs that generate the image stimuli (makepatches.m) and produce figures (showw.m, analyzewlat.m and makefigures.m)

- showw.m produces a picture of the receptive fields (subtracting feedforward weights from the ON and OFF retinal cells).
- analyzewlat.m produces various graphs and preovides information about cluster formation
- makefigures.m makes the figures in the paper. You need to read it and
manually activate the code for each figure, after having produced the necessary
data as instructed in the code for this figure. Note that figures  will look
slightly different from those in the paper, but qualitatively similar.



