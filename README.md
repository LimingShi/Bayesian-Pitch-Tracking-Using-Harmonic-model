# Bayesian Pitch Tracking Using Harmonic model

A fast pitch tracking algorithm using the harmonic model.

# How to run
Run run_white_example.m (white Gaussian noise) or run_colored_example.m (factory noise) in MATLAB

# Examples
An example under white Gaussian noise:

<p align="center">
<img src=figures/2222.png>
</p>
<center> Figure 1: Pitch estimates under 0 dB white Gaussian noise.</center>

<p align="center">
<img src=figures/1111.png>
</p>

<center>Figure 2: Pitch estimates under 0 dB factory noise.</center>



# How to cite
 Shi L, Nielsen J K, Jensen J R, et al. Bayesian Pitch Tracking Based on the Harmonic Model[J]. IEEE/ACM Transactions on Audio, Speech, and Language Processing, 2019.

# References
This fast computation of the likelihood function is based on the fast pitch estimation method proposed in

Fast fundamental frequency estimation: Making a statistically efficient estimator computationally efficient. Nielsen, Jesper Kjær; Jensen, Tobias Lindstrøm; Jensen, Jesper Rindom; Christensen, Mads Græsbøll; Jensen, Søren Holdt. In: Signal Processing, 135, 2017, pp. 188-197.

Bayesian Model Comparison With the g-Prior. Nielsen, Jesper Kjær; Christensen, Mads Græsbøll; Cemgil, Ali Taylan; Jensen, Søren Holdt. In: IEEE Transactions on Signal Processing, 62 (1), 2014, pp. 225-238.

where the source code is available in
https://github.com/jkjaer/fastF0Nls


This noise PSD tracker is based on the method proposed in

Gerkmann, T. & Hendriks, R. C. Unbiased MMSE-Based Noise Power Estimation With Low Complexity and Low Tracking Delay, IEEE Trans Audio, Speech, Language Processing, 2012, 20, 1383-1393

where the source code is available in
http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
