# SFND_Radar
SensorFusionNanoDegree Radar Signal Processing Project

## Introduction

The purpose of the project was to use MATLAB to first simulate radar signals from a moving vehicle, then process those signals to retrive the range and speed. The project follows these steps:

1. Calculate the parameters of the FMCW (Frequency Modulated Continuous Wave) waveform given requirements;
2. Set range and speed of a target and simulate its radar signal;
3. Generate the radar signature by mixing the transmitted and returned signals to generate the beat frequency;
4. Perform range FFT on the beat signal to determine the range set in step 2;
5. Perform 2D range and doppler FFT to reveal the range and velocity;
6. Perform CFAR (Constant False Alarm Rate) analysis to track the target in the 2D FFT map from step 5.

The sequence of operations is shown in the figure below

<img src="./images/project_flow.png">

The matlab code to perform the operations is included in the directory 'matlab'. The code is contained in the file 'RadarProject.mlx.' Parts of that file will be copied to this page and discussed to explain how the above steps were performed.

## FMCW Signal Design

<img src="./images/fmcw.png">

The figure above shows the characteristics of the FMCW radar signal. Over the time span 'Ts' the transmitted signal frequency is ramped up from a base carrier frequency (fc), then dropped and ramped up again. The sawtooth frequency ramps are called 'chirps.' A chirp is defined by its 'slope' of chirp frequency bandwidth 'B' over chirp time 'Ts', 'slope=B/Ts.' The range to a target is determined by the frequency difference between emitted and received frequencies, shown as 'fb' or beat frequency. The generation of the beat frequency is done by subtracting received from transmitted signal, an operation called 'mixing.'

The chirp bandwidth is related to the required range resolution; smaller the resolution, higher the bandwidth. Chirp time is related to the maximum range of a target; greater the range, longer the sweep time.

From the code:

```
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_resolution = 1; % range resolution
r_max = 200;      % maximum radar range
c = 3e8;          % speed of light

B = c/r_resolution/2;    % sweep height
Tchirp = 5.5*2*r_max/c;  % 5.5*max_round_time=5.5*(2*rMax/c)
Slope = B/Tchirp;
```

