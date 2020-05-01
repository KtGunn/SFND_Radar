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

## Radar Signature

This section explains the simulation of a moving target--intial range and speed, assumed constant is set--and the signals transmitted and received from it using our FMCW design above.

The signal of a FMCW transmitter is:

_Tx = cos (2*pi*(fc + Slope*t/2)*t_

The return signal is the same except for a shift in time, t-&tau

_Rx = cos (2*pi*(fc + Slope*(t-Tau)/2)*(t-Tau)_

Tau simly is the time is takes the signal to make the roundtrip from the tranmitter/receiver to the target

_Tau = 2*range/c_

From the code:

```
for i=1:length(t)

    % For each time stamp update the Range of the Target for constant velocity.
    t = i*Tchirp/Nr;
    r = ro + vo*t;

    % For each time sample we need update the transmitted and
    % received signal.
    Tx(i) = cos(2*pi*(fc + 0.5*Slope*t)*t);

    td = t - (2*r/c);  % signal round time to transmitter
    Rx(i) = cos(2*pi*(fc + 0.5*Slope*td)*td);

    % Now by mixing the Transmit and Receive generate the beat signal
    % This is done by element wise matrix multiplication of Transmit and

    % Receiver Signal
    Mix(i) = Tx(i).*Rx(i);  % order is immaterial

end
```
