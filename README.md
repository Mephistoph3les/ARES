# Electron CT Multiple Scattering Validation Data Analysis

This project encompasses the data analysis of electron CT multiple scattering validation data collected from 02.04.24 to 05.04.24 at ARES.

## Project Overview

The primary objective of this measurement campaign was to establish a fine calibration for multiple scattering material reconstruction. The focus was on obtaining high statistical data on multiple scattering across various known materials. A finely calibrated scattering width provides information about the material responsible for the scattering. The calibration, found to be independent of the material, characterizes the beam broadening based on the material budget.

## Steps Involved

### 1. Burj Calibration

The first step involved the utilization of a pyramid-shaped block of PEEK (Polyether ether ketone) named Burj. This block featured well-defined steps of increasing widths at 5mm intervals. Positioned within an electron pencil beam of 154.5 MeV energy, a moving stage facilitated the measurement of scattering width for different thickness steps of Burj. Each thickness was measured for 5 minutes, with data stored in separate files.

### 2. Metal Sheets Analysis

Metal sheets of aluminum and nickel, varying in thickness, were mounted on a holding mechanism. These sheets were moved by a stage to enable the electron beam to penetrate them and their overlap regions. High statistics were gathered, including an overnight scan of nickel sheets and a 1-hour scan of aluminum. Data from these scans is stored in two distinct files, necessitating separation based on the corresponding timestamps in the ARES moving stage logfile. It's important to consider not only the 5-minute measurement time per stage position but also the variable moving time of the stage to the next position.

## Conclusion

This project lays the foundation for understanding multiple scattering of electrons. The collected data promises valuable insights into material reconstruction in electron CT imaging.
