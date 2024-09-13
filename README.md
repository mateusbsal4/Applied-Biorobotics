# Applied Biorobotics

This repository contains the implementation and reports with the corresponding results for the assignments of the Applied Biorobotics course at TUM. The course covered the fundamentals of bipedal locomotion, bipedal robot control, and bio-inspired optimization.


## Assignment 2 - The Bouncing Mass Problem

The task involved deriving the equations of motion for a mass dropped onto a spring and calculating mechanical energy contributions at each time step. The solution required deriving symbolic EoMs for both flight and stance phases, and solving numerically using MATLAB. 

## Assignment 6 - The Spring-Mass Walker

The task involved deriving the equations of motion for a spring-mass walker and implementing the system in MATLAB/Simulink, with a focus on gait phases, touchdown/takeoff conditions, and the angle of attack. The report includes the resulting trajectory, forward speed, spring forces, phase plot, and mechanical energy contributions.

## Assignment 8 - Neural Network Based Controller

The task involved programming a neural network-based controller from Geng 2006 for the multibody JenaFox robot model in MATLAB/Simulink. After validating the controller with sample data, it was tested in the simulation. The report includes a summary of neural control, network architecture, adaptations due to the ISB convention, implementation challenges, and an analysis of parameters influencing the robot's average forward velocity.

## Assignment 9 - Particle Swarm Optimization

The task involved implementing Particle Swarm Optimization (PSO) to maximize the walking speed of the JenaFox multibody simulation using the neural network controller from Assignment 8, by optimizing the parameters PEA_hip, AEA_hip and GM_h (hip motor gain).
