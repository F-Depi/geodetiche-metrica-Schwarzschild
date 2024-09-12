Batchelor degree thesis in Physics from the University of Trento.

The thesis can be found
[Geodesics in Schwarzschild Metric, F. De Paoli](latex/main.pdf).

# Chapter 1
Theoretical introduction to the Schwarzschild metric and
the geodesic equations for massive and massless particles.

The content is mostly based on the books Gravity: an introduction to Einstein's
general relativity by James B. Hartle and Black Holes, White Dwarfs and Neutron
Stars by Stuart L. Shapiro and Saul A. Teukolsky. 

# Chapter 2: Simulations
The second chapter is a numerical simulation of the geodesic equations for
a massive particle.

[main.c](script/src/main.c) solves the equations for an orbit of a massive particle. To try the programm go to 
the script folder and compile it
```
cd script && make main
```
Run the script without any arguments the first time
```
./main.x
```
Follow the instructions to choose the parameters for the simulation.

The data can be plotted with the python script ```plot_orbit.py``` by giving l
and E as input arguments.

Alternatively, the script ```sim.sh``` can be used to run the simulation,
immediately plot the data and then get prompted to either save the data in a
new folder in ```data/keep``` or delete it.

```ani_orbit.py``` can be used to create an animation of the orbit.

## Some examples
**Unbound orbit**: ```./main.x 10 3.5``` \
![undound_orbit](https://github.com/user-attachments/assets/9c4519c2-1a2c-466a-96ca-5604836b8fd1)


**Bound orbit**: ```./main.x 3 -0.0032 -t 15000 -r 1000``` \
![precession](https://github.com/user-attachments/assets/c073ee3b-7f62-4403-bc6b-6c6e93f49365)


**Infall**: ```./main.x 3 0.2``` \
![infall](https://github.com/user-attachments/assets/2708764d-bca2-4a54-94aa-73bf0f893f0e)

**Spirograph**: ```./sim.sh 3 -0.005 -t 116400 -h 0.1``` (a lot of data!)
![image](https://github.com/user-attachments/assets/031b80e1-d4bf-4ecd-92b4-86765a9c6181)
