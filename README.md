[main.pdf](latex/main.pdf) file found in latex/

[Overleaf link](https://www.overleaf.com/read/mgwqrhrwhphc#dfb5b3) updated less often, close to compile timeout.

# Chapter 2: Simulations
[main.c](script/src/main.c) solves the equations for an orbit of a massive particle. To try the programm go to 
the script folder and compile it
```
cd script && make main
```
Run the script without any arguments the first time
```
./main.x
```
Follow the instructions to choose the parameters for the simulation. \
The data is saved in ```data/orbit.csv```, run the python script to automatically plot it
```
python3 ch2_plots.py
```

## Some examples
**Unbound orbit**: ```./main.x 10 3.5``` \
![undound_orbit](https://github.com/user-attachments/assets/9c4519c2-1a2c-466a-96ca-5604836b8fd1)


**Bound orbit**: ```./main.x 3 -0.0032 -t 15000 -r 1000``` \
![precession](https://github.com/user-attachments/assets/c073ee3b-7f62-4403-bc6b-6c6e93f49365)


**Infall**: ```./main.x 3 0.2``` \
![infall](https://github.com/user-attachments/assets/2708764d-bca2-4a54-94aa-73bf0f893f0e)

