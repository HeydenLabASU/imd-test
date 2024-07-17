# imd-test
Experiments to connect a running molecular dynamics simulations via the Interactive Molecular Dynamics Protocol (IMD) to MDAnalysis.

## Definitions
- the molecular dynamics simulation (we provide an example setup for a simulation of a small protein in water with GROMACS) produces data and wait for a TCP/IP socket connection to send this data to a consumer
- the consumer is our python script, which:
  - performs a handshake with the producer to determine endianess of the data
  - sends a GO signal to start the simulation
  - sets the transmission frequency to 10 (ecvery 10 time frames of the simulation)
  - receives energies and coordinates from the producer
- the consumer script currently saves the received coordinates in the MDAnalysis u.atoms.positions array and writes it to file 'imd-test.trr'
Note: the time and index of the trajectory frames are all 0 at the moment

## Requirements
- GROMACS 2018 or later
- MDAnalysis

## Compiling edited imd.cpp
- copy imd.cpp to gromacs-XXXX.X/src/gromacs/imd
- compile GROMACS again to apply changes to the source code

## Testing
- download repository
- open two terminals

Terminal 1:
- navigate to folder imd-test/GMX
- run simulation
```
chmod +x run.sh; ./run.sh
```

Terminal 2:
- navigate to imd-test
- start script:
```
python consumer_IMD.py
```
