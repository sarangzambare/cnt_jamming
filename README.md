# cnt_jamming
This project explores the jamming behaviour of carbon nanotubes in solution phase, when its concentration is uniformly increased under various initial conditions

The tubes were considered to be sphereocylinders with helical charge distribution.

The jamming was simulated in cubes with periodic boundary conditions.  

The parameters ***pressure, charge per rod*** and ***salt concentration*** were varied so as to study their impact on the jamming densities. The system used for these simulations consisted of **2116 rods**, which was large enough to guarantee that there are enough rods (more than 2) along each dimension of the box, so that the assumed periodicity is valid.

## Observations:
For a given set of parameters, the system was compressed from the isotropic low density phase, until there was no appreciable increase in the density. At this point, the density was noted and plotted.

### Pressure variation:

Five values of pressure were chosen, 50, 100, 200, 300, 500 (mPa) The density curves are plotted below, (in rods/cm^3)

![alt text](https://raw.githubusercontent.com/sarangzambare/cnt_jamming/master/png/press_var.png)

The table below summarises the density values:

![alt text](https://raw.githubusercontent.com/sarangzambare/cnt_jamming/master/png/press_var_tab.png)

Here the parameter Z is the charge per rod in mC.

### Charge per rod variation:

Five values of charge per rod were chosen, 100,200,300,400,500 (mC) The density curves are plotted below.

![alt text](https://raw.githubusercontent.com/sarangzambare/cnt_jamming/master/png/cpr_var.png)

The table below summarises the density values:

![alt text](https://raw.githubusercontent.com/sarangzambare/cnt_jamming/master/png/cpr_var_tab.png)


### Salt concentration variation:

Five values of salt concentration were used, 0.0, 2.0, 4.0, 6.0, 8.0. The density curves are plotted below. (mM)


![alt text](https://raw.githubusercontent.com/sarangzambare/cnt_jamming/master/png/sc_var.png)

The table below summarises the density values:

![alt text](https://raw.githubusercontent.com/sarangzambare/cnt_jamming/master/png/sc_var_tab.png)

## Study of clusters:

The simulations were carried out using a starting state which is sparse, that is the rods were far from jamming. The snapshot below captures the initial state:

![alt text](https://raw.githubusercontent.com/sarangzambare/cnt_jamming/master/png/unjammed.jpg)

In the simulations that were carried out, it was almost always observed that the final jammed configuration of rods comes in the form of clusters of rods, with rods of same cluster pointing in the same direction. In the snapshot below, rods of same color point along the same vector, hence its a measure of how aligned the rods are.


![alt text](https://raw.githubusercontent.com/sarangzambare/cnt_jamming/master/png/jammed.png)



An effort was made to tweak the program so that it calculates the cluster sizes. In the ***infile*** , I have introduced two parameters, namely Distance and MinScalPro.
Distance is the minimum distance at which two rods have to be, to be called as neighbours and MinScalPro is the minimum dot product of the orientations of the two rods, above which they will be counted as having the same orientation (ranges from 0 to 1).
It turns out that the number of clusters counted, and the cluster size distribution is very sensitive to the parameters Distance and MinScalPro. These parameters should be changed according to the definition of the user’s perception of a cluster.
Since its difficult to tell whether the outcome of the program is correct or not just by looking at a snapshot of rods, further work in this field was not carried out, due to lack to time. Though it would be really great if one can find out the correct set of values for the parameters Distance and MinScalPro to get correct cluster size distributions.
