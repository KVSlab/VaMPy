(post:quantities)=
# Computed quantities

The following tables list the different quantities computed by the post-processing scripts implemented in VaMPy.
In table {numref}`hemo` we present the hemodynamic indices computed by the script `compute_hemodynamic_indices.py`.

```{table} Hemodynamic indices 
:name: hemo
| Quantity                              | Abbreviation/Symbol | Definition                                                        | Unit                                     |
|---------------------------------------|---------------------|-------------------------------------------------------------------|------------------------------------------|
| Wall shear stress                     | WSS, $\tau$         | $\displaystyle \mu \frac{\partial u}{\partial n}$                 | [Pa]                                     |
| Time averaged wall shear stress       | TAWSS               | $\displaystyle \frac{1}{T}\int_0^T \left\| \tau \right\| \, d t$  | [Pa]                                      |
| Temporal wall shear stress gradient   | TWSSG               | $\displaystyle \frac{1}{T}\int_0^T \left\| \frac{\partial \tau}{\partial t} \right\| \,d t$                                                | [Pa/s] |
| Oscillatory shear index               | OSI                 | $\displaystyle \frac{1}{2}\left(1- \frac{\left\| \int_0^T \tau \,d t \right\|}{\int_0^T \left\| \tau \right\|\,d t} \right)$ | [ - ]                                     |
| Relative residence time               | RRT                 | $\displaystyle \frac{1}{(1-2\cdot \text{OSI})\cdot \text{TAWSS}}$ | [1/Pa]                                   |
| Endothelial cell activation potential | ECAP                | $\displaystyle \frac{\text{OSI}}{\text{TAWSS}}$                   | [1/Pa]                                   |
```

In table {numref}`cfd` we present the fluid mechanical metrics and simulation quantities that are computed by the script `compute_flow_and_simulation_metrics.py`.

```{table} Flow and simulation metrics
:name: cfd

| Quantity                          | Abbreviation/Symbol       | Definition                                                                              | Unit           |
|-----------------------------------|---------------------------|-----------------------------------------------------------------------------------------|----------------|
| Velocity                          | $u$                       | $\displaystyle u(x,y,z,t) = (u_x, u_y, u_z)$                                            | [m/s]      |
| Mean velocity                     | $\bar{u}, u_{\text{mean}}$ | $\displaystyle \frac{1}{T} \int_0^T u \,dt$                                             | [m/s]          |
| Turbulent velocity                | $u'$                      | $\displaystyle u - \bar{u}$                                                             | [m/s]          |
| Kinematic viscosity               | $\nu$                     | $\displaystyle \frac{\mu}{\rho}$                                                        | [m$^2$/s]      |
| Time interval                     | $T$                       | $\textit{User defined}$                                                                | [s]            |
| Time step                         | $\Delta t$                | $\displaystyle \frac{T }{ N}$                                                           | [s]            |
| Characteristic edge length        | $\Delta x, h$             | $\displaystyle   \texttt{CellDiameter(mesh)}$                                           | [m]            |
| Courant–Friedrichs–Lewy condition | CFL                       | $\displaystyle \|u\|\frac{\Delta t}{\Delta x}$                                            | [ - ]          |
| Rate of strain                    | $S_{ij}$                  | $\displaystyle  \frac{1}{2} \left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} \right) $   | [1/s]          |
| Turbulent rate of strain          | $s_{ij}$                  | $\displaystyle  \frac{1}{2} \left(\frac{\partial u'_i}{\partial x_j} + \frac{\partial u'_j}{\partial x_i} \right) $ | [1/s]          |
| Absolute rate of strain           | Strain                    | $\displaystyle  \sqrt{\langle S_{ij}, S_{ij}} \rangle $                                 | [1/s]          |
| Dissipation                       | $\mathcal{E} $            | $\displaystyle 2\nu \langle S_{ij}, S_{ij} \rangle $                                    | [m$^2$/s$^3$]  |
| Turbulent dissipation             | $\varepsilon $            | $\displaystyle 2\nu \langle s_{ij}, s_{ij} \rangle $                                    | [m$^2$/s$^3$]  |
| Kinetic energy                    | KE, $E_k$                 | $\displaystyle \frac{1}{2} \left( u_x^2 + u_y^2 + u_z^2  \right)$                       | [m$^2$/s$^2$]  |
| Turbulent kinetic energy          | TKE, $k$                  | $\displaystyle \frac{1}{2} \left( (u'_x)^2 + (u'_y)^2 + (u'_z)^2  \right)$                    | [m$^2$/s$^2$]  |
| Friction velocity                 | $u^{\star},u_\tau$        | $\displaystyle \sqrt{\nu S_{ij}}$                                                       | [m/s]          |
| Generalized length scale          | $\ell^+$                  | $\displaystyle \frac{u^\star \Delta x}{\nu}$                                            | [ - ]          |
| Generalized time scale            | $t^+$                     | $\displaystyle \frac{(u^\star)^2 \Delta t}{\nu}$                                          | [ - ]          |
| Kolmogorov length scale           | $\eta$                    | $\displaystyle \left( \frac{\nu^3}{\varepsilon} \right)^{\frac{1}{4}} $                 | [m]            |
| Kolmogorov time scale             | $\tau_\eta$               | $\displaystyle \left( \frac{\nu}{\varepsilon} \right)^{\frac{1}{2}} $                   | [s]            |
| Kolmogorov velocity scale         | $u_\eta$                  | $\displaystyle \left( \varepsilon \nu  \right)^{\frac{1}{4}} $                          | [m/s]          |

```
