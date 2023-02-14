# Pre-processing

## Meshing and boundary conditions

In this brief introduction, we present VaMPy's functionality by preparing, simulating, and post-processing models of
the [internal carotid artery](https://en.wikipedia.org/wiki/Internal_carotid_artery) and
the [left atrium](https://en.wikipedia.org/wiki/Atrium_(heart)), although the software may be readily used for other
tubular or vascular shapes.

The first step of using the Vascular Modeling Pypeline is pre-processing. The pre-processing scripts are located inside
the
`automatedPreprocessing` folder, and can be run by using the `vampy-mesh` command in the terminal. The script generates
a mesh, boundary conditions, and probes for velocity and pressure sampling. Here we will perform pre-processing for the
artery case located in the `models` folder. In this example we pass the input model path (`-i`) and mesh coarsening
factor (`-c`) command-line arguments:

``` console
$ vampy-mesh -i models/artery/artery.vtp -c 1.3
```

When complete, a render window will pop up showing the meshed model, prope points, flow rate and pressure split, as
shown in {numref}`render`. The script will save the volumetric mesh as `artery.vtu`, alongside a compressed DOLFIN mesh
in `artery.xml.gz`, used for the computational fluid dynamics (CFD) simulations. The pre-processing script will also
produce an info file and a probe file, named `artery_info.json` and `artery_probes`, respectively.

```{figure} figures/render.png
---
name: render
---
Visualization of the meshed artery model by running `vampy-mesh`.
```


