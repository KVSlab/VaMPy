# Left atrium model
The provided left atrium model is collected from the Left Atrial Segmentation Challenge 2013 [1], spesifically case b003, segmented by INRIA. 

The model meshing and probe points can be computed using automatedPreProcessing with the following command:

```
python automatedPreProcessing/automatedPreProcessing.py -m diameter -i models/atrium/atrium.vtp -c 1.3
```

Note: Boundary conditions for other than vascular models are still work in progress.

---

[1] Tobon-Gomez, Catalina; Geers, Arjan J.; Peters, Jochen; Jürgen Weese; Karen Pinto; Rashed Karim; et al. (2015): Left Atrial Segmentation Challenge 2013: MRI results. figshare. Dataset. https://doi.org/10.6084/m9.figshare.1492974.v3  

---
