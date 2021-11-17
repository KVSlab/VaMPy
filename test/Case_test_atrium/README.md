# Left atrium model
The provided left atrium model is collected from the Left Atrial Segmentation Challenge 2013 [1], spesifically case b003, segmented by INRIA. 

The model meshing and probe points can be computed using automatedPreProcessing with the following command:

```
python automatedPreProcessing/automatedPreProcessing.py -m diameter -i test/Case_test_atrium/atrium.vtp -c 1.3
```

Note: Aneurysm workflow is intended for artery models. Therefore, boundary conditions for other models will be wrong. Adapting to other vascular models, such as the atrium, is work in progress.  

---

[1] Tobon-Gomez, Catalina; Geers, Arjan J.; Peters, Jochen; Jürgen Weese; Karen Pinto; Rashed Karim; et al. (2015): Left Atrial Segmentation Challenge 2013: MRI results. figshare. Dataset. https://doi.org/10.6084/m9.figshare.1492974.v3  

---
