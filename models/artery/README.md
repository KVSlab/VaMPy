# Carotid artery model
The provided carotid artert model is collected from the Aneurisk database [1], and is a modified version of case C0015.
The copy of the full dataset is located [here](https://github.com/hkjeldsberg/AneuriskDatabase).

The model meshing and probe points can be computed using automatedPreProcessing with the following command:

```
vampy-mesh -m diameter -i models/artery/artery.vtp -c 1.3
```


---

[1] Sangalli, L. M., Secchi, P., & Vantini, S. (2014). AneuRisk65: A dataset of three-dimensional cerebral vascular geometries.

---
