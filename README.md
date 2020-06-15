# Aneurysm_workflow

This is a collection of scripts to run an aneurysm problem with Oasis. First use the automatedPreProcessing to create a mesh
and boundary conditions, then install Oasis, and you can run the problem with:

FEniCS 2019.1.0 or later available:
```
git clone https://github.com/mikaem/Oasis
cd Oasis
python setup.py install  # (or "pip install .", add "--user" if you are on a cluster)
cd ..
````
vmtk available:
```
git clone https://github.com/KVSLab/Aneurysm_workflow.git
cd Aneurysm_workflow/automatedPreProcessing
python automatedPreProcessing.py -m diameter -i ../test/Case_test_71.vtp --aneurysm False -c 1.3
cd ..
```

FEniCS 2019.1.0 or later available:
```
oasis NSfracStep problem=Artery mesh_path=test/Case_test_71.xml.gz
```
