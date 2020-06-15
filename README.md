# Aneurysm_workflow

This is a collection of scripts to run an aneurysm problem with Oasis. First use the automatedPreProcessing to create a mesh
and boundary conditions, then run the CFD simulation with Oasis, and post-process the results using the , and you can run the problem with:

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

You might run into a problem with vmtk (1.4) if using python 3. To fix this, please [follow these instructions](https://morphman.readthedocs.io/en/latest/installation.html#basic-installation) for fixing the centerline problem. For fixing mesh writing change line 263 to:
```
file = open(self.OutputFileName, 'rb')
````
and line 267 to
```
gzfile = gzip.open(self.OutputFileName, 'wb')
```
Please note that these changes are fixed in the development version of vmtk, but a new version has not been released in a long time.
