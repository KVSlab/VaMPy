# Known issues
## Anaconda compatability issue

Some users may experience errors regarding compatibility if Anaconda
already has been configured with certain channels. To resolve this issue
you can remove conflicting channels using:

``` console
$ conda config --remove channels [CHANNEL NAME]
```

Alternatively, you set your Anaconda channel priority to *flexible*,
with the following command:

``` console
$ conda config --set channel_priority flexible
```

## Issue with vmtkMeshWriter

After installing `morphMan`, you may experience an error
during pre-processing, when the `.vtu` mesh is converted and compressed
into `.xml.gz` format using the `vmtkMeshWriter` method. As a temporary
fix you will need to update the `vmtkMeshWriter` script manually to
avoid this error, located at
`/Users/[USERNAME]/miniconda3/envs/your_environment/lib/python3.10/site-packages/vmtk/vmtkmeshwriter.py`

To apply the fix, open the `vmtkmeshwriter.py` file, navigate to `line
264`, and change:

``` console
file = open(self.OutputFileName,'r')
```

to

``` console
file = open(self.OutputFileName,'rb')
```


