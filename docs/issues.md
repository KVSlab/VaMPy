# Known issues

## Anaconda compatability issue

Some users may experience errors regarding compatibility if Anaconda already has been configured with certain channels.
To resolve this issue you can remove conflicting channels using:

``` console
$ conda config --remove channels [CHANNEL NAME]
```

Alternatively, you can set your Anaconda channel priority to *flexible*, with the following command:

``` console
$ conda config --set channel_priority flexible
```
