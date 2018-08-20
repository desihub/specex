# Installing specex

## Prerequisites

This package depends only on [HARP](https://github.com/tskisner/HARP),
which provides links to `cfitsio` and `lapack`.  The command `harpconfig`
must be in your `$PATH`.


## Installing with Makefile

Ensure that HARP is installed and that the harpconfig script is in your $PATH.

Next decide where to install specex.  We will call that directory `$WHERE`.
Make all software and install it by doing:

`SPECEX_PREFIX=$WHERE make install`

In order to use the specex plugin with HARP, you need to append location of
the plugin to the search path that HARP uses by doing:

`export HARP_PLUGIN_PATH=$WHERE/lib`


## Tagging

Update the tag version with
```
make version TAG=x.y.z
```

To automatically increment a dev tag:
```
make version
```

