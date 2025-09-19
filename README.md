# Garachico eruption exercise

This is the CWL-based workflow for the exercise of Garachico to be developed on September 26, 2025. The FALL3D model is driven by GFS data in order to generate a 24-h forecast
of the Garachico eruption started at 10:00 UTC.

It can be executed using the command:

```console
cwltool workflow.cwl small_args.yml
```

for the small domain covering the island or

```console
cwltool workflow.cwl large_args.yml
```

for the large domain.

It runs the FALL3D model using a user-defined number of MPI processes and generates a png and GeoTIFF images as  an output.

# Requirements

## CWL runner

A CWL runner is required for running the CWL workflow. For example, you can install `cwltool`, a Python Open Source project maintained by the CWL community:

```console
pip install cwltool
```

## Docker container

By default the job is executed in a [Docker container][Dockerhub].
If you prefer Podman runtime for running containers, use:
```console
cwltool --podman workflow.cwl arguments.yml
```
or
```console
cwltool --udocker workflow.cwl arguments.yml
```
or
```console
cwltool --singularity workflow.cwl arguments.yml
```
for singularity/apptainer.

<!----------------------------------------------------------------------------->

[Dockerhub]: https://hub.docker.com/r/lmingari/garachico
