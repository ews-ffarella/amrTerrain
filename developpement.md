# Developpement

We will use `uv` for local developpement.

## Linting

Linting the code:

```bash
uvx ruff check --fix --show-fixes .
```

Code formatter:

```bash
uvx ruff format .
```


## APT packages

We need gdal-bin for local developpement

```bash
sudo apt install gdal-bin -yqq
```

## Running the backend using YAML

```bash
uv run amrTerrain yaml_files/3_terrain_rans.yaml
```


## Jupyter notebook

Please refer to this [page](https://docs.astral.sh/uv/guides/integration/jupyter/#using-jupyter-within-a-project) for more information.

```bash
uv add --dev ipykernel
```

Register our kernel

```bash
uv run ipython kernel install --user --env VIRTUAL_ENV $(pwd)/.venv --name="amr-terrain"
```

From there, start the server with (see [here](https://github.com/EWS-Consulting-Private/edv-jupyterhub/blob/main/docker_images/dockerfiles/Dockerfile.jupyterhub)):

```bash
uv run \
    --with jupyter \
    --with notebook \
    --with nbclassic \
    --with jupyter-server-proxy \
    --with jupyter-trame-proxy \
    --with trame \
    --with trame-vtk \
    --with trame-vuetify \
    --with trame-jupyter-extension \
    --with jupyter-resource-usage \
    --with spellchecker \
    --with jupyterlab-spreadsheet-editor \
    --with jupyterlab-quickopen \
    --with jupyterlab-cell-flash \
    --with jupyterlab-filesystem-access \
    --with jupyterlabcodetoc \
    --with jupyterlab-topbar-text \
    --with jupyterlab-logout \
    --with jupyterlab-theme-toggler \
    --with jupyterlab-lsp \
    --with jupyterlab_iframe \
    --with jupyterlab_code_formatter \
    --with jedi-language-server \
    --with black \
    --with ruff \
    --with isort \
    jupyter lab \
     --NotebookApp.port=8888  \
     --NotebookApp.ip=0.0.0.0  \
     --NotebookApp.allow_origin='*'  \
     --NotebookApp.token=''  \
     --NotebookApp.password=''  \
     --no-browser  \
     --notebook-dir="./notebooks"
```
