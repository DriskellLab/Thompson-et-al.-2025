All scRNA-seq and Stereo-seq datasets are publicly query-able using our webtools on skinregeneration.org (<link>).

Raw and processed data files utilized in this study are publicly available at NCBI GEO (<link>).

Previously-published datasets re-analyzed in this study: Glover et al., 2023 (doi.org/10.1016/j.cell.2023.01.015), Cheng et al., 2018 (doi.org/10.1016/j.celrep.2018.09.006), Solé-Boldo et al., 2020 (doi.org/10.1038/s42003-020-0922-4), and Liu et al., 2022 (doi.org/10.1016/j.devcel.2022.06.005). For NCBI GEO accession numbers for the samples we re-analyzed, please see the Methods section of our manuscript (<link>).

scRNA-seq Analysis in R:
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets 
[7] methods   base     

other attached packages:
 [1] SeuratWrappers_0.4.0        monocle3_1.3.7             
 [3] SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0
 [5] GenomicRanges_1.58.0        GenomeInfoDb_1.42.1        
 [7] IRanges_2.40.1              S4Vectors_0.44.0           
 [9] MatrixGenerics_1.18.0       matrixStats_1.4.1          
[11] shiny_1.10.0                harmony_1.2.3              
[13] Rcpp_1.0.13-1               CellChat_2.1.2             
[15] ggplot2_3.5.1               igraph_2.1.2               
[17] future_1.34.0               patchwork_1.3.0            
[19] viridisLite_0.4.2           limma_3.62.1               
[21] reticulate_1.40.0           tictoc_1.2.1               
[23] Seurat_5.1.0.9006           SeuratObject_5.0.2         
[25] sp_2.1-4                    dplyr_1.1.4                
[27] Biobase_2.66.0              BiocGenerics_0.52.0        


Stereo-seq Analysis in Python:
3.8.20
Package                       Version
----------------------------- --------------
aiohappyeyeballs              2.4.4
aiohttp                       3.10.11
aiosignal                     1.3.1
alabaster                     0.7.13
anndata                       0.9.2
annoy                         1.17.3
anyio                         4.5.2
arboreto                      0.1.6
argon2-cffi                   23.1.0
argon2-cffi-bindings          21.2.0
arrow                         1.3.0
asciitree                     0.3.3
asttokens                     2.4.1
async-lru                     2.0.4
async-timeout                 5.0.1
attrs                         24.2.0
babel                         2.16.0
backcall                      0.2.0
beautifulsoup4                4.12.3
biopython                     1.83
bleach                        6.1.0
blosc2                        2.0.0
bokeh                         2.4.3
boltons                       24.1.0
cell-bin                      1.3.4.1
certifi                       2024.8.30
cffi                          1.17.1
charset-normalizer            3.4.0
click                         8.1.7
cloudpickle                   3.1.0
colorcet                      3.1.0
coloredlogs                   15.0.1
comm                          0.2.2
contourpy                     1.1.1
ctxcore                       0.2.0
cusingler                     1.1.0
cycler                        0.12.1
Cython                        3.0.11
cytoolz                       1.0.0
dask                          2023.5.0
dask-image                    2023.3.0
datashader                    0.15.2
datashape                     0.5.2
debugpy                       1.8.8
decorator                     5.1.1
defusedxml                    0.7.1
dill                          0.3.9
distinctipy                   1.3.4
distributed                   2023.5.0
docrep                        0.3.2
docutils                      0.20.1
exceptiongroup                1.2.2
executing                     2.1.0
fastcluster                   1.2.6
fasteners                     0.19
fastjsonschema                2.20.0
fbpca                         1.0
flatbuffers                   24.3.25
fonttools                     4.55.2
fqdn                          1.5.1
frozendict                    2.4.6
frozenlist                    1.5.0
fsspec                        2024.10.0
gefpy                         1.1.20
geojson                       3.1.0
geosketch                     1.3
get-annotations               0.1.2
gtfparse                      1.2.1
h11                           0.14.0
h5py                          3.8.0
harmonypy                     0.0.6
holoviews                     1.17.1
hotspotsc                     1.1.1
httpcore                      1.0.7
httpx                         0.27.2
humanfriendly                 10.0
hvplot                        0.10.0
idna                          3.10
igraph                        0.11.5
imagecodecs                   2023.3.16
imageio                       2.31.1
imagesize                     1.4.1
importlib_metadata            8.5.0
importlib_resources           6.4.5
inflect                       7.4.0
interlap                      0.2.7
intervaltree                  3.1.0
ipykernel                     6.29.5
ipython                       8.12.3
ipywidgets                    8.1.5
isoduration                   20.11.0
jedi                          0.19.2
Jinja2                        3.1.4
joblib                        1.4.2
joypy                         0.2.6
json5                         0.9.28
jsonpointer                   3.0.0
jsonschema                    4.23.0
jsonschema-specifications     2023.12.1
jupyter                       1.1.1
jupyter_client                8.6.3
jupyter-console               6.6.3
jupyter_core                  5.7.2
jupyter-events                0.10.0
jupyter-lsp                   2.2.5
jupyter_server                2.14.2
jupyter_server_terminals      0.5.3
jupyterlab                    4.2.6
jupyterlab_pygments           0.3.0
jupyterlab_server             2.27.3
jupyterlab_widgets            3.0.13
KDEpy                         1.1.0
kiwisolver                    1.4.7
lazy_loader                   0.4
leidenalg                     0.10.2
llvmlite                      0.39.1
locket                        1.0.0
loompy                        3.0.7
louvain                       0.8.2
lxml                          5.3.0
lz4                           4.3.3
Markdown                      3.7
MarkupSafe                    2.1.5
matplotlib                    3.7.1
matplotlib-inline             0.1.7
matplotlib-scalebar           0.8.1
mistune                       3.0.2
more-itertools                10.5.0
mpmath                        1.3.0
msgpack                       1.1.0
multidict                     6.1.0
multipledispatch              1.0.0
multiprocessing_on_dill       3.5.0a4
natsort                       7.1.1
nbclient                      0.10.0
nbconvert                     7.16.4
nbformat                      5.10.4
nest-asyncio                  1.6.0
networkx                      3.1
notebook                      7.2.2
notebook_shim                 0.2.4
numba                         0.56.4
numcodecs                     0.12.1
numexpr                       2.8.6
numpy                         1.23.5
numpy_groupies                0.9.22
omnipath                      1.0.8
onnxruntime                   1.15.1
opencv-python                 4.8.0.76
overrides                     7.7.0
packaging                     24.2
pandas                        1.5.3
pandocfilters                 1.5.1
panel                         0.14.4
param                         1.13.0
parso                         0.8.4
partd                         1.4.1
patsy                         0.5.6
pexpect                       4.9.0
PhenoGraph                    1.5.7
pickleshare                   0.7.5
pillow                        10.4.0
PIMS                          0.7
pip                           24.2
pkgutil_resolve_name          1.3.10
platformdirs                  4.3.6
plotly                        5.24.1
POT                           0.9.1
prometheus_client             0.21.0
prompt_toolkit                3.0.48
propcache                     0.2.0
protobuf                      5.29.1
psutil                        6.1.0
ptyprocess                    0.7.0
pure_eval                     0.2.3
py-cpuinfo                    9.0.0
pyarrow                       17.0.0
pyCirclize                    1.6.0
pycparser                     2.22
pyct                          0.5.0
Pygments                      2.18.0
pynndescent                   0.5.13
pyparsing                     3.1.4
pyscenic                      0.12.1
python-dateutil               2.9.0.post0
python-json-logger            2.0.7
pytz                          2024.2
pyvips                        2.2.1
pyviz_comms                   3.0.3
PyWavelets                    1.4.1
PyYAML                        6.0.2
pyzmq                         26.2.0
referencing                   0.35.1
requests                      2.32.3
retrying                      1.3.4
rfc3339-validator             0.1.4
rfc3986-validator             0.1.1
rpds-py                       0.20.1
scanorama                     1.7.4
scanpy                        1.9.6
scikit-image                  0.21.0
scikit-learn                  1.3.0
scipy                         1.10.1
seaborn                       0.12.2
Send2Trash                    1.8.3
session_info                  1.0.0
setuptools                    68.2.2
shapely                       2.0.6
six                           1.17.0
slicerator                    1.1.0
slideio                       2.6.5
sniffio                       1.3.1
snowballstemmer               2.2.0
sortedcontainers              2.4.0
soupsieve                     2.6
spatialpandas                 0.4.9
Sphinx                        7.1.2
sphinxcontrib-applehelp       1.0.4
sphinxcontrib-devhelp         1.0.2
sphinxcontrib-htmlhelp        2.0.1
sphinxcontrib-jsmath          1.0.1
sphinxcontrib-qthelp          1.0.3
sphinxcontrib-serializinghtml 1.1.5
SQLAlchemy                    1.3.24
squidpy                       1.2.2
stack-data                    0.6.3
statsmodels                   0.14.1
stdlib-list                   0.10.0
stereopy                      1.3.0
sympy                         1.13.3
tables                        3.8.0
tblib                         3.0.0
tenacity                      9.0.0
terminado                     0.18.1
texttable                     1.7.0
threadpoolctl                 3.5.0
tifffile                      2023.2.3
tinycss2                      1.4.0
tomli                         2.1.0
toolz                         1.0.0
tornado                       6.4.1
tqdm                          4.65.0
traitlets                     5.14.3
typeguard                     4.4.0
types-python-dateutil         2.9.0.20241003
typing_extensions             4.12.2
umap-learn                    0.5.1
uri-template                  1.3.0
urllib3                       2.1.0
validators                    0.34.0
wcwidth                       0.2.13
webcolors                     24.8.0
webencodings                  0.5.1
websocket-client              1.8.0
wheel                         0.44.0
widgetsnbextension            4.0.13
wrapt                         1.17.0
xarray                        0.20.1
yarl                          1.15.2
zarr                          2.16.1
zict                          3.0.0
zipp                          3.20.2
