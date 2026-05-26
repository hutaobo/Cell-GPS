# -*- mode: python ; coding: utf-8 -*-

from pathlib import Path
from PyInstaller.utils.hooks import collect_submodules, copy_metadata


project_root = Path(SPECPATH)
metadata_datas = []
for package_name in ("anndata", "pyXenium", "scanpy", "h5py", "pyarrow"):
    metadata_datas += copy_metadata(package_name)

a = Analysis(
    [str(project_root / "src" / "sfplot" / "gui" / "gui_app.py")],
    pathex=[str(project_root), str(project_root / "src")],
    binaries=[],
    datas=metadata_datas,
    hiddenimports=[
        "cellgps.analysis.searcher_findee_score",
        "cellgps.preprocessing.data_processing",
        "sfplot.analysis.searcher_findee_score",
        "sfplot.preprocessing.data_processing",
        "pyXenium.io",
    ]
    + collect_submodules("pyXenium")
    + collect_submodules("pygments"),
    hookspath=[str(project_root)],
    hooksconfig={
        "matplotlib": {
            "backends": ["Agg"],
        },
    },
    runtime_hooks=[str(project_root / "my_startup_hook.py")],
    excludes=[
        "keras",
        "pytest",
        "pycirclize",
        "sfplot.TCR_distance_heatmap",
        "sfplot.binned_analysis",
        "sfplot.binned_analysis_gpu",
        "sfplot.circle_heatmap",
        "sfplot.circular_dendrogram",
        "sfplot.compute_col_dendrogram_scores",
        "sfplot.compute_cophenetic_distances_from_df_memory_opt",
        "sfplot.plotting",
        "sfplot.sfplot",
        "sfplot.tbc_analysis",
        "sfplot.tbc_analysis_serial",
        "sfplot.visium_preprocesssing",
        "sfplot.xenium_preprocessing",
        "tensorflow",
        "torch",
    ],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="cellgps",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    icon=str(project_root / "cell_gps_icon.ico"),
)
