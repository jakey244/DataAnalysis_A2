from pathlib import Path
from typing import Callable, List, Union

import pandas as pd
import xarray as xr

from amoc_analysis import analysis
from amoc_analysis.analysis import apply_defaults

# Default list of RAPID data files
RAPID_DEFAULT_SOURCE = "https://rapid.ac.uk/sites/default/files/rapid_data/"
RAPID_TRANSPORT_FILES = ["moc_transports.nc"]
RAPID_DEFAULT_FILES = ["moc_transports.nc"]

# Inline metadata dictionary
RAPID_METADATA = {
    "description": "RAPID 26N transport estimates dataset",
    "project": "RAPID-AMOC 26°N array",
    "web_link": "https://rapid.ac.uk/rapidmoc",
    "note": "Dataset accessed and processed via xarray",
}

# File-specific metadata placeholder
RAPID_FILE_METADATA = {
    "moc_transports.nc": {
        "data_product": "RAPID layer transport time series",
    },
}


@apply_defaults(RAPID_DEFAULT_SOURCE, RAPID_DEFAULT_FILES)
def read_rapid(
    source: Union[str, Path, None],
    file_list: Union[str, list[str]],
    transport_only: bool = True,
    data_dir: Union[str, Path, None] = None,
    redownload: bool = False,
) -> list[xr.Dataset]:
    """Load the RAPID transport dataset from a URL or local file path into an xarray.Dataset.

    Parameters
    ----------
    source : str, optional
        URL or local path to the NetCDF file(s).
        Defaults to the RAPID data repository URL.
    file_list : str or list of str, optional
        Filename or list of filenames to process.
        If None, will attempt to list files in the source directory.
    transport_only : bool, optional
        If True, restrict to transport files only.
    data_dir : str, Path or None, optional
        Optional local data directory.
    redownload : bool, optional
        If True, force redownload of the data.

    Returns
    -------
    list of xr.Dataset
        The loaded xarray datasets with basic inline metadata.

    Raises
    ------
    ValueError
        If the source is neither a valid URL nor a directory path.
    FileNotFoundError
        If no valid NetCDF files are found in the provided file list.

    """
    if file_list is None:
        file_list = RAPID_DEFAULT_FILES
    if transport_only:
        file_list = RAPID_TRANSPORT_FILES
    if isinstance(file_list, str):
        file_list = [file_list]

    local_data_dir = Path(data_dir) if data_dir else analysis.get_default_data_dir()
    local_data_dir.mkdir(parents=True, exist_ok=True)

    datasets = []

    for file in file_list:
        if not file.lower().endswith(".nc"):
            continue

        download_url = (
            f"{source.rstrip('/')}/{file}" if analysis._is_valid_url(source) else None
        )

        file_path = analysis.resolve_file_path(
            file_name=file,
            source=source,
            download_url=download_url,
            local_data_dir=local_data_dir,
            redownload=redownload,
        )

        try:
            ds = xr.open_dataset(file_path)
        except Exception as e:
            raise FileNotFoundError(f"Failed to open NetCDF file: {file_path}: {e}")

        file_metadata = RAPID_FILE_METADATA.get(file, {})
        analysis.safe_update_attrs(
            ds,
            {
                "source_file": file,
                "source_path": str(file_path),
                **RAPID_METADATA,
                **file_metadata,
            },
        )
        if "time" in ds.dims or "time" in ds.coords:
            ds = ds.rename({"time": "TIME"})

        datasets.append(ds)

    if not datasets:
        raise FileNotFoundError(f"No valid RAPID NetCDF files found in {file_list}")

    return datasets


def _get_reader(array_name: str) -> Callable:
    """Return the reader function for the given array name.

    Parameters
    ----------
    array_name : str
        The name of the observing array.

    Returns
    -------
    function
        Reader function corresponding to the given array name.

    Raises
    ------
    ValueError
        If an unknown array name is provided.

    """
    readers = {
        "rapid": read_rapid,
    }
    try:
        return readers[array_name.lower()]
    except KeyError:
        raise ValueError(
            f"Unknown array name: {array_name}. Valid options are: {list(readers.keys())}",
        )


def load_sample_dataset(array_name: str = "rapid") -> xr.Dataset:
    """Load a sample dataset for quick testing.

    Currently supports:
    - 'rapid' : loads the 'RAPID_26N_TRANSPORT.nc' file

    Parameters
    ----------
    array_name : str, optional
        The name of the observing array to load. Default is 'rapid'.

    Returns
    -------
    xr.Dataset
        A single xarray Dataset from the sample file.

    Raises
    ------
    ValueError
        If the array_name is not recognised.

    """
    if array_name.lower() == "rapid":
        sample_file = "moc_transports.nc"
        datasets = load_dataset(
            array_name=array_name,
            file_list=sample_file,
            transport_only=True,
        )
        if not datasets:
            raise FileNotFoundError(
                f"No datasets were loaded for sample file: {sample_file}",
            )
        return datasets[0]

    raise ValueError(
        f"Sample dataset for array '{array_name}' is not defined. "
        "Currently only 'rapid' is supported.",
    )


def load_dataset(
    array_name: str,
    source: str = None,
    file_list: Union[str, List[str], None] = None,
    transport_only: bool = True,
    data_dir: Union[str, Path, None] = None,
    redownload: bool = False,
) -> List[xr.Dataset]:
    """Load raw datasets from a selected AMOC observing array.

    Parameters
    ----------
    array_name : str
        The name of the observing array to load. Options are:
        - 'rapid' : RAPID 26N array
    source : str, optional
        URL or local path to the data source.
        If None, the reader-specific default source will be used.
    file_list : str or list of str, optional
        Filename or list of filenames to process.
        If None, the reader-specific default files will be used.
    transport_only : bool, optional
        If True, restrict to transport files only.
    data_dir : str, optional
        Local directory for downloaded files.
    redownload : bool, optional
        If True, force redownload of the data.

    Returns
    -------
    list of xarray.Dataset
        List of datasets loaded from the specified array.

    Raises
    ------
    ValueError
        If an unknown array name is provided.

    """
    reader = _get_reader(array_name)
    datasets = reader(
        source=source,
        file_list=file_list,
        transport_only=transport_only,
        data_dir=data_dir,
        redownload=redownload,
    )

    _summarise_datasets(datasets, array_name)

    return datasets


def _summarise_datasets(datasets: List[xr.Dataset], array_name: str) -> None:
    """Print a summary of loaded datasets."""
    summary_lines = []
    summary_lines.append(f"Summary for array '{array_name}':")
    summary_lines.append(f"Total datasets loaded: {len(datasets)}\n")

    for idx, ds in enumerate(datasets, start=1):
        summary_lines.append(f"Dataset {idx}:")

        # Filename from metadata
        source_file = ds.attrs.get("source_file", "Unknown")
        summary_lines.append(f"  Source file: {source_file}")

        # Time coverage
        time_var = ds.get("TIME")
        if time_var is not None:
            time_start = pd.to_datetime(time_var.values[0]).strftime("%Y-%m-%d")
            time_end = pd.to_datetime(time_var.values[-1]).strftime("%Y-%m-%d")
            summary_lines.append(f"  Time coverage: {time_start} to {time_end}")
        else:
            summary_lines.append("  Time coverage: TIME variable not found")

        # Dimensions
        summary_lines.append("  Dimensions:")
        for dim, size in ds.sizes.items():
            summary_lines.append(f"    - {dim}: {size}")

        # Variables
        summary_lines.append("  Variables:")
        for var in ds.data_vars:
            shape = ds[var].shape
            summary_lines.append(f"    - {var}: shape {shape}")

        summary_lines.append("")  # empty line between datasets

    summary = "\n".join(summary_lines)

    # Print to console
    print(summary)


def save_dataset(
    ds: xr.Dataset,
    output_file: Union[str, Path],
    delete_existing: bool = False,
    prompt_user: bool = True,
) -> bool:
    """Save an xarray Dataset to a NetCDF file.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset to save.
    output_file : str or Path
        Path where the file will be saved.
    delete_existing : bool, optional
        If True, overwrite existing files without prompting.
    prompt_user : bool, optional
        If True, ask user before overwriting existing files.

    Returns
    -------
    bool
        True if file was saved successfully, False if skipped.

    """
    output_path = Path(output_file)
    
    if output_path.exists():
        if delete_existing:
            output_path.unlink()
        elif prompt_user:
            response = input(f"File {output_path} exists. Overwrite? (y/n): ")
            if response.lower() not in ['y', 'yes']:
                print("Save cancelled.")
                return False
            output_path.unlink()
        else:
            print(f"File {output_path} exists and delete_existing=False. Skipping save.")
            return False
    
    # Ensure directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save the dataset
    ds.to_netcdf(output_path)
    print(f"Dataset saved to {output_path}")
    return True