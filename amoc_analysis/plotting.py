from pathlib import Path
from typing import Any, Tuple, Union

import matplotlib.pyplot as plt
import xarray as xr
from pandas import DataFrame
from pandas.io.formats.style import Styler


def plot_monthly_transport(
    ds: xr.Dataset, var: str = "moc_mar_hc10"
) -> Tuple[Any, Any]:
    """Plot original and monthly averaged transport time series.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with a time dimension and a transport variable.
    var : str, optional
        Name of the variable to plot. Default is "moc_mar_hc10".

    Returns
    -------
    tuple
        Figure and axis objects from matplotlib.
    """
    here = Path(__file__).resolve().parent
    style_file = here / "amoc_analysis.mplstyle"
    if style_file.exists():
        plt.style.use(style_file)

    da = ds[var]
    ds_monthly = ds.resample(TIME="ME").mean()

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(ds.TIME, da, color="grey", alpha=0.5, linewidth=0.5, label="Original")
    ax.plot(
        ds_monthly.TIME,
        ds_monthly[var],
        color="red",
        linewidth=1.0,
        label="Monthly Avg",
    )
    ax.axhline(0, color="black", linestyle="--", linewidth=0.5)

    ax.set_title("RAPID 26°N - AMOC Transport")

    # Use variable attributes if present
    label = da.attrs.get("long_name", var)
    units = da.attrs.get("units", "")
    ax.set_ylabel(f"{label} [{units}]" if units else label)
    ax.set_xlabel("Time")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend()
    plt.tight_layout()

    return fig, ax


def show_variables(data: Union[str, xr.Dataset]) -> Styler:
    """
    Processes an xarray Dataset or a netCDF file, extracts variable information,
    and returns a styled DataFrame with details about the variables.

    Parameters
    ----------
    data : str or xr.Dataset
        The input data, either a file path to a netCDF file or an xarray Dataset.

    Returns
    -------
    pandas.io.formats.style.Styler
        A styled DataFrame containing variable information including dimensions,
        names, units, and comments.
    """
    if isinstance(data, str):
        print(f"Information is based on file: {data}")
        dataset = xr.open_dataset(data)
        variables = dataset.variables
    elif isinstance(data, xr.Dataset):
        print("Information is based on xarray Dataset")
        variables = data.variables
    else:
        raise TypeError("Input data must be a file path (str) or an xarray Dataset")

    info = {}
    for i, key in enumerate(variables):
        var = variables[key]
        if isinstance(data, str):
            dims = var.dims[0] if len(var.dims) == 1 else "multiple"
            units = var.attrs.get("units", "")
            comment = var.attrs.get("comment", "")
        else:
            dims = var.dims[0] if len(var.dims) == 1 else "multiple"
            units = var.attrs.get("units", "")
            comment = var.attrs.get("comment", "")

        info[i] = {
            "name": key,
            "dims": dims,
            "units": units,
            "comment": comment,
            "standard_name": var.attrs.get("standard_name", ""),
            "dtype": str(var.dtype),
        }

    vars_df = DataFrame(info).T

    # Clean up dimensions display
    dims = vars_df.dims
    dims[dims.str.startswith("str")] = "string"
    vars_df["dims"] = dims

    vars_styled = (
        vars_df.sort_values(["dims", "name"])
        .reset_index(drop=True)
        .loc[:, ["dims", "name", "units", "comment", "standard_name", "dtype"]]
        .set_index("name")
        .style
    )

    return vars_styled


def show_attributes(data: Union[str, xr.Dataset]) -> DataFrame:
    """
    Processes an xarray Dataset or a netCDF file, extracts attribute information,
    and returns a DataFrame with details about the attributes.

    Parameters
    ----------
    data : str or xr.Dataset
        The input data, either a file path to a netCDF file or an xarray Dataset.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing attribute names, values, and data types.
    """
    if isinstance(data, str):
        print(f"Information is based on file: {data}")
        dataset = xr.open_dataset(data)
        attributes = dataset.attrs.keys()
        get_attr = lambda key: dataset.attrs[key]
    elif isinstance(data, xr.Dataset):
        print("Information is based on xarray Dataset")
        attributes = data.attrs.keys()
        get_attr = lambda key: data.attrs[key]
    else:
        raise TypeError("Input data must be a file path (str) or an xarray Dataset")

    info = {}
    for i, key in enumerate(attributes):
        dtype = type(get_attr(key)).__name__
        info[i] = {"Attribute": key, "Value": get_attr(key), "DType": dtype}

    attrs_df = DataFrame(info).T

    return attrs_df


def plot_time_series(
    ds: xr.Dataset, 
    var: str, 
    title: str = None,
    ylabel: str = None,
    color: str = "blue",
    figsize: Tuple[int, int] = (12, 6)
) -> Tuple[Any, Any]:
    """Plot a simple time series from an xarray Dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with a time dimension and the variable to plot.
    var : str
        Name of the variable to plot.
    title : str, optional
        Plot title. If None, uses variable's long_name or the variable name.
    ylabel : str, optional
        Y-axis label. If None, uses variable's long_name and units.
    color : str, optional
        Line color. Default is "blue".
    figsize : tuple, optional
        Figure size as (width, height). Default is (12, 6).

    Returns
    -------
    tuple
        Figure and axis objects from matplotlib.
    """
    da = ds[var]
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(ds.TIME, da, color=color, linewidth=1.0)
    
    # Set title
    if title is None:
        title = da.attrs.get("long_name", var)
    ax.set_title(title)
    
    # Set ylabel
    if ylabel is None:
        label = da.attrs.get("long_name", var)
        units = da.attrs.get("units", "")
        ylabel = f"{label} [{units}]" if units else label
    ax.set_ylabel(ylabel)
    
    ax.set_xlabel("Time")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    return fig, ax