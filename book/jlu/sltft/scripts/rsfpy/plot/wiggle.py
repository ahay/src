import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Union
import warnings

def wiggle(
    data: Union[np.ndarray, "Rsfarray"],
    *,
    d1: Optional[float] = None, o1: Optional[float] = None,
    d2: Optional[float] = None, o2: Optional[float] = None,
    min1: Optional[float] = None, max1: Optional[float] = None,
    min2: Optional[float] = None, max2: Optional[float] = None,
    label1: Optional[str] = None, label2: Optional[str] = None,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    ncolor: str = 'blue',
    pcolor: str = 'red',
    lcolor: str = 'black',
    zplot: float = 1.0,
    transp: Optional[bool] = False,
    yreverse: Optional[bool] = False,
    xreverse: Optional[bool] = False,
    **plot_params: Optional[dict]
) -> plt.Axes:
    """
    Plot 2D wiggle image.

    Parameters
    ----------
    data : Union[np.ndarray, "Rsfarray"]
        Input data array to be plotted.
    d1, o1, d2, o2 : Optional[float]
        Sampling intervals and origins for the axes.
    min1, max1, min2, max2 : Optional[float]
        Clipping limits for the data.
    label1, label2 : Optional[str]
        Labels for the axes.
    title : Optional[str]
        Title of the plot.
    ax : Optional[plt.Axes]
        Matplotlib axes to plot on.
    ncolor, pcolor, lcolor : str
        Colors for the negative, positive, and line plots.
    zplot : float
        Vertical scaling factor.
    transp : Optional[bool]
        Whether to apply transparency to the plot.
    yreverse, xreverse : Optional[bool]
        Whether to reverse the y or x axis.
    **plot_params : Optional[dict]
        Additional parameters for the plot.

    Returns
    -------
    plt.Axes
        The axes object with the wiggle plot.
    """
    # Check dimensions
    data = np.squeeze(data)
    if data.ndim < 2:
        raise ValueError("Input data must be at least 2D.")
    elif data.ndim > 2:
        warnings.warn("Got data dimensions > 2, use first slice.")
        data = data.reshape(data.shape[0], data.shape[1], -1)[:, :, 0]

    defaults = {
        'clip': None,
        'bias': 0,
        'allpos': False,
        'pclip': 99,
        'linewidth': 0.5,
        'interpolate': False
    }
    params = {**defaults, **plot_params}

    # Axes 创建
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    data = data if transp else data.T
    # 数据属性（Rsfarray 支持）
    if hasattr(data, "d1"):
        d1 = d1 if d1 is not None else getattr(data, "d1", None)
        d2 = d2 if d2 is not None else getattr(data, "d2", None)
        o1 = o1 if o1 is not None else getattr(data, "o1", None)
        o2 = o2 if o2 is not None else getattr(data, "o2", None)
        if label1 is None and hasattr(data, "label_unit"):
            label1 = data.label_unit(axis=1)
        if label2 is None and hasattr(data, "label_unit"):
            label2 = data.label_unit(axis=0)
    if d1 is None or d1 == 0: d1 = 1.
    if d2 is None or d2 == 0: d2 = 1.
    if o1 is None: o1 = 0.
    if o2 is None: o2 = 0.

    n1, n2 = data.shape

    if hasattr(data, "axis"):
        axis1, axis2 = data.axis([0, 1])
    else:
        axis1 = np.arange(n1) * d1 + o1
        axis2 = np.arange(n2) * d2 + o2



    # 坐标范围
    if min1 is None: min1 = o1
    if max1 is None: max1 = o1 + d1 * (n1-1)
    if min2 is None: min2 = o2
    if max2 is None: max2 = o2 + d2 * (n2-1)

    t = axis1 
    x_positions = axis2 # 道位置

    # ---- 振幅裁剪 ----
    clip = params['clip']
    bias = params['bias']
    pclip = params['pclip']
    allpos = params['allpos']

    if clip is None:
        if pclip is not None:
            if hasattr(data, "pclip"):
                clip = data.pclip(pclip)
            else:
                clip = np.percentile(np.abs(data), pclip)
        else:
            clip = np.max(np.abs(data))

    clip = abs(clip)
    if allpos:
        vmin, vmax = 0, clip
    else:
        vmin, vmax = bias - clip, bias + clip

    # ---- zplot 缩放系数 ----
    if zplot == 0:
        zplot = 1.0
    # 默认最大幅度
    default_amp = (max2 - min2) / n2 / 2
    scale = default_amp * abs(zplot) / clip

    # ---- 绘制每道 ----
    for i in range(n2):
        if not x_positions[i] >= min2 and x_positions[i] <= max2:
            continue
        trace = data[:, i] - bias
        wiggle_x = x_positions[i] + trace * scale
        # 画线
        ax.plot(wiggle_x, t, color=lcolor, linewidth=params['linewidth'])

        # 填充正值区
        pos_mask = trace > 0
        ax.fill_betweenx(t, x_positions[i], wiggle_x, where=pos_mask, facecolor=pcolor, interpolate=params['interpolate'])

        # 填充负值区
        neg_mask = trace < 0
        ax.fill_betweenx(t, x_positions[i], wiggle_x, where=neg_mask, facecolor=ncolor, interpolate=params['interpolate'])

    max2 + trace

    if label1:
        ax.set_xlabel(label1)
    if label2:
        ax.set_ylabel(label2)
    if title:
        ax.set_title(title)
    if yreverse:
        ax.set_ylim(max1, min1)
    else:
        ax.set_ylim(min1, max1)
    if xreverse:
        ax.set_xlim(max2, min2)
    else:
        ax.set_xlim(min2, max2)

    if plot_params.get('show', True):
        plt.show()
    return ax
