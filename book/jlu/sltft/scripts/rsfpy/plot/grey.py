import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Union
import warnings

def grey(
    data: Union[np.ndarray, "Rsfarray"],
    *,
    d1: Optional[float] = None, o1: Optional[float] = None,
    d2: Optional[float] = None, o2: Optional[float] = None,
    min1: Optional[float] = None, max1: Optional[float] = None,
    min2: Optional[float] = None, max2: Optional[float] = None,
    label1: Optional[str] = None, label2: Optional[str] = None,
    title: Optional[str] = None,
    ax: Optional[plt.Axes] = None,
    cax: Optional[plt.Axes] = None, 
    colorbar: Optional[bool] = False,
    transp: Optional[bool] = True,
    yreverse: Optional[bool] = True,
    xreverse: Optional[bool] = False,
    **plot_params: Optional[dict]
) -> plt.Axes:
    """
    Plot 2D image using imshow

    Parameters
    ----------
    data : Union[np.ndarray, "Rsfarray"]
        The input data to plot.
    d1, o1, d2, o2 : Optional[float]
        The spatial sampling and origin for the axes.
    min1, max1, min2, max2 : Optional[float]
        The minimum and maximum values for the axes.
    label1, label2 : Optional[str]
        The labels for the axes.
    title: Optional[str]
        The title of the plot.
    clip : Optional[float]
        The clipping value for the data.
    pclip: Optional[float]
        The percentile clipping value for the data.
    bias: Optional[float]
        The bias value for the data.
    allpos: Optional[bool]
        Whether to use all positive values for the color scale.
    ax : Optional[plt.Axes]
        The axes to plot on. If None, a new figure and axes will be created.
    colorbar: Optional[bool]
        Whether to display a colorbar.
    transp: Optional[bool]
        Transpose data.
    yreverse: Optional[bool]
        Reverse Y Axis (default is True)
    xreverse: Optional[bool]
        Reverse X Axis (default is False)
    cax : Optional[plt.Axes]
        The colorbar axes. If None, a new colorbar will be created.

    Returns
    -------
    plt.Axes
        The axes object with the image plot.
    """
    # Check dimensions
    data = np.squeeze(data)
    if data.ndim < 2:
        raise ValueError("Input data must be at least 2D.")
    elif data.ndim > 2:
        warnings.warn("Got data dimensions > 2, use first slice.")
        data = data.reshape(data.shape[0], data.shape[1], -1)[:, :, 0]

    defaults = {
        'cmap': 'gray',
        'aspect': 'auto',
        'origin': 'upper',
        'interpolation': 'nearest',
        'vmin': None,
        'vmax': None,
        'clip': None,
        'bias': None,
        'allpos': False,
        'pclip': 99
    }
    params = {**defaults, **plot_params}

    # ---- 1. Axes 创建 ----
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    data = data if transp else data.T
    # ---- 2. 数据属性获取 ----
    if hasattr(data, "d1"):
        d1 = d1 if d1 is not None else getattr(data, "d1", None)
        d2 = d2 if d2 is not None else getattr(data, "d2", None)
        o1 = o1 if o1 is not None else getattr(data, "o1", None)
        o2 = o2 if o2 is not None else getattr(data, "o2", None)
        if label1 is None and hasattr(data, "label_unit"):
            label1 = data.label_unit(axis=0)
        if label2 is None and hasattr(data, "label_unit"):
            label2 = data.label_unit(axis=1)
    if d1 is None or d1 == 0: d1 = 1.
    if d2 is None or d2 == 0: d2 = 1.
    if o1 is None: o1 = 0.
    if o2 is None: o2 = 0.

    # ---- 3. extent 计算 ----
    ny, nx = data.shape
    if yreverse:
        params['origin'] = 'upper'
        extent1 = (o1 + d1 * (ny-1), o1)
    else:
        params['origin'] = 'lower'
        extent1 = (o1, o1 + d1 * (ny-1))
    if xreverse:
        extent2 = (o2 + d2 * (nx-1), o2)
    else:
        extent2 = (o2, o2 + d2 * (nx-1))
    extent = extent2 + extent1
    if min1 is None: min1 = o1
    if max1 is None: max1 = o1 + d1 * (ny-1)
    if min2 is None: min2 = o2
    if max2 is None: max2 = o2 + d2 * (nx-1)

    # ---- 4. 计算 vmin / vmax ----
    vmin = params['vmin']
    vmax = params['vmax']
    clip = params['clip']
    bias = params['bias']
    allpos = params['allpos']
    pclip = params['pclip']

    if bias is None:
        bias = 0

    # 优先级 1: pclip
    if pclip is not None and (clip is None) and (vmin is None and vmax is None):
        if hasattr(data, "pclip"):
            # Rsfarray 自带方法
            clip = data.pclip(pclip)
        else:
            clip = np.percentile(np.abs(data), pclip)
        if allpos:
            vmin, vmax = 0, clip
        else:
            vmin, vmax = bias - clip, bias + clip


    # 优先级 2: clip
    if clip is not None:
        clip = abs(clip)
        if allpos:
            vmin, vmax = 0, clip
        else:
            vmin, vmax = bias - clip, bias + clip

    # 优先级 3: vmax/vmin 输入，但前面 clip/pclip 都没走时
    elif vmin is not None or vmax is not None:
        if vmin is None and vmax is not None:
            vmin = -vmax
        elif vmax is None and vmin is not None:
            vmax = -vmin
        if vmin > vmax:
            vmin, vmax = vmax, vmin
        # clip/bias 未触发，不调整 bias

    # 最终赋值
    imshow_kwargs = {k: v for k, v in params.items()
                     if k in ['cmap', 'origin', 'interpolation','aspect']}
    imshow_kwargs.update({'vmin': vmin, 'vmax': vmax})

    # ---- 5. 绘制 ----
    im = ax.imshow(data, extent=extent, **imshow_kwargs)

    # ---- 6. colorbar ----
    if colorbar:
        if cax is not None:
            fig.colorbar(im, cax=cax)
        else:
            fig.colorbar(im, ax=ax)

    # ---- 7. 坐标轴标签 ----
    if label1:
        ax.set_ylabel(label1)
    if label2:
        ax.set_xlabel(label2)
    if title:
        ax.set_title(title)

    ax.set_xlim(min2, max2)
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
