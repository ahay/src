import numpy as np
import matplotlib.pyplot as plt
import warnings
import matplotlib.ticker as mticker
from typing import Optional, Union
from .grey import grey
from .wiggle import wiggle
import matplotlib.transforms as transforms
# from mpl_toolkits.axisartist.floating_axes import \
    # GridHelperCurveLinear
# from mpl_toolkits.axisartist.grid_finder import MaxNLocator
import mpl_toolkits.axisartist as artist
from matplotlib.ticker import MaxNLocator as MaxNLocator1

def grey3flat(
    data: Union[np.ndarray, "Rsfarray"],
    *,
    frame1: int = 0, frame2: int = 0, frame3: int = 0,
    point1: float = 0.75, point2: float = 0.5,
    d1: Optional[float] = None, o1: Optional[float] = None,
    d2: Optional[float] = None, o2: Optional[float] = None,
    d3: Optional[float] = None, o3: Optional[float] = None,
    label1: Optional[str] = None, label2: Optional[str] = None, label3: Optional[str] = None,
    title: Optional[str] = None,
    vmin: Optional[float] = None, vmax: Optional[float] = None,
    clip: Optional[float] = None, pclip: Optional[float] = 99,
    bias: Optional[float] = None, allpos: bool = False,
    cmap: str = "seismic",
    colorbar: bool = False,
    ax: Optional[plt.Axes] = None,
    cax: Optional[plt.Axes] = None,
    wig: bool = False,
    **plot_params
) -> plt.Figure:
    """
    三视图灰度/彩色绘制，框架与 grey 一致
    """
    data = np.squeeze(data)
    # ---- 1. 数据检查 ----
    if data.ndim < 3:
        raise ValueError(f"Input data must be at least 3D, got shape {data.shape}")
    elif data.ndim > 3:
        warnings.warn("Got data dimensions > 3, use first cube.")
        data = data.reshape(data.shape[0], data.shape[1], data.shape[2], -1)[:, :, :, 0]

    ny, nx, nz = data.shape
    frame1 = min(max(0, int(frame1)), ny - 1)
    frame2 = min(max(0, int(frame2)), nx - 1)
    frame3 = min(max(0, int(frame3)), nz - 1)

    # ---- 2. 从 Rsfarray 读取属性 ----
    if hasattr(data, "d1"):
        d1 = d1 if d1 is not None else getattr(data, "d1", 1.0)
        d2 = d2 if d2 is not None else getattr(data, "d2", 1.0)
        d3 = d3 if d3 is not None else getattr(data, "d3", 1.0)
        o1 = o1 if o1 is not None else getattr(data, "o1", 0.0)
        o2 = o2 if o2 is not None else getattr(data, "o2", 0.0)
        o3 = o3 if o3 is not None else getattr(data, "o3", 0.0)
        axis1, axis2, axis3 = data.axis([0,1,2])
        if label1 is None and hasattr(data, "label_unit"):
            label1 = data.label_unit(axis=0)
        if label2 is None and hasattr(data, "label_unit"):
            label2 = data.label_unit(axis=1)
        if label3 is None and hasattr(data, "label_unit"):
            label3 = data.label_unit(axis=2)
    else:
        d1 = 1.0 if d1 is None else d1
        d2 = 1.0 if d2 is None else d2
        d3 = 1.0 if d3 is None else d3
        o1 = 0.0 if o1 is None else o1
        o2 = 0.0 if o2 is None else o2
        o3 = 0.0 if o3 is None else o3
        axis1 = np.linspace(o1, o1 + d1 * ny, ny)
        axis2 = np.linspace(o2, o2 + d2 * nx, nx)
        axis3 = np.linspace(o3, o3 + d3 * nz, nz)

    # ---- 3. 数据切片 ----
    if hasattr(data, "window"):
        slice1 = data.window(n3=1, f3=frame3)
        slice2 = data.window(n2=1, f2=frame2)
        slice3 = data.window(n1=1, f1=frame1)
    else:
        slice1 = data[:, :, frame3]  # n1-n2
        slice2 = data[:, frame2, :]  # n1-n3
        slice3 = data[frame1, :, :]  # n2-n3

    # ---- 4. 颜色范围计算（与 grey 一致）----
    allslice = np.concatenate([slice1.ravel(), slice2.ravel(), slice3.ravel()])
    if bias is None:
        bias = 0
    if pclip is not None and clip is None and vmin is None and vmax is None:
        clip_val = np.percentile(np.abs(allslice), pclip)
        if allpos:
            vmin, vmax = 0, clip_val
        else:
            vmin, vmax = bias - clip_val, bias + clip_val
    if clip is not None:
        clip = abs(clip)
        if allpos:
            vmin, vmax = 0, clip
        else:
            vmin, vmax = bias - clip, bias + clip
    elif vmin is not None or vmax is not None:
        if vmin is None: vmin = -vmax
        if vmax is None: vmax = -vmin
        if vmin > vmax: vmin, vmax = vmax, vmin

    # ---- 5. 创建 Figure 和布局 ----
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    ax_main = ax
    ax_main.axis("off")

    if title: 
        height = 0.9
    else:
        height = 1

    if colorbar and cax is None and not wig:
        # Stole colorbar axes
        width = 0.9
    else: width = 1.

    # 三个 inset_axes
    ax1 = ax_main.inset_axes([0.0, 0.0, point2 * width, point1 * height])
    ax2 = ax_main.inset_axes([point2*width, 0.0, (1 - point2) * width, point1 * height], sharey=ax1)
    ax3 = ax_main.inset_axes([0.0, point1 * height, point2 * width, (1 - point1) * height], sharex=ax1)

    # # ---- 6. extent 计算 ----
    # extent1 = (o2, o2 + d2 * nx, o1 + d1 * ny, o1)  # slice1
    # extent2 = (o3, o3 + d3 * nz, o1 + d1 * ny, o1)  # slice2
    # extent3 = (o2, o2 + d2 * nx, o3, o3 + d3 * nz)  # slice3

    # ---- 7. 绘制 ----
    if wig:
        wigparas = {
            'zplot': plot_params.get('zplot', 1),
            'ncolor': plot_params.get('ncolor', 'blue'),
            'pcolor': plot_params.get('pcolor', 'red'),
            'lcolor': plot_params.get('lcolor', 'black'),
        }
        wiggle(slice1, ax=ax1, clip=clip, show=False, transp=True, yreverse=True, **wigparas)
        wiggle(slice2, ax=ax2, clip=clip, show=False, transp=True, yreverse=True, **wigparas)
        wiggle(slice3, ax=ax3, clip=clip, show=False, transp=False, yreverse=False, **wigparas)
        colorbar=False
    else:
        grey(slice1, ax=ax1, cmap=cmap, vmin=vmin, vmax=vmax, show=False)
        grey(slice2, ax=ax2, cmap=cmap, vmin=vmin, vmax=vmax, show=False)
        grey(slice3, ax=ax3, cmap=cmap, vmin=vmin, vmax=vmax, show=False, transp=False, yreverse=False)

        im1 = ax1.images[0]
        im2 = ax2.images[0]
        im3 = ax3.images[0]



    ax2.yaxis.set_visible(False)
    ax3.xaxis.set_visible(False)

    lcol = plot_params.get('framelinecol', 'blue' if cmap == 'grey' else 'black')
    ax1.hlines(y=axis1[frame1], xmin=axis2[0], xmax=axis2[-1], color=lcol)
    ax1.vlines(x=axis2[frame2], ymin=axis1[0], ymax=axis1[-1], color=lcol)
    ax2.hlines(y=axis1[frame1], xmin=axis3[0], xmax=axis3[-1], color=lcol)
    ax2.vlines(x=axis3[frame3], ymin=axis1[0], ymax=axis1[-1], color=lcol)
    ax3.hlines(y=axis3[frame3], xmin=axis2[0], xmax=axis2[-1], color=lcol)
    ax3.vlines(x=axis2[frame2], ymin=axis3[0], ymax=axis3[-1], color=lcol)

    ax2.text(x=axis3[-1]+d3, y=axis1[frame1], s=f"{axis1[frame1]:.3g}", ha="left", va="center", color=lcol, rotation=90)
    ax2.text(x=axis3[frame3], y=axis1[0]-d1, s=f"{axis3[frame3]:.3g}", ha="center", va="bottom", color=lcol, rotation=0)
    ax3.text(x=axis2[frame2], y=axis3[-1]+d3, s=f"{axis2[frame2]:.3g}", ha="center", va="bottom", color=lcol, rotation=0)

    ticks1 = ax1.get_yticks()
    ticks2 = ax3.get_yticks()
    all_ticks = np.concatenate([ticks1, ticks2])
    diffs = np.diff(np.unique(np.round(all_ticks, 10)))
    min_step = np.min(np.abs(diffs[diffs > 0]))
    decimals = max(0, -int(np.floor(np.log10(min_step))))
    fmt_str = f"%.{decimals}f"
    ax3.yaxis.set_major_formatter(mticker.FormatStrFormatter(fmt_str))
    ax1.yaxis.set_major_formatter(mticker.FormatStrFormatter(fmt_str))

    # ---- 8. colorbar ----
    if colorbar:
        if cax is not None:
            fig.colorbar(im1, cax=cax)
        else:
            cwidth = (1 - width)/3
            cax = ax_main.inset_axes([width+cwidth, 0.0, cwidth, height])
            fig.colorbar(im1, cax=cax)

    if title:
        ax_title = ax_main.inset_axes([0.0, height, 1, 1 - height])
        ax_title.text(0.5, 0.5, title, ha="center", va="center", fontsize=12)
        ax_title.axis("off")
    if plot_params.get('show', True):
        plt.show()

    return fig


def grey3cube(
    data: Union[np.ndarray, "Rsfarray"],
    *,
    frame1: int = 0, frame2: int = 0, frame3: int = 0,
    point1: float = 0.75, point2: float = 0.5,
    d1: Optional[float] = None, o1: Optional[float] = None,
    d2: Optional[float] = None, o2: Optional[float] = None,
    d3: Optional[float] = None, o3: Optional[float] = None,
    label1: Optional[str] = None, label2: Optional[str] = None, label3: Optional[str] = None,
    title: Optional[str] = None,
    vmin: Optional[float] = None, vmax: Optional[float] = None,
    clip: Optional[float] = None, pclip: Optional[float] = 99,
    bias: Optional[float] = None, allpos: bool = False,
    cmap: str = "seismic",
    colorbar: bool = False,
    ax: Optional[plt.Axes] = None,
    cax: Optional[plt.Axes] = None,
    **plot_params
) -> plt.Figure:
   
     # ---- 1. 数据检查 ----
    data = np.squeeze(data)
    if data.ndim < 3:
        raise ValueError(f"Input data must be at least 3D, got shape {data.shape}")
    elif data.ndim > 3:
        warnings.warn("Got data dimensions > 3, use first cube.")
        data = data.reshape(data.shape[0], data.shape[1], data.shape[2], -1)[:, :, :, 0]


    ny, nx, nz = data.shape
    frame1 = min(ny-1, frame1)
    frame2 = min(nx-1, frame2)
    frame3 = min(nz-1, frame3)
    if hasattr(data, "d1"):
        d1 = d1 if d1 is not None else getattr(data, "d1", 1.0)
        d2 = d2 if d2 is not None else getattr(data, "d2", 1.0)
        d3 = d3 if d3 is not None else getattr(data, "d3", 1.0)
        o1 = o1 if o1 is not None else getattr(data, "o1", 0.0)
        o2 = o2 if o2 is not None else getattr(data, "o2", 0.0)
        o3 = o3 if o3 is not None else getattr(data, "o3", 0.0)
        axis1, axis2, axis3 = data.axis([0,1,2])
        if label1 is None and hasattr(data, "label_unit"):
            label1 = data.label_unit(axis=0)
        if label2 is None and hasattr(data, "label_unit"):
            label2 = data.label_unit(axis=1)
        if label3 is None and hasattr(data, "label_unit"):
            label3 = data.label_unit(axis=2)
    else:
        d1 = 1.0 if d1 is None else d1
        d2 = 1.0 if d2 is None else d2
        d3 = 1.0 if d3 is None else d3
        o1 = 0.0 if o1 is None else o1
        o2 = 0.0 if o2 is None else o2
        o3 = 0.0 if o3 is None else o3
        axis1 = np.linspace(o1, o1 + d1 * ny, ny)
        axis2 = np.linspace(o2, o2 + d2 * nx, nx)
        axis3 = np.linspace(o3, o3 + d3 * nz, nz)
    # ---- 3. 数据切片 ----
    if hasattr(data, "window"):
        slice1 = data.window(n3=1, f3=frame3)
        slice2 = data.window(n2=1, f2=frame2)
        slice3 = data.window(n1=1, f1=frame1)
    else:
        slice1 = data[:, :, frame3]  # n1-n2
        slice2 = data[:, frame2, :]  # n1-n3
        slice3 = data[frame1, :, :]  # n2-n3
    allslice = np.concatenate([slice1.ravel(), slice2.ravel(), slice3.ravel()])
    if bias is None:
        bias = 0
    if pclip is not None and clip is None and vmin is None and vmax is None:
        clip_val = np.percentile(np.abs(allslice), pclip)
        if allpos:
            vmin, vmax = 0, clip_val
        else:
            vmin, vmax = bias - clip_val, bias + clip_val
    if clip is not None:
        clip = abs(clip)
        if allpos:
            vmin, vmax = 0, clip
        else:
            vmin, vmax = bias - clip, bias + clip
    elif vmin is not None or vmax is not None:
        if vmin is None: vmin = -vmax
        if vmax is None: vmax = -vmin
        if vmin > vmax: vmin, vmax = vmax, vmin



    n1tic, n2tic, n3tic = plot_params.get("n1tic", 5), plot_params.get("n2tic", 5), plot_params.get("n3tic", 5)

    
    amax1, amax2, amax3 = axis1[-1], axis2[-1], axis3[-1]
    len1, len2, len3 = amax1 - o1, amax2 - o2, amax3 - o3
    if len1==0:
        len1 = 1
        o1 = 0
        amax1 = 1
    if len2==0:
        len2 = 1
        o2 = 0
        amax2 = 1
    if len3==0:
        len3 = 1
        o3 = 0
        amax3 = 1
    extents = [
        [o2, amax2, amax1, o1],
        [0, len2, len3, 0],
        [0, len3, 0, len1],
    ]

    axis_xlims = [
        [o2, amax2],
        [o2, amax2 + len2 / point2 * (1 - point2)],
        [o3, amax3]
    ]
    axis_ylims = [
        [amax1, o1],
        [o3, amax3],
        [o1, amax1 + len1 / point1 * (1 - point1)]
    ]


    axshear = (1 - point2) / len3 * len2 / point2
    axshift = -(1 - point2) / len3 * len2 / point2 * o3
    dxshear = (1 - point2) / len3
    dxshift = 0
    ayshear = (1 - point1) / len3 * len1 / point1
    ayshift = -(1 - point1) / len3 * len1 / point1 * o3
    dyshear = (1 - point1) / len3
    dyshift = 0

    l11 = (amax1 - axis1[frame1]) / len1 * point1
    l12 = (axis2[frame2] - o2) / len2 * point2

    loff1 = (axis3[frame3] - o3) / len3 * (1 - point1)
    loff2 = (axis3[frame3] - o3) / len3 * (1 - point2)

    dtrans1 = transforms.Affine2D(
        np.array([[point2 / len2, dxshear, 0],
         [0, 1 / len3, 0],
         [0, 0, 1]])
    )
    dtrans2 = transforms.Affine2D(
        np.array([[1 / len3, 0, 0],
         [dyshear, point1 / len1, 0],
         [0, 0, 1]])
    )
    atrans1 = transforms.Affine2D(
        np.array([[1, axshear, axshift],
         [0, 1, 0],
         [0, 0, 1]])
    )
    atrans2 = transforms.Affine2D(
        np.array([[1, 0, 0],
         [ayshear, 1, ayshift],
         [0, 0, 1]])
    )
    grid_helper = artist.floating_axes.GridHelperCurveLinear(atrans1,
                                        extremes=[o2, amax2, o3, amax3],
                                        grid_locator1=artist.grid_finder.MaxNLocator(nbins=n2tic),
                                        grid_locator2=artist.grid_finder.MaxNLocator(nbins=n3tic))
    grid_helper1 = artist.floating_axes.GridHelperCurveLinear(atrans2,
                                         extremes=[o3, amax3, o1, amax1],
                                         grid_locator1=artist.grid_finder.MaxNLocator(nbins=n1tic),
                                         grid_locator2=artist.grid_finder.MaxNLocator(nbins=n1tic)
                                         )

    topmargin = 0.05


    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    axbase0 = ax

    axbase0.axis('off')
    axbase0.set_label("Base axes")

    if title:
        title_pos = [0, 1 - topmargin, 1, topmargin]
        ax_title = axbase0.inset_axes(title_pos)
        ax_title.set_frame_on(False)
        ax_title.get_xaxis().set_visible(False)
        ax_title.get_yaxis().set_visible(False)
        ax_title.set_xlim([0, 1])
        ax_title.set_ylim([0, 1])
        ax_title.text(0.5, 1, title, va='bottom', ha='center')
        axbase0.add_child_axes(ax_title)
    else:
        topmargin = 0.
    hei_ax = 1 - topmargin
    

    if colorbar and cax is None:
        axbase = axbase0.inset_axes(bounds=[0, 0, 0.85, 1 - topmargin], facecolor="none")
        axbar = axbase0.inset_axes(bounds=[0.90, 0, 0.05, 1 - topmargin], facecolor="none")
        axbase.set_label("Base axes with bar")
    else:
        axbase = axbase0
    axbase.axis('off')

    ax1 = axbase.inset_axes(bounds=[0, 0, point2, point1*hei_ax], facecolor="none")
    ax1.set_label("Axes 1 [n1, n2]")

    im0 = ax1.imshow(slice1, aspect="auto", extent=extents[0],
                     cmap=cmap, vmin=vmin, vmax=vmax)

    ax1.yaxis.set_major_locator(MaxNLocator1(nbins=n1tic))
    ax1.set_xlim(axis_xlims[0])
    ax1.set_ylim(axis_ylims[0])

    ax2 = axbase.inset_axes(bounds=[0, point1*hei_ax, 1, (1 - point1)*hei_ax],
                            axes_class=artist.Axes, grid_helper=grid_helper,
                            facecolor="none")
    ax2.set_label("Axes 3 [n1, n3]")

    im1 = ax2.imshow(slice3.T, aspect="auto", extent=extents[1],
                     cmap=cmap, vmin=vmin, vmax=vmax)
 
    im1.set_transform(dtrans1 + ax2.transAxes)

    #
    ax2.set_ylim(axis_ylims[1])
    ax2.set_xlim(axis_xlims[1])
    # ax2.axis["t1"] = ax2.new_floating_axis(0,-5)

    ax3 = axbase.inset_axes(bounds=[point2, 0, 1 - point2, hei_ax],
                            axes_class=artist.Axes, grid_helper=grid_helper1,
                            facecolor="none")
    ax3.set_label("Axes 2 [n2, n3]")

    im2 = ax3.imshow(slice2, aspect="auto", extent=extents[2],
                     cmap=cmap, vmin=vmin, vmax=vmax)
    im2.set_transform(dtrans2 + ax3.transAxes)
    #
    ax3.set_xlim(axis_xlims[2])
    ax3.set_ylim(axis_ylims[2])

    ax2.axis[:].major_ticklabels.set(visible=False)
    ax2.axis[:].major_ticks.set(visible=False)

  

    
    ax2.axis["left"].major_ticks.set(tick_out=True, visible=True)
    ax2.axis["right"].set(visible=False)
    ax2.axis["left"].major_ticklabels.set(visible=True)



    ax3.axis[:].major_ticklabels.set(visible=False)
    ax3.axis[:].major_ticks.set(visible=False)
    ax3.axis["right"].major_ticks.set(visible=False)

    # Labels
    ax1.set_xlabel(label2)
    ax1.set_ylabel(label1)

    ax2.axis["left"].label.set(text=label3)


    # Colorbar
    if colorbar and cax is None:
        # cbar = fig.colorbar(im0, cax=axbar, label=bartitle, orientation=barpos)
        cbar = fig.colorbar(im0, cax=axbar, cmap=cmap, label=plot_params.get('bartitle', ''),)

    # Indicating lines

    l21 = point1 + loff1
    l22 = point2 + loff2
    l221 = l12 +  (1-point2)

    l31 = l11 + 1 - point1
    l32 = point1 + loff1

    lcol = plot_params.get('framelinecol', 'blue' if cmap == 'grey' else 'black')

    ax1.hlines(l11*hei_ax,0,point2, color=lcol,transform=axbase.transAxes)
    ax1.vlines(l12,0,point1*hei_ax, color=lcol, transform=axbase.transAxes)
    ax2.plot([loff2,l22],[l21*hei_ax,l21*hei_ax], color=lcol, transform=axbase.transAxes)
    ax2.plot([l12,l221],[point1*hei_ax,1*hei_ax], color=lcol, transform=axbase.transAxes)
    ax3.vlines(1, (1-point1)*hei_ax, 1*hei_ax, color=lcol, transform=axbase.transAxes)
    ax3.vlines(l22, loff1*hei_ax, l32*hei_ax, color=lcol,  transform=axbase.transAxes)
    ax3.plot([point2,1],[l11*hei_ax,l31*hei_ax], color=lcol, transform=axbase.transAxes)

    # Indicating labels
    lab1 = str(axis2[frame2])
    lab11 = "%.2f" % axis2[frame2]
    lab1 = lab11 if len(lab11) < len(lab1) else lab1
    lab2 = str(axis1[frame1])
    lab21 = "%.2f" % axis1[frame1]
    lab2 = lab21 if len(lab21) < len(lab2) else lab2
    lab3 = str(axis3[frame3])
    lab31 = "%.2f" % axis3[frame3]
    lab3 = lab31 if len(lab31) < len(lab3) else lab3
    axbase.text(l221, 1*hei_ax, "%s" % lab1, va='bottom',
                ha='center', color=lcol)
    axbase.text(1, l31*hei_ax, "%s" % lab2, va='center',
                ha='left', color=lcol, rotation=-90)
    axbase.text(l22, loff1*hei_ax,
                "%s" % lab3, va='top', ha='left', color=lcol)


def grey3(*args, **kargs):
    if kargs.get('flat', False):
        return grey3flat(*args, **kargs)
    else:
        return grey3cube(*args, **kargs)