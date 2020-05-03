import sys
import os
import neo
import quantities as pq
import elephant
from scipy.interpolate import interp1d, interp2d
from scipy.ndimage import gaussian_filter
from scipy.stats import zscore
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse


def stretch_to_framerate(t_start, t_stop, num_frames, frame_rate=None):
    if args.frame_rate is None:
        return np.arange(num_frames, dtype=int)
    else:
        new_num_frames = (t_stop.rescale('s').magnitude
                        - t_start.rescale('s').magnitude) \
                        * args.frame_rate
        return np.linspace(0, num_frames-1, new_num_frames).astype(int)


def AnalogSignal2ImageSequence(block):
    for seg_count, segment in enumerate(block.segments):
        block.segments[seg_count].imagesequences = []
        for asig_count, asig in enumerate(segment.analogsignals):
            asig_array = asig.as_array()
            dim_t, dim_channels = asig_array.shape

            if 'x_coords' not in asig.array_annotations\
                or 'y_coords' not in asig.array_annotations:
                print('AnalogSignal {} in Segment {} has no spatial Information '\
                      .format(asig_count, seg_count)\
                    + ' as array_annotations "x_coords" "y_coords", skip.')
                break

            coords = np.array([(x,y) for x,y in zip(asig.array_annotations['x_coords'],
                                                    asig.array_annotations['y_coords'])],
                              dtype=float)

            if len(coords) != dim_channels:
                raise IndexError("Number of channels doesn't agree with "\
                               + "number of coordinates!")

            dim_x = np.max(asig.array_annotations['x_coords']) + 1
            dim_y = np.max(asig.array_annotations['y_coords']) + 1

            image_data = np.empty((dim_t, dim_x, dim_y), dtype=asig.dtype)
            image_data[:] = np.nan

            for channel in range(dim_channels):
                x, y = coords[channel]
                x, y = int(x), int(y)
                image_data[:, x, y] = asig_array[:, channel]

            spatial_scale = asig.annotations['spatial_scale']

            # array_annotations = {}
            # for k, v in asig.array_annotations.items():
            #     array_annotations[k] = v.reshape((dim_x, dim_y))

            imgseq = neo.ImageSequence(image_data=image_data,
                                       units=asig.units,
                                       dtype=asig.dtype,
                                       # t_start=asig.t_start, # NotImplementedError
                                       sampling_rate=asig.sampling_rate,
                                       name=asig.name,
                                       description=asig.description,
                                       file_origin=asig.file_origin,
                                       # array_annotations=array_annotations,
                                       **asig.annotations)

            block.segments[seg_count].imagesequences.append(imgseq)
    return block


def plot_frame(frame, cmap=None, vmin=None, vmax=None,
               colorbar=True, ax=None, interpolate='cubic'):
    if ax is None:
        fig, ax = plt.subplots()
    if cmap is None:
        cmap = plot.cm.gray
    else:
        cmap = plt.get_cmap(cmap)

    if interpolate is not None:
        dim_x, dim_y = frame.shape
        f = interp2d(np.arange(dim_x), np.arange(dim_y), np.nan_to_num(frame), kind=interpolate)
        frame = f(np.arange(0, dim_x-1, .1), np.arange(0, dim_y-1, .1))

    img = ax.imshow(frame, interpolation='nearest',
                    cmap=cmap, vmin=vmin, vmax=vmax)
    if colorbar:
        plt.colorbar(img, ax=ax)
    ax.axis('image')
    ax.set_xticks([])
    ax.set_yticks([])
    return img


def smooth_frames(frames, sigma):
    # replace nan sites by median
    if np.isfinite(frames).any():
        # assume constant nan sites over time
        nan_sites = np.where(np.bitwise_not(np.isfinite(frames[0])))
        if np.iscomplexobj(frames):
            frames[:,nan_sites[0],nan_sites[1]] = np.nanmedian(np.real(frames)) \
                                                + np.nanmedian(np.imag(frames))*1j
        else:
            frames[:,nan_sites[0],nan_sites[1]] = np.nanmedian(frames)
    else:
        nan_sites = None

    # apply gaussian filter
    if np.iscomplexobj(frames):
        frames = gaussian_filter(np.real(frames), sigma=sigma, mode='nearest') \
               + gaussian_filter(np.imag(frames), sigma=sigma, mode='nearest')*1j
    else:
        frames = gaussian_filter(frames, sigma=sigma, mode='nearest')

    # set nan sites back to nan
    if nan_sites is not None:
        if np.iscomplexobj(frames):
            frames[:,nan_sites[0],nan_sites[1]] = np.nan + np.nan*1j
        else:
            frames[:,nan_sites[0],nan_sites[1]] = np.nan

    return frames


def plot_figure(frame_id, segment, t_window, cortex_img):
    fig = plt.figure(figsize=(11,9))
    gs = mpl.gridspec.GridSpec(nrows=3, ncols=5, width_ratios=[1,1,1,1,.1], height_ratios=[1,.85,1])
    gs.update(hspace=0.25)
    ax = [fig.add_subplot(gs[0, 0:]),
          fig.add_subplot(gs[2, 0:]),
          fig.add_subplot(gs[1, 0]),
          fig.add_subplot(gs[1, 1]),
          fig.add_subplot(gs[1, 2]),
          fig.add_subplot(gs[1, 3]),
          fig.add_subplot(gs[1, 4]),
         ]

    colors = ['k',
              '#4c72b0', # Delta
              '#c44e52', # Theta
              '#dd8452', # Alpha
              '#55a868'  # Beta
             ]
    alpha = .7

    raw_asig = segment.analogsignals[0]
    t = raw_asig.times[frame_id].rescale('s')
    tw_start = max([raw_asig.t_start.rescale('s'), t-t_window/2])
    tw_stop = tw_start + t_window
    if tw_stop > raw_asig.t_stop:
        tw_stop = raw_asig.t_stop
        tw_start = tw_stop - t_window

    # Frames
    for i, imgseq in enumerate(segment.imagesequences[1:]):
        frames = smooth_frames(imgseq.as_array(), .6)
        frame = np.angle(frames[frame_id])
        img = plot_frame(frame,
                         cmap='twilight',
                         ax=ax[i+2],
                         colorbar=False,
                         vmin=-np.pi,
                         vmax=np.pi,
                         interpolate='linear'
                        )
        ax[i+2].set_title(imgseq.name, weight='bold')
    ax[2].set_ylabel(r'$\longleftarrow\; 4mm\; \longrightarrow$')
    ax[2].set_xlabel(r'$\longleftarrow\; 4mm\; \longrightarrow$')

    # Colorbar
    ax[-1].set_axis_off()
    cax = ax[-1].inset_axes([0,0,.8,1])
    cbar = plt.colorbar(img, ax=ax[-1], cax=cax,
                        ticks=[-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
    cbar.set_label(label='Phase', weight='bold')
    cax.yaxis.labelpad = -10
    cbar.ax.set_yticklabels([r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'])

    # PSD
    xlim = (0.95, 35)
    ylim = (0.5 , 50000)
    ax[0].set_xlim(xlim)
    ax[0].set_ylim((ylim))
    ax[0].set_xlabel('Frequency [Hz]', weight='bold')
    ax[0].xaxis.set_label_coords(.94,-.04)
    ax[0].set_ylabel('Power Spectral Density', weight='bold')
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].yaxis.set_ticks_position('left')
    ax[0].xaxis.set_ticks_position('bottom')

    freqs, psd = elephant.spectral.welch_psd(raw_asig.time_slice(tw_start, tw_stop),
                                         freq_res=1/t_window,
                                         overlap=.7)

    xvalues = np.linspace(xlim[0], xlim[1], 200)
    avg_psd = np.mean(psd, axis=0)
    avg_psd_f = interp1d(freqs, avg_psd, kind='cubic')
    ax[0].loglog(xvalues, avg_psd_f(xvalues), linewidth=3,
                 color=colors[0], label='Average')
    for channel_psd in psd:
        f = interp1d(freqs, channel_psd, kind='slinear')
        ax[0].loglog(xvalues, f(xvalues), alpha=0.15,
                     color=colors[0], label='Channels')

    xfreqs = [1,5,10,15,20,25,30]
    ax[0].set_xticks(xfreqs)
    ax[0].set_xticklabels(xfreqs)
    ax[0].tick_params(axis="x",direction="in", pad=-15)
    handles, labels = ax[0].get_legend_handles_labels()

    ax[0].legend(handles[:2], labels[:2], frameon=False, fontsize=11.5,
                 loc=3, bbox_to_anchor=(.03, .05, 0.3, 0.3))

    # Freq areas
    trans = lambda a, x, y: fig.transFigure.inverted().transform(ax[a].transData.transform([x, y]))
    def add_polygon(ax2, x1, x3, color=None):
        ax1 = 0
        y1 = ylim[0]
        y2 = 0
        x2 = 1
        x4 = 89
        tx1, ty1 = trans(ax1, x1, y1)
        tx2, ty2 = trans(ax2, x2, y2)
        tx3, ty3 = trans(ax2, x4, y2)
        tx4, ty4 = trans(ax1, x3, y1)
        plyg = mpl.patches.Polygon(np.array([[tx1, tx2, tx3, tx4],
                                             [ty1, ty2, ty3, ty4]]).T,
                                   facecolor=color, zorder=0, alpha=alpha)
        fig.add_artist(plyg)

    # for i, (name, (highpass, lowpass)) in enumerate(freq_bands.items()):
    for i, asig in enumerate(segment.analogsignals[1:]):
        highpass = asig.annotations['highpass_freq']
        lowpass = asig.annotations['lowpass_freq']
        x = np.linspace(highpass, lowpass, 50)
        ax[0].fill_between(x, .5, avg_psd_f(x), color=colors[i+1], alpha=alpha)
        add_polygon(i+2, highpass, lowpass, color=colors[i+1])

    # Brain inset
    inax = ax[0].inset_axes([.63,.6,.5,.5])
    imgplot = inax.imshow(cortex_img)
    inax.set_axis_off()

    # Signal
    for i, asig in enumerate(block.segments[0].analogsignals):
        asig = asig.time_slice(tw_start, tw_stop)
        ax[1].plot(asig.times.rescale('s'),
                   zscore(np.nanmean(np.real(asig.as_array()), axis=1)),
                   alpha=alpha, label=asig.name, color=colors[i])
    ax[1].axvline(t, color='k', ls=':')
    xticks = np.linspace((t-t_window/2).magnitude,
                         (t+t_window/2).magnitude, 11)[1::2]
    xtl = np.round(xticks - t.magnitude, decimals=1)
    xticklabels = [f'{xtl[0]}',
                   f'{xtl[1]}',
                   '',
                   f'+{xtl[3]}',
                   f'+{xtl[4]}']
    ax[1].set_xticks(xticks)
    ax[1].set_xticklabels(xticklabels, fontstyle='italic')
    ax[1].set_title('{:.2f}'.format(np.round(t, decimals=2)),
                    y=.98, fontsize=11)
    ax[1].set_xlim((t-t_window/2).magnitude, (t+t_window/2).magnitude)
    ax[1].set_ylim((-3,3))
    ax[1].set_ylabel('Relative Signal', weight='bold')
    ax[1].set_xlabel('Time [s]', weight='bold')
    ax[1].xaxis.labelpad = 0
    ax[1].xaxis.set_ticks_position('top')
    ax[1].xaxis.set_label_position('top')
    ax[1].xaxis.set_label_coords(.97, 1.04)
    ax[1].tick_params(axis="x", direction="in", pad= -18)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[1].yaxis.set_ticks_position('left')
    legend = ax[1].legend(loc=6, fontsize=11.5)
    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('white')
    return fig


if __name__ == '__main__':
    CLI = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    CLI.add_argument("--data", nargs='?', type=str, required=True)
    CLI.add_argument("--frame_folder",  nargs='?', type=str, required=True)
    CLI.add_argument("--frame_rate",  nargs='?', type=int, required=True)
    CLI.add_argument("--t_window",  nargs='?', type=float, required=True)
    CLI.add_argument("--cortex_img",  nargs='?', type=str, required=True)
    args = CLI.parse_args()

    with neo.NixIO(args.data) as nio:
        block = nio.read_block()

    block = AnalogSignal2ImageSequence(block)

    raw_asig = block.segments[0].analogsignals[0]

    frame_ids = stretch_to_framerate(t_start=raw_asig.t_start,
                                     t_stop=raw_asig.t_stop,
                                     num_frames=len(raw_asig),
                                     frame_rate=args.frame_rate*pq.Hz)

    cortex_img = mpl.image.imread(args.cortex_img)

    if not os.path.exists(args.frame_folder):
        os.makedirs(args.frame_folder)

    for frame_id in frame_ids:
        plot_figure(frame_id,
                    segment=block.segments[0],
                    t_window=args.t_window*pq.s,
                    cortex_img=cortex_img)
        plt.savefig(os.path.join(args.frame_folder,
                                 'frame_{}.png'.format(str(frame_id).zfill(5))
                                 )
                    )
        plt.close()
