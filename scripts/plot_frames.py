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


def stretch_to_framerate(t_start, t_stop, times, num_frames, frame_rate=None):
    if frame_rate is None:
        frame_ids = np.arange(num_frames, dtype=int)
    else:
        new_num_frames = (t_stop.rescale('s').magnitude
                        - t_start.rescale('s').magnitude) \
                        * args.frame_rate
        frame_ids = np.linspace(0, num_frames-1, int(new_num_frames), dtype=int)

    id_offset = np.argmax(t_start <= times.rescale('s'))
    frame_ids += id_offset
    return frame_ids


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


def smooth_frames(frames, sigma):
    ndim = frames.ndim
    if ndim == 2:
        frames = frames[np.newaxis,:,:]
    # replace nan sites by median
    if np.isfinite(frames).all():
        nan_sites = None
    else:
        # assume constant nan sites over time
        nan_sites = np.where(np.bitwise_not(np.isfinite(frames[0])))
        if np.iscomplexobj(frames):
            frames[:,nan_sites[0],nan_sites[1]] = np.nanmedian(np.real(frames)) \
                                                + np.nanmedian(np.imag(frames))*1j
        else:
            frames[:,nan_sites[0],nan_sites[1]] = np.nanmedian(frames)

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

    if ndim == 2:
        return frames[0,:,:]
    else:
        return frames

## Plotting functions

def plot_array(frame, cmap=plt.cm.gray, vmin=None, vmax=None, ax=None, zorder=3):
    if ax is None:
        fig, ax = plt.subplots()
    if cmap is None:
        cmap = plot.cm.gray
    else:
        cmap = plt.get_cmap(cmap)

    img = ax.imshow(frame, interpolation='nearest',
                    cmap=cmap, vmin=vmin, vmax=vmax, zorder=zorder)
    ax.axis('image')
    ax.set_xticks([])
    ax.set_yticks([])
    return img


def interpolate_frame(frame, resolution=(100,100), kind='cubic'):
    dim_x, dim_y = frame.shape
    if np.iscomplexobj(frame):
        re_frame = interpolate_frame(np.real(frame),
                                     resolution=resolution,
                                     kind=kind)
        im_frame = interpolate_frame(np.imag(frame),
                                     resolution=resolution,
                                     kind=kind)
        return re_frame + 1j*im_frame
    else:
        f = interp2d(np.arange(dim_x), np.arange(dim_y), np.nan_to_num(frame), kind=kind)
        frame = f(np.linspace(0, dim_x-1, resolution[0]),
                  np.linspace(0, dim_y-1, resolution[1]))
    return frame


def plot_psd(asig, freq_res, xlim, ylim, ax, color='k', interpolate=None):
    freqs, psd = elephant.spectral.welch_psd(asig, freq_res=freq_res, overlap=.7)
    avg_psd = np.mean(psd, axis=0)

    if interpolate is None:
        handle = ax.loglog(psd, freqs, linewidth=3, color=color)
    else:
        fvalues = np.linspace(ylim[0], ylim[1], 200)
        avg_psd_f = interp1d(freqs, avg_psd, kind=interpolate)
        handle = ax.loglog(avg_psd_f(fvalues), fvalues, linewidth=3,
                              color=color, label='Average')

    ax.set_xlim(xlim)
    ax.set_ylim((ylim))
    ax.set_ylabel('Frequency [Hz]', weight='bold')
    ax.set_axis_off()
    ax.annotate('', xy=(0, -0.01), xycoords='axes fraction', xytext=(.4, -0.01),
                arrowprops=dict(arrowstyle="->", color='k'))
    ax.text(.04, -.04, 'Log Spectral Power', weight='bold',
            transform=ax.transAxes, fontsize=12)
    ax.annotate('', xy=(0, -0.01), xycoords='axes fraction', xytext=(0, .2),
                arrowprops=dict(arrowstyle="<-", color='k'))
    ax.text(-.06, .01, 'Log Frequency', weight='bold', rotation=90,
            transform=ax.transAxes, fontsize=12)

    return avg_psd_f


def plot_raw_signal(asig, t_start, t_stop, color, ax):
    start_i = np.argmax(asig.times >= t_start)
    stop_i  = np.argmax(asig.times >= t_stop)+1
    plot_signal = zscore(np.nanmean(np.real(asig.as_array()), axis=1))
    handle = ax.plot(asig.times.rescale('s')[start_i:stop_i],
                     plot_signal[start_i:stop_i],
                     linewidth=3, alpha=.2, label=asig.name, color=color)

    ax.xaxis.set_ticks_position('bottom')
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_yticks([])
    ax.annotate('', xy=(.8, -0.01), xycoords='axes fraction', xytext=(1, -0.01),
                arrowprops=dict(arrowstyle="<-", color='k'))
    ax.text(.96, -.04, 'Time', weight='bold', ha='right',
               transform=ax.transAxes, fontsize=12)
    ax.annotate('', xy=(1, -0.01), xycoords='axes fraction', xytext=(1, .2),
                arrowprops=dict(arrowstyle="<-", color='k'))
    ax.text(1.01, 0.04, 'Signal', weight='bold', rotation=90,
               transform=ax.transAxes, fontsize=12)
    return handle


def add_polygon(ax1, ax2, y1, y4, bend, res, color=None, alpha=1):
    trans = lambda a, x, y: plt.gcf().transFigure.inverted().transform(a.transData.transform([x, y]))

    x1 = bend
    x2 = 0.345*res - 0.5 # (50, 0.335/ 100, 34)
    y2 = res
    x3 = x2
    y3 = 0
    x4 = x1
    tx1, ty1 = trans(ax1, x1, y1)
    tx2, ty2 = trans(ax2, x2, y2)
    tx3, ty3 = trans(ax2, x3, y3)
    tx4, ty4 = trans(ax1, x4, y4)
    plyg = mpl.patches.Polygon(np.array([[tx1, tx2, tx3, tx4],
                                         [ty1, ty2, ty3, ty4]]).T,
                               facecolor=color, zorder=1, alpha=alpha)
    plt.gcf().add_artist(plyg)
    return None


def plot_events(event, ax, xlim):
    event_keys = ['WS-ON', 'CUE-ON', 'GO-ON', 'SR', 'HEplat-ON','RW-ON']
    event_labels = ['Start Cue', 'Grip-Type Cue', 'Grip-Force Cue', 'Grasping', 'Holding', 'Reward']
    for key, label in zip(event_keys, event_labels):
        index = np.argmax(event.labels == key)
        event_t = event.times[index].rescale('s')
        if xlim[0] < event_t.magnitude < xlim[1]:
            ax.plot(event_t, 0, marker=5, markersize=15, color='0.2', linestyle=None)
            ax.text(event_t-0.02*pq.s, 0, label, va='center', ha='right', fontsize=12, backgroundcolor='w')
    return None


def plot_colorbar(ax, img, bbox):
    ax.set_axis_off()
    cax = ax.inset_axes(bbox)
    cbar = plt.colorbar(img, ax=ax, cax=cax, orientation='horizontal',
                        ticks=[-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
    cbar.set_label(label='Phase', weight='bold', fontsize=12)
    cax.xaxis.labelpad = -40
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.set_xticklabels([r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$'], fontsize=11)
    return cbar


def plot_figure(frame_id, raw_asig, freq_asigs, freq_imgseqs, events,
                t_window, cortex_img):
    fig = plt.figure(figsize=(16,9))
    gs = mpl.gridspec.GridSpec(nrows=5, ncols=3, height_ratios=[.07,1,1,1,1], width_ratios=[.8,.9,1.5])
    gs.update(hspace=0.18)
    gs.update(wspace=-0.4)
    ax = [fig.add_subplot(gs[0:, 0]), # freq
          fig.add_subplot(gs[0:, 2]), # signal box
          fig.add_subplot(gs[1, 2]), # signals
          fig.add_subplot(gs[2, 2]),
          fig.add_subplot(gs[3, 2]),
          fig.add_subplot(gs[4, 2]),
          fig.add_subplot(gs[0, 1]), # colorbar
          fig.add_subplot(gs[1, 1]), # frames
          fig.add_subplot(gs[2, 1]),
          fig.add_subplot(gs[3, 1]),
          fig.add_subplot(gs[4, 1]),
         ]

    colors = ['k',
              '#DED93E', # Delta
              '#8BCD50', # Theta
              '#1D741B', # Alpha
              '#8FA01F'  # Beta
             ]
    alpha = .7

    t = raw_asig.times[frame_id].rescale('s')
    tw_start = max([raw_asig.t_start.rescale('s'), t-t_window/2])
    tw_stop = tw_start + t_window
    if tw_stop > raw_asig.t_stop:
        tw_stop = raw_asig.t_stop
        tw_start = tw_stop - t_window

    # Raw Signal
    plot_raw_signal(raw_asig,
                    t_start=tw_start,
                    t_stop=tw_stop,
                    color=colors[0],
                    ax=ax[1])

    ax[1].axvline(t, color='k', ls=':', zorder=0)
    ax[1].text(t, 3.01, '{:.2f}'.format(np.round(t, decimals=2)), ha='center', fontsize=12)
    xticks = np.linspace((t-t_window/2).magnitude,
                         (t+t_window/2).magnitude, 11)[3::2]
    xtl = np.round(xticks - t.magnitude, decimals=1)
    xticklabels = [f'- {abs(xtl[0])} s', '0.0 s', f'+ {xtl[2]} s', f'+ {xtl[3]} s']
    ax[1].set_xticks(xticks)
    ax[1].set_xticklabels(xticklabels, fontstyle='italic', fontsize=11)
    ax[1].tick_params(axis="x",direction="in", pad=-15)
    plt.setp(ax[1].get_xticklabels(), backgroundcolor="white")

    ax[1].set_xlim((t-t_window/3).magnitude, (t+t_window/2).magnitude)
    ax[1].set_ylim((-3,3))

    # Frequency Component Signals
    start_i = np.argmax(freq_asigs[0].times >= tw_start)
    stop_i  = np.argmax(freq_asigs[0].times >= tw_stop)+1
    times = freq_asigs[0].times.rescale('s')[start_i:stop_i]
    for i, asig in enumerate(freq_asigs):
        ax[5-i].patch.set_visible(False)
        ax[5-i].axison = False
        plot_signal = zscore(np.nanmean(np.real(asig.as_array()), axis=1))
        ax[5-i].plot(times,
                     plot_signal[start_i:stop_i],
                     linewidth=1.5, label=asig.name, color=colors[i+1])
        ax[5-i].set_xlim((t-t_window/3).magnitude, (t+t_window/2).magnitude)
        ax[5-i].set_ylim((-4, 4))

    color_patches = [mpl.patches.Patch(color=c) for c in colors[1:]]
    grey_patch = mpl.patches.Patch(color='0.8')
    legend = ax[2].legend([grey_patch, tuple(color_patches)],
                          ['Full Signal', 'Frequency Components'],
                          markerscale=42, numpoints=1,
                          fontsize=12, loc='upper right',
                          framealpha=1,
                          handler_map={tuple: mpl.legend_handler.HandlerTuple(ndivide=None)},
                          bbox_to_anchor=(1,1.3))
    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('white')

    # Electrode Array Frames
    resolution = (30,30)
    for i, imgseq in enumerate(freq_imgseqs):
        frame = imgseq.as_array()[frame_id]
        frame = smooth_frames(frame, .6)
        frame = interpolate_frame(frame, resolution=resolution, kind='linear')
        frame = np.angle(frame)
        img = plot_array(frame,
                         cmap='twilight',
                         ax=ax[10-i],
                         vmin=-np.pi,
                         vmax=np.pi,
                        )
        ax[10-i].patch.set_visible(False)
    ax[10].set_xlabel(r'$\longleftarrow\; 4mm\; \longrightarrow$', fontsize=11)

    # Colorbar
    plot_colorbar(ax[6], img=img, bbox=[.27,0,.5,1])

    # PSD
    avg_psd_f = plot_psd(asig=raw_asig.time_slice(tw_start, tw_stop),
                         freq_res=1/t_window,
                         ax=ax[0],
                         xlim=(4*10**3, 0.1),
                         ylim=(0.95, 35),
                         interpolate='cubic',
                         color='0.2')

    # Freq areas
    bend = 1 # bend from fill_between to polygon
    plygns = []

    for i, asig in enumerate(freq_asigs):
        highpass = asig.annotations['highpass_freq']
        lowpass = asig.annotations['lowpass_freq']
        x = np.linspace(highpass, lowpass, 50)
        ax[0].fill_betweenx(x, bend, avg_psd_f(x), color=colors[i+1], alpha=0.9)
        add_polygon(ax[0], ax[10-i], highpass, lowpass, bend=bend, res=resolution[1], color=colors[i+1], alpha=0.9)

        logdiff = (np.log10(lowpass)-np.log10(highpass))
        pos = 10**(np.log10(highpass) + logdiff/2)
        ax[0].text(bend*1, pos, asig.name, rotation=90, va='center', ha='right',
                   weight='bold', fontsize=12, zorder=1)

        ax[0].text(bend*1.7, highpass*1.01, f'{highpass} Hz', ha='right', fontsize=11, alpha=.7)
    ax[0].text(bend*1.7, lowpass*1.01, f'{lowpass} Hz', ha='right', fontsize=11, alpha=.7)

    # Brain inset
    inax = ax[0].inset_axes([-.05,.65,.5,.5])
    imgplot = inax.imshow(cortex_img, zorder=3)
    inax.set_axis_off()

    plot_events(event=events, ax=ax[1],
                xlim=((t-t_window/4).magnitude, (t+t_window/2).magnitude))
    return fig


if __name__ == '__main__':
    CLI = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    CLI.add_argument("--data", nargs='?', type=str, required=True)
    CLI.add_argument("--frame_folder",  nargs='?', type=str, required=True)
    CLI.add_argument("--frame_rate",  nargs='?', type=int, required=True)
    CLI.add_argument("--t_window",  nargs='?', type=float, required=True)
    CLI.add_argument("--t_start",  nargs='?', type=float, required=True)
    CLI.add_argument("--t_stop",  nargs='?', type=float, required=True)
    CLI.add_argument("--cortex_img",  nargs='?', type=str, required=True)
    args = CLI.parse_args()

    with neo.NixIO(args.data) as nio:
        block = nio.read_block()

    block = AnalogSignal2ImageSequence(block)

    raw_asig = block.segments[0].analogsignals[0]

    num_avail_frames = len(raw_asig.time_slice(args.t_start*pq.s,
                                               args.t_stop*pq.s))
    frame_ids = stretch_to_framerate(t_start=args.t_start*pq.s,
                                     t_stop=args.t_stop*pq.s,
                                     times=raw_asig.times.rescale('s'),
                                     num_frames=num_avail_frames,
                                     frame_rate=args.frame_rate*pq.Hz)

    cortex_img = mpl.image.imread(args.cortex_img)
    trial_events = block.segments[0].events[2].time_slice(args.t_start*pq.s,
                                                          args.t_stop*pq.s)

    if not os.path.exists(args.frame_folder):
        os.makedirs(args.frame_folder)

    for i, frame_id in enumerate(frame_ids):
        plot_figure(frame_id,
                    raw_asig=raw_asig,
                    freq_asigs=block.segments[0].analogsignals[1:],
                    freq_imgseqs=block.segments[0].imagesequences[1:],
                    events=trial_events,
                    t_window=args.t_window*pq.s,
                    cortex_img=cortex_img)

        plt.savefig(os.path.join(args.frame_folder,
                                 'frame_{}.png'.format(str(i).zfill(5))
                                 ),
                    bbox_inches='tight')
        plt.close()
