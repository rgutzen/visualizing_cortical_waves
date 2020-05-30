import sys
# sys.path.insert(0, '../multielectrode_grasp/code/python-odml/')
# sys.path.insert(0, '../multielectrode_grasp/code/reachgraspio/')
import os
import neo
import quantities as pq
import reachgraspio
import numpy as np
import argparse


def load_block(session_path, odml_dir, channels, t_start, t_stop):
    session = reachgraspio.ReachGraspIO(session_path, odml_directory=odml_dir)
    return session.read_block(nsx_to_load=2,
                              n_starts=t_start*pq.s, n_stops=t_stop*pq.s,
                              channels=range(channels[0], channels[1]+1),
                              units=None,
                              load_events=True, load_waveforms=False,
                              scaling='voltage',
                              correct_filter_shifts=True)


def merge_analogsingals(asigs):
    min_length = np.min([len(asig.times) for asig in asigs])
    max_length = np.max([len(asig.times) for asig in asigs])
    if min_length != max_length:
        print('Warning: the length of the analog signals differs '\
            + 'between {} and {} '.format(min_length, max_length)\
            + 'All signals will be cut to the same length and merged '\
            + 'into one AnalogSignal object.')

    if len(np.unique([asig.sampling_rate for asig in asigs])) > 1:
        print([asig.sampling_rate for asig in asigs])
        raise ValueError('The AnalogSignal objects have different '\
                       + 'sampling rates!')

    asig_array = np.zeros((min_length, len(asigs)))

    for channel_number, asig in enumerate(asigs):
        asig_array[:, channel_number] = np.squeeze(asig.as_array()[:min_length])

    merged_asig = neo.AnalogSignal(asig_array*asigs[0].units,
                                sampling_rate=asigs[0].sampling_rate,
                                t_start=asigs[0].t_start)
    for key in asigs[0].annotations.keys():
        try:
            merged_asig.array_annotations[key] = np.array([a.annotations[key] for a in asigs])
        except:
            e = sys.exc_info()[0]
            print('can not merge annotation ', key, e)
    return merged_asig


def channel2coords(channels, gridsize=(10,10)):
    channels -= 1  # starting from 0
    x_coords = np.floor_divide(channels, gridsize[0])
    y_coords = np.remainder(channels, gridsize[1])
    return x_coords, y_coords


if __name__ == '__main__':
    CLI = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    CLI.add_argument("--data", nargs='?', type=str, required=True)
    CLI.add_argument("--odml",  nargs='?', type=str, required=True)
    CLI.add_argument("--output",  nargs='?', type=str, required=True)
    CLI.add_argument("--channels",  nargs='+', type=int, default=[1,97])
    CLI.add_argument("--t_start", nargs='?', type=float, default=0)
    CLI.add_argument("--t_stop", nargs='?', type=float, default=10)
    CLI.add_argument("--grid_size", nargs='+', type=int, default=[10,10])
    args = CLI.parse_args()

    block = load_block(session_path=os.path.splitext(args.data)[0],
                       odml_dir=os.path.dirname(args.odml),
                       channels=args.channels,
                       t_start=args.t_start,
                       t_stop=args.t_stop)

    asig = merge_analogsingals(block.segments[0].analogsignals)
    x_coords, y_coords = channel2coords(
                                asig.array_annotations['connector_aligned_id'],
                                args.grid_size)
    asig.array_annotations['x_coords'] = x_coords
    asig.array_annotations['y_coords'] = y_coords
    asig.annotations['spatial_scale'] = block.annotations['electrodes_pitch']
    asig.name = 'Raw Signal'

    block.segments[0].analogsignals = [asig]

    # remove neo elements for compability
    del block.annotations['subject_birthday']
    block.channel_indexes = []

    # save data as nix
    with neo.NixIO(args.output) as nio:
        nio.write_block(block)
