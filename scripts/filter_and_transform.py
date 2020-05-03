import sys
import os
import neo
import quantities as pq
import elephant
import numpy as np
import json
import argparse


if __name__ == '__main__':
    CLI = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    CLI.add_argument("--data", nargs='?', type=str, required=True)
    CLI.add_argument("--output",  nargs='?', type=str, required=True)
    CLI.add_argument("--freq_bands",  nargs='?', type=str, required=True)
    CLI.add_argument("--filter_order", nargs='?', type=int, default=2)
    args = CLI.parse_args()

    with neo.NixIO(args.data) as nio:
        block = nio.read_block()

    raw_asig = block.segments[0].analogsignals[0]
    del raw_asig.annotations['nix_name']

    freq_bands = json.loads(args.freq_bands)

    for i, (name, (highpass, lowpass)) in enumerate(freq_bands.items()):
        asig = elephant.signal_processing.butter(raw_asig,
                                                 highpass_freq=highpass*pq.Hz,
                                                 lowpass_freq=lowpass*pq.Hz,
                                                 order=2,
                                                 filter_function='filtfilt')
        hilbert_signal = elephant.signal_processing.hilbert(asig)
        hilbert_asig = neo.AnalogSignal(signal=hilbert_signal,
                                        units='dimensionless',
                                        dtype=np.complex,
                                        t_start=asig.t_start,
                                        sampling_rate=asig.sampling_rate,
                                        name=name,
                                        highpass_freq=highpass,
                                        lowpass_freq=lowpass,
                                        **asig.annotations)

        hilbert_asig.array_annotations = asig.array_annotations
        block.segments[0].analogsignals.append(hilbert_asig)

    with neo.NixIO(args.output) as nio:
        nio.write_block(block)
