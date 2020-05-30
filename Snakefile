import os
import json

configfile: 'config.yaml'

output_folder = os.path.join(os.getcwd(), 'output')
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

rgio_path = 'multielectrode_grasp/code/reachgraspio/'
odml_path = 'multielectrode_grasp/code/python-odml/'
session_id = 'i140703-001'

rule all:
    input:
        os.path.join(output_folder, 'frames.mp4')

rule get_data:
    output:
        lfp_data = '{path}.ns2',
        spike_data = '{path}-03.nev',
        metadata = '{path}.odml'
    params:
        dir = lambda wildcards: os.path.dirname(wildcards.path),
        lfp_data = lambda wildcards, output: os.path.basename(output.lfp_data),
        spike_data = lambda wildcards, output: os.path.basename(output.spike_data),
        metadata = lambda wildcards, output: os.path.basename(output.metadata)
    shell:
        """
        rm -rf multielectrode_grasp/
        git clone git@gin.g-node.org:/INT/multielectrode_grasp.git
        git-annex init
        cd {params.dir}
        git-annex get {params.lfp_data}
        git-annex get {params.spike_data}
        git-annex get {params.metadata}
        """

rule load_data:
    input:
        lfp_data = os.path.join('multielectrode_grasp',
                                'datasets', f'{session_id}.ns2'),
        spike_data = os.path.join('multielectrode_grasp',
                                  'datasets', f'{session_id}-03.nev'),
        odml = os.path.join('multielectrode_grasp',
                            'datasets', f'{session_id}.odml'),
        script = 'scripts/load_data.py'
    params:
        session = session_id,
        t_start = config['START_TIME'] - config['TIME_WINDOW']/2,
        t_stop  = config['STOP_TIME']  + config['TIME_WINDOW']/2,
        channels = [1, 96], # incl.
        grid_size = [10, 10],
        add_paths = f"export PYTHONPATH='{rgio_path}:{odml_path}:$PYTHONPATH'"
    output:
        data = os.path.join(output_folder, 'parsed_data.nix')
    shell:
        """
        {params.add_paths}
        python {input.script} --data {input.lfp_data} \
                              --odml {input.odml} \
                              --output {output.data} \
                              --channels {params.channels} \
                              --t_start {params.t_start} \
                              --t_stop {params.t_stop} \
                              --grid_size {params.grid_size}
        """

rule filter_and_transform:
    input:
        data = os.path.join('{path}', 'parsed_data.nix'),
        script = 'scripts/filter_and_transform.py'
    params:
        freq_bands = json.dumps({"Delta": config['DELTA_FREQ'],
                                 "Theta": config['THETA_FREQ'],
                                 "Alpha": config['ALPHA_FREQ'],
                                 "Beta" : config['BETA_FREQ']}),
        filter_order = 2
    output:
        data = os.path.join('{path}', 'processed_data.nix')
    shell:
        """
        python {input.script} --data {input.data} \
                              --output {output.data} \
                              --freq_bands '{params.freq_bands}' \
                              --filter_order {params.filter_order}
        """


rule plot_movie_frames:
    input:
        data = os.path.join('{path}', 'processed_data.nix'),
        cortex_img = 'images/cortex_location_croped.png',
        script = "scripts/plot_frames.py"
    params:
        frame_rate = config['FRAME_RATE'],
        t_window = config['TIME_WINDOW'],
        t_start = config['START_TIME'],
        t_stop = config['STOP_TIME']
    output:
        frame_folder = directory(os.path.join('{path}', 'frames'))
    shell:
        """
        python {input.script} --data "{input.data}" \
                              --frame_folder "{output.frame_folder}" \
                              --frame_rate {params.frame_rate} \
                              --t_window {params.t_window} \
                              --t_start {params.t_start} \
                              --t_stop {params.t_stop} \
                              --cortex_img {input.cortex_img}
        """

rule plot_movie:
    input:
        os.path.join('{path}', 'frames')
    output:
        os.path.join('{path}', 'frames.mp4')
    params:
        frame_path = lambda wildcards, input: os.path.join(input[0],
                                                           'frame_%05d.png'),
        quality = config['VIDEO_QUALITY'],
        scale_x = config['SCALE_X'],
        scale_y = config['SCALE_Y'],
        bitrate = config['BITRATE'],
        fps = config['FPS']
    shell:
        """
        ffmpeg -y -framerate {params.fps} \
               -i "{params.frame_path}" \
               -crf {params.quality} \
               -vb {params.bitrate} \
               -vcodec libx264 \
               -vf scale={params.scale_x}:{params.scale_y} \
               "{output}"
        """
