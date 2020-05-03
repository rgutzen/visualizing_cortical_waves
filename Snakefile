import os
import json

output_folder = os.path.join(os.getcwd(), 'output')
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

rgio_path = 'multielectrode_grasp/code/reachgraspio/'
odml_path = 'multielectrode_grasp/code/python-odml/'
neo_path = 'multielectrode_grasp/code/python-neo/'

print(os.path.join(output_folder, 'parsed_data.nix'))

rule get_data:
    output:
        data = directory('multielectrode_grasp/')
    params:
        session = 'i140703-001'
    shell:
        """
        git clone git@gin.g-node.org:/INT/multielectrode_grasp.git
        git-annex get multielectrode_grasp/datasets/{params.seesion}.nev
        git-annex get multielectrode_grasp/datasets/{params.session}.ns2
        git-annex get multielectrode_grasp/datasets/{params.session}.odml
        """

rule load_data:
    input:
        data = os.path.join('multielectrode_grasp', 'datasets', 'i140703-001.ns2'),
        odml = os.path.join('multielectrode_grasp/', 'datasets', 'i140703-001.odml'),
        script = 'scripts/load_data.py'
    params:
        t_start = 0,  # s
        t_stop = 2,  # s
        channels = [1, 96], # incl.
        grid_size = [10, 10],
        add_paths = f"export PYTHONPATH='{rgio_path}:{odml_path}:$PYTHONPATH'"
    output:
        data = os.path.join(output_folder, 'parsed_data.nix')
    shell:
        """
        export PYTHONPATH='{rgio_path}:{odml_path}:$PYTHONPATH'
        python {input.script} --data {input.data} \
                              --odml {input.odml} \
                              --output {output.data} \
                              --channels {params.channels} \
                              --t_start {params.t_start} \
                              --t_stop {params.t_stop} \
                              --grid_size {params.grid_size}
        """

rule filter_and_transform:
    input:
        data = rules.load_data.output.data,
        script = 'scripts/filter_and_transform.py'
    params:
        freq_bands = json.dumps({"Delta": ( 1,  4),
                                 "Theta": ( 4,  8),
                                 "Alpha": ( 8, 12),
                                 "Beta" : (12, 25)}),
        filter_order = 2
    output:
        data = os.path.join(output_folder, 'processed_data.nix')
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
        cortex_img = 'cortex_location.png',
        script = "scripts/plot_frames.py"
    params:
        frame_rate = 400,  # Hz
        t_window = 1  # s
    output:
        frame_folder = directory(os.path.join('{path}', 'frames'))
    shell:
        """
        python {input.script} --data "{input.data}" \
                              --frame_folder "{output.frame_folder}" \
                              --frame_rate {params.frame_rate} \
                              --t_window {params.t_window} \
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
        quality = 5,
        scale_x = 720,
        scale_y = 720,
        bitrate = '20M',
        fps = 20
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
