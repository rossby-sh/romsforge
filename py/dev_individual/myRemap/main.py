import sys
sys.path.append('C:/Users/ust21/shjo/projects/myROMS/remap/')
from remap_module import process_files
import argparse
import os
import yaml


def main():
    parser = argparse.ArgumentParser(description="Remap 15km ROMS output to 5km grid using weights defined in a config file.")
    parser.add_argument('--config', required=True, help='Path to config YAML file')
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    os.makedirs(config['output_dir'], exist_ok=True)

    weight_paths = {
        'rho': config['weight_rho'],
        'u': config['weight_u'],
        'v': config['weight_v']
    }

    process_files(
        avg_dir=config['avg_dir'],
        grid01=config['grid01'],
        grid02=config['grid02'],
        weight_paths=weight_paths,
        output_dir=config['output_dir']
    )


if __name__ == '__main__':
    main()