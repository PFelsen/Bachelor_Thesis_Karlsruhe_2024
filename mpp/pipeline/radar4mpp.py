import argparse
import json
import os


def update_metadata_json(tag):
    with open('metadata.json', 'r') as file:
        data = json.load(file)
        data['descriptiveMetadata']['title'] = 'Mpp {}'.format(tag)
        data['descriptiveMetadata']['descriptions']['description'][0]['value'] = \
            "Release notes: https://git.scc.kit.edu/mpp/mpp/-/tags/{}".format(tag)
        data['descriptiveMetadata']['software']['softwareType'][0]\
            ['softwareName'][0]['softwareVersion'] = tag

    os.remove('metadata.json')
    with open('metadata.json', 'w') as file:
        json.dump(data, file, indent=2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser('radar4mpp')

    parser.add_argument('-t', '--tag', type=str, default='X.X.X', help='tag name')

    args = parser.parse_args()

    update_metadata_json(args.tag)
