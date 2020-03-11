import os

_FLAGSTAT_CATEGORIES = [
    'Total',
    'Secondary',
    'Supplementary',
    'Duplicates',
    'Mapped',
    'Paired in sequencing',
    'Read1',
    'Read2',
    'Properly paired',
    'With itself and mate mapped',
    'Singletons',
    'With mate mapped to a different chr',
    'With mate mapped to a different chr (mapQ>=5)'
]

FLAGSTAT_TABLE_HEADER = ['Sample']
for category in _FLAGSTAT_CATEGORIES:
    FLAGSTAT_TABLE_HEADER.append(f"{category} (QC passed)")
    FLAGSTAT_TABLE_HEADER.append(f"{category} (QC failed)")


class FlagstatParser:

    def __init__(self, file_path):
        self.file_path = file_path

    def _parse(self):
        parsed_list = []
        with open(self.file_path, 'r') as file:
            for line in file:
                elements = line.rstrip().split('+')
                qc_passed = elements[0].strip()
                qc_failed = elements[1].split(' ')[1].strip()
                parsed_list.append(qc_passed)
                parsed_list.append(qc_failed)
        self._parsed_list = parsed_list

    def get_list(self):
        if getattr(self, '_parsed_list', None) is None:
            self._parse()
        return self._parsed_list


def generate_mapping_flagstat(flagstat_file_paths, output_file_path, sep=','):
    with open(output_file_path, "w") as output_file:
        print(sep.join(FLAGSTAT_TABLE_HEADER), file=output_file)
        for file_path in flagstat_file_paths:
            sample_name = os.path.basename(file_path).split('.flagstat')[0]
            print(sep.join([sample_name] + FlagstatParser(file_path).get_list()), file=output_file)
