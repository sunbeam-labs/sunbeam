import csv
import re
import string
from pathlib import Path


class SampleList:
    def __init__(
        self, fp: Path = None, paired_end: bool = True, format_str: str = None
    ):
        self.paired_end = paired_end
        self.samples = {}

        if fp:
            if fp.is_dir():
                self.samples = self.load_from_dir(fp, format_str)
            elif fp.is_file():
                self.samples = self.load_from_file(fp)

        if self.samples:
            self.check()

    def check(self):
        """
        Check the sample list for duplicates, missing values, and invalid values.
        """
        seen = set()
        for sample, data in self.samples.items():
            if sample in seen:
                raise ValueError(f"Duplicate sample name found: {sample}")
            seen.add(sample)
            # Check for disallowed characters in sample names ("/", " ", ",")
            if re.search(r"[ /,]", sample):
                raise ValueError(f"Invalid characters in sample name: {sample}")

            if "r1" not in data or (self.paired_end and "r2" not in data):
                raise ValueError(f"Missing read files for sample: {sample}")

            # Custom function for checking file existence (handle local and remote data sources)

    def load_from_file(self, fp: Path) -> dict[str, dict[str, str]]:
        """
        Load a sample list from a file.
        """
        samples = {}
        with open(fp, "r") as f:
            reader = csv.DictReader(f, delimiter=",", fieldnames=["sample", "r1", "r2"])
            for row in reader:
                samples[row["sample"]] = {"r1": row["r1"]}
                if self.paired_end:
                    samples[row["sample"]]["r2"] = row["r2"]

        return samples

    def load_from_dir(
        self, fp: Path, format_str: str = None
    ) -> dict[str, dict[str, str]]:
        """
        Load a sample list from a directory.
        """
        samples = {}
        fnames = [
            f for f in fp.iterdir() if f.is_file() and f.name.endswith(".fastq.gz")
        ]
        if len(fnames) == 0:
            raise ValueError("No gzipped FASTQ files found in the directory.")

        pattern = (
            self.format_string_to_regex(format_str)
            if format_str
            else self.guess_format_string(fnames)
        )
        print(f"Guessed format string: {pattern}")
        for f in fnames:
            match = re.match(pattern, f.name)
            if not match:
                raise ValueError(
                    f"Filename {f.name} does not match format string {pattern}"
                )

            sample = match.group("sample")
            if sample not in samples:
                samples[sample] = {"r1": None, "r2": None}

            if "rp" in match.groupdict():
                if match.group("rp") == "1":
                    samples[sample]["r1"] = f
                else:
                    samples[sample]["r2"] = f
            else:
                samples[sample]["r1"] = f

        return samples

    def guess_format_string(self, fnames: list[Path]) -> re.Pattern:
        patterns = [
            re.compile(r"(?P<sample>.+)_R(?P<rp>[12])(_[A-Za-z0-9]+)*\.fastq\.gz"),
            re.compile(r"(?P<sample>.+)_(?P<rp>[12])(_[A-Za-z0-9]+)*\.fastq\.gz"),
            re.compile(
                r"(?P<sample>.+)_S\d+_L\d+_R(?P<rp>[12])_001\.fastq\.gz"
            ),  # Illumina-style
            re.compile(r"(?P<sample>.+)_R(?P<rp>[12])_.*\.fastq\.gz"),
            re.compile(r"(?P<sample>.+)_(?P<rp>[12])_.*\.fastq\.gz"),
            re.compile(r"(?P<rp>[12])_(?P<sample>.+)\.fastq\.gz"),
            re.compile(r"R(?P<rp>[12])_(?P<sample>.+)\.fastq\.gz"),
        ]

        if not self.paired_end:
            pattern = re.compile(r"(?P<sample>.+)\.fastq\.gz")
            if all(pattern.match(f.name) for f in fnames):
                return pattern
        else:
            for pattern in patterns:
                if all(pattern.match(f.name) for f in fnames):
                    return pattern

        raise ValueError(f"Could not guess format string from filenames: {fnames}")

    @staticmethod
    def format_string_to_regex(format_str: str) -> re.Pattern:
        """
        Convert a format string like '{sample}_R{rp}.fastq.gz' to a regex pattern.
        """
        regex = ""
        formatter = string.Formatter()

        for literal_text, field_name, format_spec, conversion in formatter.parse(
            format_str
        ):
            regex += re.escape(literal_text)
            if field_name == "sample":
                regex += r"(?P<sample>.+?)"
            elif field_name == "rp":
                regex += r"(?P<rp>[12])"
            elif field_name:  # fallback
                regex += rf"(?P<{field_name}>.+?)"

        return re.compile(regex)

    def to_file(self, fp: Path):
        """
        Write the sample list to a file.
        """
        with open(fp, "w") as f:
            writer = csv.DictWriter(f, delimiter=",", fieldnames=["sample", "r1", "r2"])
            for sample, data in self.samples.items():
                row = {"sample": sample, "r1": data["r1"]}
                if self.paired_end:
                    row["r2"] = data["r2"]
                writer.writerow(row)

    def generate_subset(self, func: callable) -> "SampleList":
        """
        Generate a subset of the sample list based on a function. The function takes three args:
        sample_name, r1, and r2, and should return True if the sample should be included in the subset.
        """
        subset = SampleList()
        for sample, data in self.samples.items():
            if func(sample, data["r1"], data.get("r2")):
                subset.samples[sample] = data
        return subset

    def get_samples(self):
        """
        Getting samples in a format that is backwards compatible
        Convert "r1" and "r2" to "1" and "2", make sure all fields are strings
        """
        samples = {}
        for sample, data in self.samples.items():
            samples[sample] = {k: str(v) for k, v in data.items()}
            if "r1" in samples[sample]:
                samples[sample]["1"] = samples[sample].pop("r1")
            if "r2" in samples[sample]:
                samples[sample]["2"] = samples[sample].pop("r2")

        return samples
