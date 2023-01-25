* [ ] I have run `bash tests/test.sh` on a local deployment and the tests passed successfully
* [ ] If this adds a new output file, I have added a check to tests/targets.txt
* [ ] If this fixes a bug, I have added or modified an appropriate test in tests/test.sh
* [ ] If this adds or modifies a rule that uses FASTQ files, the input accepts gzipped FASTQ and outputs gzipped FASTQ, or marks uncompressed FASTQ output as `temp`

If this is for a release:
* [ ] I have updated documentation
* [ ] I have updated the hardcoded version at the top of `install.sh` to match what this release's version will be
* [ ] I have created a release archive that will be attached to this release