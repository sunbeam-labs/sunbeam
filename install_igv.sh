#!/bin/bash
set -e

# Quick install of IGV into the sunbeam/local directory.
# For use in Sunbeam the IGV path should be specified in the mapping section of
# the configuration file as 'igv_fp'.
install_igv() {
    DIR=$(readlink -f $(dirname $BASH_SOURCE))/local
    IGV_VER=2.3.68
    wget http://data.broadinstitute.org/igv/projects/downloads/IGV_${IGV_VER}.zip
    unzip IGV_${IGV_VER}.zip -d $DIR
    ln -s IGV_$IGV_VER $DIR/IGV
    # A symlink will confuse igv.sh so I'm using a wrapper script instead
    echo -e "#!/bin/bash\ncd $DIR/IGV && ./igv.sh \$@" > $DIR/IGV/igv
    chmod +x $DIR/IGV/igv
    export PATH=$DIR/IGV:$PATH
    command -v igv >/dev/null 2>&1 || { echo "IGV still isn't on the path, try installing manually"; exit 1; }
}

command -v igv >/dev/null 2>&1 || { echo "IGV not installed, installing now"; install_igv; }
