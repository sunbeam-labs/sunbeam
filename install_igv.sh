#!/bin/bash
set -e

# Quick install of IGV into the working directory.
# $DIR/IGV will still need to be separately added to PATH for this to work.
install_igv() {
    IGV_VER=2.3.92
    wget http://data.broadinstitute.org/igv/projects/downloads/IGV_${IGV_VER}.zip
    unzip IGV_${IGV_VER}.zip
    ln -s IGV_$IGV_VER IGV
    DIR=$(pwd)
    # A symlink will confuse igv.sh so I'm using a wrapper script instead
    echo -e "#!/bin/bash\ncd $DIR/$IGV_DIR && igv.sh \$@" > IGV/igv
    chmod +x IGV/igv
    export PATH=$DIR/IGV:$PATH
    command -v igv >/dev/null 2>&1 || { echo "IGV still isn't on the path, try installing manually"; exit 1; }
}

command -v igv >/dev/null 2>&1 || { echo "IGV not installed, installing now"; install_igv; }
