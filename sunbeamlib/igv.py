"""
Helper functions for interfacing with the Integrative Genomics Viewer (IGV).

Dependencies:
 * IGV
 * xvfb
 * xdotool (for the socket-based method)
"""

from pathlib import Path
from snakemake import shell
import socket
import subprocess
import tempfile
import time
import os

def render(genome, bams, imagefile, seqID=None, igv_fp="igv", method="script"):
        """ Render an alignment to an image, given a genome and bam files.

        genome: path to a fasta file
        bams: list of path to a sorted, indexed bam file
        imagefile: path to the image to save
        seqID: (optional) sequence identifier to load from genome
        igv_fp: (optional) path to IGV executable
        method: (optional) method for controlling IGV process: "script" or "socket"

        The image file may be smaller than expected.  See
        igv_render_socket_nonblocking() for an attempt to enlarge the window
        before saving the image.
        """
        input_paths = [str(Path(bam).resolve()) for bam in bams]
        genome_path = str(Path(genome).resolve())
        output_path = str( Path('.').resolve() / Path(imagefile) )
        genome_cmd = 'genome ' + genome_path
        if seqID:
            genome_cmd += "\ngoto %s" % seqID
        igvcommands = ['new',
            genome_cmd,
            'load ' + ','.join(input_paths),
            'collapse',
            'snapshot ' + output_path,
            'exit']
        if method == "script":
            igvprefs = ["GENOME_LIST=;%s" % genome_path,
                "DEFAULT_GENOME_KEY=%s" % genome_path]
            _control_script(igvcommands, igv_fp, igvprefs)
        elif method == "socket":
            _control_socket(igvcommands, igv_fp)
        else:
            raise ValueError("method should be 'script' or 'socket', not '%s'" % method)

def _control_script(igvcommands, igv_fp, igvprefs):
        igvscript = tempfile.NamedTemporaryFile()
        igvscript.writelines(map(lambda x: bytes(x+'\n', 'ascii'), igvcommands))
        igvscript.flush()
        # If previous genome files listed in IGV's preferences are no longer
        # available, IGV will throw a null pointer exception at startup and
        # batch commands will fail.  So, we'll use a preferences override file
        # to list the genome file used here.  (This is probably an IGV bug.  We
        # should see if it happens in the latest release.)
        # I've also tried setting IGV.Bounds in an attempt to make the window
        # larger, but it doesn't seem to have any effect.
        igvprefsfile = tempfile.NamedTemporaryFile()
        igvprefstext = map(lambda x: bytes(x+'\n', 'ascii'), igvprefs)
        igvprefsfile.writelines(igvprefstext)
        igvprefsfile.flush()
        shell("xvfb-run -s '-screen 1 1920x1080x24' %s -o %s -b %s" % (igv_fp, igvprefsfile.name, igvscript.name))

def _control_socket(igvcommands, igv_fp):
        # Start up IGV.  Use a port between 10000 and the max available, based
        # on the PID of this process.  (TODO is using this pid safe?)
        port = 10000 + os.getpid()%(2**16-10000)
        xauth = "/tmp/xauth-%d" % os.getpid()
        igvproc = subprocess.Popen(["xvfb-run", "-l", "-f" , xauth, "-s", "-screen 1 1920x1080x24", str(igv_fp), "-p", str(port)])

        # Connect to running IGV
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        while True:
            try:
                s.connect(('localhost', port))
                break
            except ConnectionRefusedError:
                time.sleep(0.5)

        # Based on http://unix.stackexchange.com/questions/5999/ :
        # This should make the window as large as the virtual X display, but in
        # practice my screenshots aren't going over 1280 x 1296.
        display = ":99"
        shell("DISPLAY="+display+" XAUTHORITY="+xauth+" xdotool search --onlyvisible --name IGV windowsize --sync 100% 100%")

        # Generate screenshot
        s.sendall(bytes('\n'.join(igvcommands), 'ascii'))
        s.close()
        igvproc.wait()
