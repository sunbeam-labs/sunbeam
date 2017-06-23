"""
Helper functions for interfacing with the Integrative Genomics Viewer (IGV).

Dependencies:
 * IGV
 * xvfb
 * xdotool (for the socket-based method)
"""

from pathlib import Path
from snakemake import shell
from Bio import SeqIO
import socket
import subprocess
import tempfile
import time
import os

def render(genome, bams, imagefile, seqID=None, igv_fp="igv", method="script", igv_prefs=None):
        """ Render an alignment to an image, given a genome and bam files.

        genome: path to a fasta file
        bams: list of path to a sorted, indexed bam file
        imagefile: path to the image to save
        seqID: (optional) sequence identifier to load from genome
        igv_fp: (optional) path to IGV executable
        method: (optional) method for controlling IGV process: "script" or "socket"
        igv_prefs: (optional) dictionary of IGV preferences to override

        The image file may be smaller than expected.  See
        igv_render_socket_nonblocking() for an attempt to enlarge the window
        before saving the image.
        """
        igv_prefs = igv_prefs or {}
        genome = str(genome)
        imagefile = str(imagefile)
        input_paths = [str(Path(str(bam)).resolve()) for bam in bams]
        genome_path = str(Path(genome).resolve())
        output_path = str( Path('.').resolve() / Path(imagefile) )
        # build a "seqID:1:length" string to force IGV to display the full
        # segment.  If no segment was given it will default to the first.
        goto_locus = _seq_length(genome, seqID)
        goto_locus = "%s:%s-%s" % (goto_locus[0], '1', goto_locus[1])
        igvcommands = ['new',
            'genome ' + genome_path,
            'goto ' + goto_locus,
            'load ' + ','.join(input_paths),
            'collapse',
            'snapshot ' + output_path,
            'exit']
        if method == "script":
            # If previous genome files listed in IGV's preferences are no
            # longer available, IGV will throw a null pointer exception at
            # startup and batch commands will fail.  So, we'll use a
            # preferences override file to list the genome file used here.
            # (This is probably an IGV bug.  We should see if it happens in the
            # latest release.) I've also tried setting IGV.Bounds in an attempt
            # to make the window larger, but it doesn't seem to have any
            # effect.
            igv_prefs["GENOME_LIST"] = ";" + genome_path
            igv_prefs["DEFAULT_GENOME_KEY"] = genome_path
            _control_script(igvcommands, igv_fp, igv_prefs)
        elif method == "socket":
            _control_socket(igvcommands, igv_fp, igv_prefs)
        else:
            raise ValueError("method should be 'script' or 'socket', not '%s'" % method)

def _control_script(igvcommands, igv_fp, igv_prefs):
        igvscript = tempfile.NamedTemporaryFile()
        igvscript.writelines(map(lambda x: bytes(x+'\n', 'ascii'), igvcommands))
        igvscript.flush()
        igvprefsfile = _write_prefs(igv_prefs)
        shell("xvfb-run -a -s '-screen 1 1920x1080x24' %s -o %s -b %s" % (igv_fp, igvprefsfile.name, igvscript.name))

def _control_socket(igvcommands, igv_fp, igv_prefs):
        igvprefsfile = _write_prefs(igv_prefs)
        # Start up IGV.  Use a port between 10000 and the max available, based
        # on the PID of this process.  (TODO is using this pid safe?)
        port = 10000 + os.getpid()%(2**16-10000)
        xauth = "/tmp/xauth-%d" % os.getpid()
        xvfb_cmdline = ["xvfb-run", "-a", "-l", "-f" , xauth, "-s", "-screen 1 1920x1080x24"]
        igv_cmdline = [ str(igv_fp), "-p", str(port), "-o", igvprefsfile.name ]
        igvproc = subprocess.Popen(xvfb_cmdline + igv_cmdline)

        # Connect to running IGV
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        while True:
            try:
                s.connect(('localhost', port))
                break
            except ConnectionRefusedError:
                time.sleep(0.5)

        # Figure out what X11 display the IGV process is using.  It should be
        # the second child of the initial xfvb-run process (the first child
        # being Xvfb, I think.)
        with open("/proc/%s/task/%s/children" % (igvproc.pid, igvproc.pid)) as f:
            child_pid = f.read().split()[1]
        with open("/proc/%s/environ" % child_pid) as f:
            env_vars = [env_var.split('=', 1) for env_var in f.read().split('\x00') ]
            env_vars = {env_var[0]: env_var[1] for env_var in env_vars if len(env_var)==2}
        display = env_vars['DISPLAY']
        # Based on http://unix.stackexchange.com/questions/5999/ :
        # This should make the window as large as the virtual X display, but in
        # practice my screenshots aren't going over 1280 x 1296.
        shell("DISPLAY="+display+" XAUTHORITY="+xauth+" xdotool search --onlyvisible --name IGV windowsize --sync 100% 100%")

        # Generate screenshot
        s.sendall(bytes('\n'.join(igvcommands), 'ascii'))
        s.close()
        igvproc.wait()

def _write_prefs(igv_prefs):
        igvprefsfile = tempfile.NamedTemporaryFile()
        for k,v in igv_prefs.items():
            igvprefsfile.write(bytes("%s=%s\n" % (k, v), 'ascii'))
        igvprefsfile.flush()
        return igvprefsfile

def _seq_length(fasta_fp, seqID=None):
    """Give the ID and sequence length of the first (or specified) sequence. """
    record_lengths = [[r.id, len(r.seq)] for r in SeqIO.parse(fasta_fp, "fasta")]
    if not seqID:
        return(record_lengths[0])
    record_length_match = [r for r in record_lengths if r[0] == seqID]
    return(record_length_match[0])
