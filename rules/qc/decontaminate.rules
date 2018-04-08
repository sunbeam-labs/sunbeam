import subprocess
from sunbeamlib.decontam import get_mapped_reads
from collections import OrderedDict

rule all_decontam:
    input:
        expand(
            str(QC_FP/'decontam'/'{sample}_{rp}.fastq.gz'),
            sample=Samples.keys(), rp=Pairs)

ruleorder: build_host_index > build_genome_index
        
rule build_host_index:
    input:
        str(Cfg['qc']['host_fp']/'{host}.fasta')
    output:
        str(Cfg['qc']['host_fp']/'{host}.fasta.amb')
    params:
        host='{host}',
        index_fp=str(Cfg['qc']['host_fp'])
    shell:
        "cd {Cfg[qc][host_fp]} && bwa index {input}"

rule align_to_host:
    input:
        reads = expand(
            str(QC_FP/'cleaned'/'{{sample}}_{rp}.fastq.gz'),
            rp = Pairs),
        index = str(Cfg['qc']['host_fp']/'{host}.fasta.amb')
    output:
        temp(str(QC_FP/'decontam'/'intermediates'/'{host}'/'{sample}.bam'))
    threads:
        Cfg['qc']['threads']
    params:
        sam = temp(str(QC_FP/'decontam'/'intermediates'/'{host}'/'{sample}.sam')),
        index_fp = str(Cfg['qc']['host_fp'])
    shell:
        """
        bwa mem -M -t {threads} \
        {params.index_fp}/{wildcards.host}.fasta \
        {input.reads} -o {params.sam} && \
        samtools view -bSF4 {params.sam} > {output} && \
        rm {params.sam}
        """

rule get_mapped_reads:
    input:
        str(QC_FP/'decontam'/'intermediates'/'{host}'/'{sample}.bam')
    output:
        ids = str(QC_FP/'decontam'/'intermediates'/'{host}'/'{sample}.ids'),
    params:
        pct_id =  Cfg['qc']['pct_id'],
        frac = Cfg['qc']['frac']
    run:
        with open(output.ids, 'w') as out:
            last = None
            for read_id in get_mapped_reads(input[0], params.pct_id, params.frac):
                if read_id == last:
                    continue
                else:
                    out.write(read_id + '\n')
                    last = read_id

rule aggregate_reads:
    input:
        expand(
            str(QC_FP/'decontam'/'intermediates'/'{host}'/'{{sample}}.ids'),
            host=HostGenomes.keys())
    output:
        temp(str(QC_FP/'decontam'/'intermediates'/'{sample}_hostreads.ids')),
    run:
        if len(input) == 0:
            shell("touch {output}")
        else:
            shell("cat {input} > {output}")

rule filter_reads:
    input:
        hostreads = str(QC_FP/'decontam'/'intermediates'/'{sample}_hostreads.ids'),
        reads = str(QC_FP/'cleaned'/'{sample}_{rp}.fastq.gz'),
        hostids = expand(str(QC_FP/'decontam'/'intermediates'/'{host}'/'{{sample}}.ids'), host=HostGenomes.keys())
    output:
        reads = str(QC_FP/'decontam'/'{sample}_{rp}.fastq.gz'),
        log = str(QC_FP/'log'/'decontam'/'{sample}_{rp}.txt')
    run:
        original = int(str(subprocess.getoutput(
            "zcat {} | wc -l".format(input.reads))).strip())//4
        host = int(subprocess.getoutput(
            "cat {} | wc -l".format(input.hostreads)).strip())
        nonhost = int(original-host)
        shell("""
        gzip -dc {input.reads} | \
        rbt fastq-filter {input.hostreads} | \
        gzip > {output.reads}
        """)
        
        hostdict = OrderedDict()
        for hostid in input.hostids:
            hostname = os.path.basename(os.path.dirname(hostid))
            hostcts = int(subprocess.getoutput("cat {} | wc -l".format(hostid)).strip())
            hostdict[hostname] = hostcts
        
        with open(output.log, 'w') as log:
            log.write("{}\n".format("\t".join(list(hostdict.keys()) + ["host","nonhost"] )))
            log.write("{}\n".format("\t".join( map(str, list(hostdict.values()) + [host, nonhost]) )))
