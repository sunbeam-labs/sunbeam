from sunbeam.bfx.reports import (
    make_fastqc_report,
    parse_fastqc_quality,
)

def test_make_fastqc_report(tmp_path, DATA_DIR):
    input_report_fps = [
        DATA_DIR / "fastqc" / "SHORT_1_fastqc" / "fastqc_data.txt",
        DATA_DIR / "fastqc" / "SHORT_2_fastqc" / "fastqc_data.txt"
    ]
    output_report_fp = tmp_path / "fastqc_quality.tsv"
    log_fp = tmp_path / "log.txt"
    make_fastqc_report(input_report_fps, output_report_fp, log_fp)
    with open(DATA_DIR / "fastqc_quality.tsv") as f:
        expected = f.read()
    with open(output_report_fp) as f:
        assert f.read() == expected

def test_make_fastqc_report_nosamples(tmp_path):
    output_report_fp = tmp_path / "fastqc_quality.tsv"
    log_fp = tmp_path / "log.txt"
    make_fastqc_report([], output_report_fp, log_fp)
    with open(output_report_fp) as f:
        assert f.read() == "Samples\n"
    
        
def test_parse_fastqc_quality(DATA_DIR):
    input_fp = DATA_DIR / "fastqc" / "SHORT_1_fastqc" / "fastqc_data.txt"
    with open(input_fp) as f:
        observed = list(parse_fastqc_quality(f))
    assert observed[0] == ("1", "17.0")
    assert observed[1] == ("2", "17.0")
