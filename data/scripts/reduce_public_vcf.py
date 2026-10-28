"""Shrink the annotated 1000 Genomes VCF so that it fits in the GitHub repository.

All variant records and all 51 samples are kept. Only the fields that DoAggregate and the data-preparation scripts
use, plus a few common annotations, are kept:
  INFO   AC, AF, AN, DP, Func/Gene/ExonicFunc (refGene and ensGene), CADD_phred, REVEL_score, CLNSIG,
         gnomad40_AF, gnomad40_AF_nfe
  FORMAT GT, AD, DP, GQ (PL and the other per-sample fields are dropped)
Tool command lines (##GATKCommandLine, ##bcftools_*, ...) are removed from the header.

Usage:
    python3 reduce_public_vcf.py annovar.vcf.hg38_multianno.vcf annovar.vcf.hg38_multianno.reduced.vcf.gz
Optional --start BYTE --max-seconds S process the input in pieces; each piece is appended as a gzip member
(a concatenation of gzip members is a valid gzip file) and the next start offset is printed.
"""
import argparse
import gzip
import time

KEEP_INFO = ["AC", "AF", "AN", "DP", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "Func.ensGene",
             "Gene.ensGene", "ExonicFunc.ensGene", "CADD_phred", "REVEL_score", "CLNSIG", "gnomad40_AF",
             "gnomad40_AF_nfe"]
KEEP_FORMAT = ["GT", "AD", "DP", "GQ"]
KEEP_HEADER = ("##fileformat", "##FILTER", "##ALT", "##contig", "##reference", "##fileDate", "##phasing")


def reduce_header(line):
    if line.startswith("##INFO=<ID="):
        return line if line[11:].split(",", 1)[0] in KEEP_INFO else None
    if line.startswith("##FORMAT=<ID="):
        return line if line[13:].split(",", 1)[0] in KEEP_FORMAT else None
    if line.startswith("#CHROM") or line.startswith(KEEP_HEADER):
        return line
    return None


def reduce_record(line):
    f = line.rstrip("\n").split("\t")
    info = [kv for kv in f[7].split(";") if kv.split("=", 1)[0] in KEEP_INFO]
    f[7] = ";".join(info) if info else "."
    fmt = f[8].split(":")
    idx = [fmt.index(k) for k in KEEP_FORMAT if k in fmt]
    f[8] = ":".join(fmt[i] for i in idx)
    for s in range(9, len(f)):
        vals = f[s].split(":")
        f[s] = ":".join(vals[i] if i < len(vals) else "." for i in idx)
    return "\t".join(f) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf")
    ap.add_argument("out")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--max-seconds", type=float, default=float("inf"))
    a = ap.parse_args()
    t0 = time.time()
    mode = "wb" if a.start == 0 else "ab"
    with open(a.vcf, "rb") as fin, open(a.out, mode) as raw, gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9) as out:
        fin.seek(a.start)
        buf = []
        for bline in fin:
            line = bline.decode()
            if line.startswith("#"):
                r = reduce_header(line)
                if r:
                    buf.append(r)
            else:
                buf.append(reduce_record(line))
            if len(buf) >= 5000:
                out.write("".join(buf).encode())
                buf = []
                if time.time() - t0 > a.max_seconds:
                    print("NEXT", fin.tell())
                    return
        out.write("".join(buf).encode())
    print("DONE")


if __name__ == "__main__":
    main()
