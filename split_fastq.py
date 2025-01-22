from argparse import ArgumentParser
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import gzip

"""
在处理大文件时速度较慢
"""

class SplitFastq:
    
    def __init__(self, in_file: Path, out_file: Path, cutoff: int):
        self.in_file = in_file
        self.out_file = out_file
        self.cutoff = cutoff
        self.out_dir = self.out_file.parent

    def __call__(self, *args, **kwds):
        self.out_dir.mkdir(parents=True, exist_ok=True)
        with gzip.open(self.in_file, "rt") as fi, gzip.open(self.out_file, "wt") as fo:
    # 读取输入的FASTQ文件，处理每条记录
            for record in SeqIO.parse(fi, "fastq"):
                # 获取序列和质量值的前50个
                seq = record.seq[:50]
                qual = record.letter_annotations["phred_quality"][:50]
                
                # 清空 letter_annotations
                record.letter_annotations.clear()
                
                # 更新记录的序列和质量值
                record.seq = Seq(seq)
                record.letter_annotations["phred_quality"] = qual
                
                # 将更新后的记录写入输出文件
                SeqIO.write(record, fo, "fastq")


def main():
    ap = ArgumentParser()
    ap.add_argument("-i", "--in_file", type=str, required=True, help="输入文件")
    ap.add_argument("-o", "--out_file", type=str, required=True, help="输出文件")
    ap.add_argument("-c", "--cutoff", type=int, required=True, help="阈值")
    ag = ap.parse_args()
    SplitFastq(
        in_file=Path(ag.in_file).absolute(),
        out_file=Path(ag.out_file).absolute(),
        cutoff=ag.cutoff,
    )()


if __name__ == "__main__":
    main()
