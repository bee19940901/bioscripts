import gzip
from argparse import ArgumentParser
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt

class PlotFqLen:
    
    def __init__(self, in_fq_gz: Path, out_xls: Path, out_pdf: Path, out_png: Path):
        self.in_fq_gz = in_fq_gz
        self.out_xls = out_xls
        self.out_dir = out_xls.parent
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.out_pdf = out_pdf
        self.out_png = out_png
        
    def __call__(self, *args, **kwds):
        # 读取.gz格式的FASTQ文件
        with gzip.open(self.in_fq_gz, "rt") as handle:
            read_lengths = [len(record.seq) for record in SeqIO.parse(handle, "fastq")]
        
        # 绘制直方图
        plt.figure(figsize=(10, 6))
        plt.hist(read_lengths, bins=50, color="blue", alpha=0.7, edgecolor="black")

        plt.title("Distribution of Read Lengths")
        plt.xlabel("Read Length")
        plt.ylabel("Frequency")
        
        # 保存为PDF和PNG
        plt.savefig(self.out_pdf)
        plt.savefig(self.out_png)

        # 输出序列长度分布到Excel文件
        import pandas as pd
        length_df = pd.DataFrame(read_lengths, columns=["Length"])
        length_df.to_csv(self.out_xls, sep="\t", header=True, index=False)

# 示例使用：
if __name__ == "__main__":
    parser = ArgumentParser(description="Plot read length distribution from a .fq.gz file")
    parser.add_argument("-i", "--in_fq_gz", type=str, help="Input .fq.gz file")
    parser.add_argument("-o", "--out_xls", type=str, help="Output Excel file for lengths")
    parser.add_argument("-f", "--out_pdf", type=str, help="Output PDF file for plot")
    parser.add_argument("-g", "--out_png", type=str, help="Output PNG file for plot")
    args = parser.parse_args()
    PlotFqLen(
        Path(args.in_fq_gz).absolute, 
        Path(args.out_xls).absolute, 
        Path(args.out_pdf).absolute, 
        Path(args.out_png).absolute
    )()
    
