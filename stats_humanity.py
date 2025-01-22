import re
from argparse import ArgumentParser
from multiprocessing import Pool
from pathlib import Path
from subprocess import run
import pandas as pd


class StatsHumanity:

    def __init__(self, in_dir: Path, out_dir: Path):
        self.in_dir = in_dir
        self.out_dir = out_dir
        self.in_files = list(in_dir.glob("*_1.fq.gz"))
        self.samples = [re.sub(r"_1\.fq\.gz$", "", in_file.name) for in_file in self.in_files]
        self.out_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def filter_humanity(sample: str, in_dir: Path, out_dir: Path):
        """处理单个样本的过滤工作"""
        cmd = (
            f"bwa mem -t 32 /data/NGS/tngs_project/pipelineV2.5/database/hg19.fa "
            f"{in_dir.joinpath(f'{sample}_1.fq.gz')} "
            f"{in_dir.joinpath(f'{sample}_2.fq.gz')} | "
            f"samtools view -@ 16 -f 12 -h > "
            f"{out_dir.joinpath(f'{sample}.sam')} && "
            f"samtools fastq -1 {out_dir.joinpath(f'{sample}_1.unmapped.fq')} "
            f"-2 {out_dir.joinpath(f'{sample}_2.unmapped.fq')} "
            f"{out_dir.joinpath(f'{sample}.sam')} && "
            f"seqkit stats {in_dir.joinpath(f'{sample}_1.fq.gz')} "
            f"> {out_dir.joinpath(f'{sample}_1.stats.csv')} && "
            f"seqkit stats {in_dir.joinpath(f'{sample}_2.fq.gz')} "
            f"> {out_dir.joinpath(f'{sample}_2.stats.csv')} && "
            f"seqkit stats {out_dir.joinpath(f'{sample}_1.unmapped.fq')} "
            f"> {out_dir.joinpath(f'{sample}_1.unmapped.stats.csv')} && "
            f"seqkit stats {out_dir.joinpath(f'{sample}_2.unmapped.fq')} "
            f"> {out_dir.joinpath(f'{sample}_2.unmapped.stats.csv')} \n"
        )
        print(cmd)
        run(cmd, shell=True)

    def __multi_filter_humanity(self):
        """并行过滤原始数据中人源序列"""
        with Pool(4) as pool:
            for sample in self.samples:
                pool.apply_async(self.filter_humanity, args=(sample, self.in_dir, self.out_dir, ))
            pool.close()
            pool.join()
        return self

    @staticmethod
    def __stats_humanity(stats_file: Path):
        """计算单个文件的序列数"""
        return int(pd.read_csv(stats_file, sep=r"\s+", thousands=",").iloc[0, 3])

    def __multi_stats_humanity(self):
        """统计人源比例"""
        records = []
        for sample in self.samples:
            # 计算每个样本的比例
            humanity = (
                self.__stats_humanity(self.out_dir.joinpath(f"{sample}_1.unmapped.stats.csv")) +
                self.__stats_humanity(self.out_dir.joinpath(f"{sample}_2.unmapped.stats.csv"))
            ) / (
                self.__stats_humanity(self.out_dir.joinpath(f"{sample}_1.stats.csv")) +
                self.__stats_humanity(self.out_dir.joinpath(f"{sample}_2.stats.csv"))
            )

            records.append({"sample": sample, "humanity": f"{(1-humanity)*100:.2f}%"})

        # 将统计结果保存为 CSV 文件
        pd.DataFrame(records).to_csv(self.out_dir.joinpath("humanity.stats.csv"), sep="\t", header=True, index=False)
        return self

    def __call__(self, *args, **kwargs):
        """调用处理流程"""
        self.__multi_filter_humanity()
        self.__multi_stats_humanity()


if __name__ == "__main__":
    ap = ArgumentParser()
    ap.add_argument("-i", "--in_dir", type=str, required=True, help="输入目录")
    ap.add_argument("-o", "--out_dir", type=str, required=True, help="输出目录")
    ag = ap.parse_args()

    StatsHumanity(
        in_dir=Path(ag.in_dir).absolute(),
        out_dir=Path(ag.out_dir).absolute()
    )()
