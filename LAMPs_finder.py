import re
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
from Bio import SeqIO

JJBJ = {
    ("a",): "a",
    ("t",): "t",
    ("g",): "g",
    ("c",): "c",
    ("a", "g"): "r",
    ("c", "t"): "y",
    ("c", "g"): "s",
    ("a", "t"): "w",
    ("g", "t"): "k",
    ("a", "c"): "m",
    ("c", "g", "t"): "b",
    ("a", "g", "t"): "d",
    ("a", "c", "t"): "h",
    ("a", "c", "g"): "v",
    ("a", "c", "g", "t"): "n"
}


class LAMPsFinder:
    """
    根据多序列比结果寻找LAMP引物
    """

    def __init__(self, in_fas: Path, out_fas: Path, window: int, step: int, min_len: int):
        self.__out_fas = out_fas
        self.__in_fas = in_fas
        self.__window = window
        self.__step = step
        self.__min_len = min_len

    @staticmethod
    def lamps_finder(seg: str, min_len: int = 16) -> set[str]:
        """
        寻找器
        :param seg: 序列区间
        :param min_len: 最小片段长度
        :return:
        """
        results = []
        # 拆分片段
        flg = [_ for _ in re.split(r"[^atcg]{2,}", seg) if _]
        # 处理每一个片段
        for f in flg:
            # 寻找断点
            wds = [i for i, s in enumerate(f) if s not in "atcg"]
            # 如果断点数量大于1
            if len(wds) > 1:
                for n, _ in enumerate(wds):
                    # 第一个断点
                    if n == 0:
                        results.append(f[:wds[n + 1]])
                    # 最后一个断点
                    elif n == len(wds) - 1:
                        results.append(f[wds[n - 1] + 1:])
                    # 中间的断点
                    else:
                        results.append(f[wds[n - 1] + 1:wds[n + 1]])
            else:
                # 如果断点数量小于等于1直接保留
                results.append(f)
        outs = set()
        # 对结果进行修剪
        for result in results:
            if len(result) < min_len:
                continue
            else:
                if re.search(r"[^atcg]", result):
                    left, right = re.split(r"[^atcg]", result)
                    if len(left) < min_len and len(right) < min_len:
                        outs.add(result)
                    elif len(left) < min_len <= len(right):
                        outs.add(right)
                    elif len(left) >= min_len > len(right):
                        outs.add(left)
                    else:
                        outs.add(right)
                        outs.add(left)
                else:
                    outs.add(result)
        return outs

    def __call__(self, *args, **kwargs):
        self.__out_fas.parent.mkdir(exist_ok=True, parents=True)
        fas_dict = {}
        with open(self.__in_fas, "r", encoding="utf-8") as fr:
            for record in SeqIO.parse(fr, "fasta"):
                fas_dict[record.id] = list(record.seq)
        fas_df = pd.DataFrame(fas_dict)
        ref = "".join(fas_df.apply(lambda x: JJBJ.get(tuple(sorted(list(set(x.tolist())))), "n"), axis=1).tolist())
        outputs = {}
        for i in range(0, len(ref) - self.__window, self.__step):
            outputs[(i, i + self.__window)] = self.lamps_finder(ref[i: i + self.__window])
        outputs[(len(ref) - self.__window, len(ref))] = self.lamps_finder(ref[len(ref) - self.__window:len(ref)])
        seqs = []
        with open(self.__out_fas, "w") as fw:
            for k, v in outputs.items():
                if len(v) >= 6:
                    if v not in seqs:
                        seqs.append(v)
                        fw.write(
                            f">ALIGN: {k[0]}-{k[1]} "
                            f"SEGMENTS: {','.join(v)}\n{ref[k[0]:k[1]].upper()}\n"
                        )


def main():
    ap = ArgumentParser()
    ap.add_argument("-i", "--in_fas", type=str, required=True, help="多序列比对结果，FASTA格式")
    ap.add_argument("-o", "--out_fas", type=str, required=True, help="相对保守区域，FASTA格式")
    ap.add_argument("-w", "--window", type=int, required=False, default=400, help="滑窗大小，默认400bp")
    ap.add_argument("-s", "--step", type=int, required=False, default=20, help="步长，默认20bp")
    ap.add_argument("-m", "--min_len", type=int, required=False, default=16, help="最小引物长度，16bp")
    ag = ap.parse_args()
    LAMPsFinder(
        in_fas=Path(ag.in_fas).absolute(),
        out_fas=Path(ag.out_fas).absolute(),
        window=int(ag.window),
        step=int(ag.step),
        min_len=int(ag.min_len)
    )()


if __name__ == "__main__":
    main()
