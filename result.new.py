import re

from pathlib import Path
from argparse import ArgumentParser

import pandas as pd
import numpy as np


SCRIPT = Path(__file__).absolute()
SCRIPT_DIR = SCRIPT.parent
cyto19 = pd.read_csv(f'{SCRIPT_DIR}/../database/UCSC_hg19/cytoBand.txt', sep="\t", header=None)
chrL19 = pd.read_csv(f'{SCRIPT_DIR}/../database/UCSC_hg19/chr_length.txt', sep="\t")
cyto38 = pd.read_csv(f'{SCRIPT_DIR}/../database/hg38/cytoBand.txt', sep="\t", header=None)
chrL38 = pd.read_csv(f'{SCRIPT_DIR}/../database/hg38/chr_length.txt', sep="\t")
chrL38_dict = chrL38.set_index('chr')['length'].to_dict()
mosR = np.arange(0, 1.1, 0.1)
mosMdel = np.log2(2 - mosR) - 1
mosMdup = np.log2(2 + mosR) - 1
mosMdely = np.log2(1 - mosR) - 1
mosMdupy = np.log2(1 + mosR) - 1
names = [str(i) + "%" for i in range(0, 101, 10)]
mosMdel = pd.Series(dict(zip(names, mosMdel)))
mosMdup = pd.Series(dict(zip(names, mosMdup)))
mosMdely = pd.Series(dict(zip(names, mosMdely)))
mosMdely["100%"] = -10
mosMdupy = pd.Series(dict(zip(names, mosMdupy)))


def mosRatio(x: float) -> str:
    print(f"##### {x} ####")
    if x > mosMdup.max() or x < mosMdel.min():
        return "100%"
    else:
        if x >= 0:
            return mosMdup[mosMdup < x].index[-1]
        else:
            return mosMdel[mosMdel > x].index[-1]


def mosRatioy(x: float) -> str:
    if x > mosMdupy.max() or x < mosMdely.min():
        return "100%"
    else:
        if x > -1:
            return mosMdupy[mosMdupy < x].index[-1]
        else:
            return mosMdely[mosMdely > x].index[-1]


def outpd(x):
    d = pd.read_csv(x, sep="\t")
    med = d.iloc[:, 4].apply(lambda i: -11 if pd.isnull(i) else i)
    length = d.apply(lambda i: i[2] - i[1], axis=1)
    mos = med.apply(mosRatio)
    chrn = med.apply(lambda i: 2 ** (i + 1))
    mosy = med.apply(mosRatioy)
    return pd.DataFrame(
        {
            "Chr": d.iloc[:, 0],
            "Start": d.iloc[:, 1],
            "End": d.iloc[:, 2],
            "Length": length,
            "Median": med,
            "Mosiac": mos,
            "ChrNum": chrn,
            "Mosiacy": mosy,
        }
    )


def gender(d: pd.DataFrame) -> str:

    # 首先处理Y染色体
    m = d.loc[d["Chr"] == "chrY", "Median"]
    # 如果Y染色体存在分段
    if len(m) > 1:
        # 计算每一段的染色体条数
        n = 2 ** (m + 1)
        # 如果每一段的条数都小于0.3，则认为染色体不存在
        if (n < 0.3).all():
            genY = ""
        # 有一段大于0.3则报Y,因为分段后是不会报整体Y染色体重复的，所以不能报YY 或 YYY
        else:
            genY = "Y"
    # 如果Y染色体不分段，直接进行计算， 后面能报 YY
    elif len(m) == 1:
        n = 2 ** (m.values[0] + 1)
        # 条数小于0.3， 则认为不存在
        if n < 0.3:
            genY = ""
        # 条数在0.3到1.8之间，则认为一条
        elif 0.3 <= n < 1.8:
            genY = "Y"
        # 大于1.8则认为是多条
        else:
            genY = "Y" * round(n)
    # 完全不存在Y染色体
    else:
        genY = ""

    # 再处理X染色体
    m = d.loc[d["Chr"] == "chrX", "Median"]
    # 如果X染色体存在分段
    if len(m) > 1:
        # 计算每个分段的染色体条数
        n = 2 ** (m + 1)
        # 如果存在Y染色体，则不会报X染色体整体重复，不能报XXY(因为后面不会报X染色体整体重复，会冲突)，只能报XY，或者X缺失
        if genY:
            # 如果X每一段都小于0.3， 则认为X整体是缺失了的
            if (n < 0.3).all():
                genX = ""
            # 如果只要有一段X大于0.3 则认为X单条存在
            else:
                genX = "X"
        # 如果Y染色体不存在
        else:
            # 如果X每一段都小于0.3， 则认为X整体是缺失了的
            if (n < 0.3).all():
                genX = ""
            # 如果每一段的条数都小于1.8， 则认为是X，
            elif (n < 1.8).all():
                genX = "X"
            # 只要有一条大于1.8 则认为是XX(因为后面不会报X染色体整体重复，会冲突)，所以不能报两条以上，例如XXX 等
            else:
                genX = "XX"
    # 如果X染色体不分段， 后面就能报 "", X，XX, XXX, 等
    elif len(m) == 1:
        # 计算染色体条数
        n = 2 ** (m.values[0] + 1)
        # 条数小于0.3， 则认为不存在
        if n < 0.3:
            genX = ""
        # 条数在0.3到1.7之间，则认为一条
        elif 0.3 <= n <= 1.7:
            genX = "X"
        # 大于1.8则认为是多条
        else:
            genX = "X" * round(n)
    # 如果X染色体完全不存在，报 ""
    else:
        genX = ""

    if genY:
        gen = genX + genY
    elif genX == "X":
        gen = "XO"
    else:
        gen = genX
    return gen


def insert_com(x): return f"{x:,}"


def whoV_single(xx: pd.Series) -> str:
    pr = ""
    rf = ""
    r = ""
    cn = xx[0].replace("chr", "")
    if xx[7] in mosMdely.index[3:8]:
        if xx[4] > -1:
            pr = "+"
            rf = " (mos,x2,~"
        if xx[4] < -1:
            pr = "-"
            rf = " (mos,x0,~"
        r = pr + cn + rf + xx[7] + ")"
    if xx[7] in mosMdely.index[8:]:
        if xx[4] > -1:
            pr = "+"
        if xx[4] < -1:
            pr = "-"
        r = pr + cn
    return r


def whoV_double(xx: pd.Series) -> str:
    pr = ""
    rf = ""
    r = ""
    cn = xx[0].replace("chr", "")
    if xx[5] in mosMdel.index[3:8]:
        if xx[4] > 0:
            pr = "+"
            rf = " (mos,x3,~"
        if xx[4] < 0:
            pr = "-"
            rf = " (mos,x1,~"
        r = pr + cn + rf + xx[5] + ")"
    if xx[5] in mosMdel.index[8:]:
        if xx[4] > 0:
            pr = "+"
        if xx[4] < 0:
            pr = "-"
        if 2 ** (xx[4] + 1) > 3.5:
            rf = f" (x{round(2 ** (xx[4] + 1), 1)})"
        r = pr + cn + rf
    return r


def get_whoV(xx: pd.Series, g: str):
    r = ""
    cn = xx[0].replace("chr", "")
    if g == "XO":
        if cn == "X" or cn == "Y":
            pass
        else:
            r = whoV_double(xx)
    if "Y" in g:
        if cn == "X" or cn == "Y":
            r = whoV_single(xx)
        else:
            r = whoV_double(xx)
    if g == "XX" or g == "XXX" or g == "XXXX":
        if cn == "Y":
            pass
        else:
            r = whoV_double(xx)
    return r


def parV_double(xx: pd.Series, v: int) -> str:
    cyto = cyto38 if v == 38 else cyto19
    pr, rf, r = "", "", ""
    cn = xx[0].replace("chr", "")
    a = xx[1]
    b = xx[2]
    i = cyto.loc[cyto.iloc[:, 0] == xx[0], :]
    i.columns = ["chr", "start", "end", "site", "flag"]
    loc1 = i.loc[i["start"] < a, "site"].iloc[-1]
    if b >= i.iloc[-1, 2]:
        b = i.iloc[-1, 2]
        loc2 = i.iloc[-1, 3]
    else:
        loc2 = i.loc[i["end"] > b, "site"].iloc[0]
    startloc = insert_com(a)
    endloc = insert_com(b)
    if xx[5] in mosMdel.index[3:8]:
        if xx[4] > 0:
            pr = "dup"
            rf = " (mos,x3,~"
        else:
            pr = "del"
            rf = " (mos,x1,~"
        r = f"{pr}({cn})({loc1}{loc2}), chr{cn}:{startloc}-{endloc}{pr}{rf}{xx[5]})"
    elif xx[5] in mosMdel.index[8:]:
        if xx[4] > 0:
            pr = "dup"
        else:
            pr = "del"
        r = f"{pr}({cn})({loc1}{loc2}), chr{cn}:{startloc}-{endloc}{pr}"
    return r


def parV_single(xx: pd.Series, v: int) -> str:
    cyto = cyto38 if v == 38 else cyto19
    pr, rf, r = "", "", ""
    cn = xx[0].replace("chr", "")
    a = xx[1]
    b = xx[2]
    i = cyto.loc[cyto.iloc[:, 0] == xx[0], :]
    i.columns = ["chr", "start", "end", "site", "flag"]
    loc1 = i.loc[i["start"] < a, "site"].iloc[-1]
    if b >= i.iloc[-1, 2]:
        b = i.iloc[-1, 2]
        loc2 = i.iloc[-1, 3]
    else:
        loc2 = i.loc[i["end"] > b, "site"].iloc[0]
    startloc = insert_com(a)
    endloc = insert_com(b)
    if xx[7] in mosMdely.index[3:8]:
        if xx[4] > -1:
            pr = "dup"
            rf = " (mos,x2,~"
        else:
            pr = "del"
            rf = " (mos,x0,~"
        r = f"{pr}({cn})({loc1}{loc2}), chr{cn}:{startloc}-{endloc}{pr}{rf}{xx[7]})"
    if xx[7] in mosMdely.index[8:]:
        if xx[4] > -1:
            pr = "dup"
        else:
            pr = "del"
        r = f"{pr}({cn})({loc1}{loc2}), chr{cn}:{startloc}-{endloc}{pr}"
    return r


def get_parV(xx: pd.Series, g: str, v: int):
    r = ""
    cn = xx[0].replace("chr", "")
    if g == "XO":
        # 由于不知道是X缺失还是Y缺失，所以无法对性染色体进行处理
        if cn == "X" or cn == "Y":
            pass
        else:
            r = parV_double(xx, v)

    if "Y" in g:
        if cn == "X" or cn == "Y":
            r = parV_single(xx, v)
        else:
            r = parV_double(xx, v)

    if re.fullmatch(r"^X{2,}$", g):
        if cn == "Y":
            # 这种情况下已经判定Y染色体不存在，不会对Y染色体在再做任何处理
            pass
        else:
            r = parV_double(xx, v)
    return r


def chrV(d: pd.DataFrame, v: int) -> str:
    gen = gender(d)
    dp = None
    duplicate_values = d.loc[d.iloc[:, 0].duplicated(), "Chr"]
    # 有分段的染色体
    # 如果有性染色体进行分段，则性染色体数量按照性别中的来
    if len(duplicate_values) > 0:
        dw = d.loc[~d["Chr"].isin(duplicate_values.values), :]
        dp = d.loc[d["Chr"].isin(duplicate_values.values), :]
        # 获取分段的常染色体
        dn = dp.loc[~dp["Chr"].isin(["chrX", "chrY"]), :]
        # 如果chrX 分段 且 chrY分段
        if "chrX" in duplicate_values.values and "chrY" in duplicate_values.values:
            chrntotal = len(np.unique(dn.iloc[:, 0])) * 2 + \
                        np.sum(np.round(dw.iloc[:, 6])) + gen.count("X") + gen.count("Y")
        # 如果只有X染色体分段：
        elif "chrX" in duplicate_values.values and "chrY" not in duplicate_values.values:
            chrntotal = len(np.unique(dn.iloc[:, 0])) * 2 + np.sum(np.round(dw.iloc[:, 6])) + gen.count("X")
        # 如果只有Y染色体分段
        elif "chrY" in duplicate_values.values and "chrX" not in duplicate_values.values:
            chrntotal = len(np.unique(dn.iloc[:, 0])) * 2 + np.sum(np.round(dw.iloc[:, 6])) + gen.count("Y")
        # 如果只有常染色体分段
        else:
            chrntotal = len(np.unique(dn.iloc[:, 0])) * 2 + np.sum(np.round(dw.iloc[:, 6]))
    else:
        dw = d
        chrntotal = np.sum(np.round(dw.iloc[:, 6]))
    chrntotal = int(chrntotal)
    whoV = dw.apply(get_whoV, args=(gen, ), axis=1)
    whoV = ", ".join(i for i in whoV.tolist() if i)

    if dp is None:
        res = ", ".join(str(i) for i in [chrntotal, gen, whoV] if i)
    else:
        parV = dp.apply(get_parV, args=(gen, v, ), axis=1)
        parV = ", ".join(i for i in parV.tolist() if i)
        res = ", ".join(str(i) for i in [chrntotal, gen, whoV, parV] if i)
    res = re.sub(r"[\s,]+$", "", res)
    return res


def whoV_gz_single(xx: pd.Series) -> str:

    r = ""
    cn = xx[0].replace("chr", "")

    if xx[7] in mosMdely.index[3:8]:

        if xx[4] > -1:
            pr = "+"
            rf = " (x2,mos,~"
            r = pr + cn + rf + xx[7] + ")"

        if xx[4] < -1:
            pr = "-"
            rf = " (x0,mos,~"
            r = pr + cn + rf + xx[7] + ")"

    if xx[7] in mosMdely.index[8:]:

        if xx[4] > -1:
            pr = "+"
            rf = " (x2)"
            r = pr + cn + rf
        
        if xx[4] > -1 and 2 ** (xx[4] + 1) > 2.45:
            pr = "+"
            rf = f" (x{round(2 ** (xx[4] + 1), 1)})"
        
        if xx[4] < -1:
            pr = "-"
            rf = " (x0)"
            r = pr + cn + rf

    return r


def whoV_gz_double(xx: pd.Series) -> str:
    pr = ""
    rf = ""
    r = ""
    cn = xx[0].replace("chr", "")

    if xx[5] in mosMdel.index[3:8]:
        if xx[4] > 0:
            pr = "+"
            rf = " (x3,mos,~"
        else:
            pr = "-"
            rf = " (x1,mos,~"
        r = pr + cn + rf + xx[5] + ")"

    if xx[5] in mosMdel.index[8:]:
        if xx[4] > 0:
            pr = "+"
            rf = " (x3)"
        if xx[4] > 0 and 2 ** (xx[4] + 1) > 3.45:
            pr = "+"
            rf = f" (x{round(2 ** (xx[4] + 1), 1)})"
        if xx[4] < 0:
            pr = "-"
            rf = " (x1)"
        r = pr + cn + rf
    return r


def get_whoV_gz(xx: pd.Series, g: str) -> str:

    r = ""
    cn = xx[0].replace("chr", "")

    if g == "XO":

        if cn == "X" or cn == "Y":
            pass

        else:
            r = whoV_gz_double(xx)

    if "Y" in g:

        if cn == "X" or cn == "Y":
            r = whoV_gz_single(xx)

        else:
            r = whoV_gz_double(xx)

    # 如果出现 XX XXX XXXX 等性别 按照 XX 处理
    if re.fullmatch(r"^X{2,}$", g):

        if cn == "Y":
            pass

        else:
            r = whoV_gz_double(xx)

    return r


def parV_gz_single(xx: pd.Series, v: int) -> str:
    cyto = cyto38 if v == 38 else cyto19
    pr = ""
    rf = ""
    r = ""
    cn = xx[0].replace("chr", "")
    a = xx[1]
    b = xx[2]
    i = cyto.loc[cyto.iloc[:, 0] == xx[0], :]
    i.columns = ["chr", "start", "end", "site", "flag"]
    loc1 = i.loc[i["start"] < a, "site"].iloc[-1]
    if b >= i.iloc[-1, 2]:
        b = i.iloc[-1, 2]
        loc2 = i.iloc[-1, 3]
    else:
        loc2 = i.loc[i["end"] > b, "site"].iloc[0]
    if loc1[0] == loc2[0]:
        loc = loc1[0]
    else:
        loc = f"{loc1[0]}{loc2[0]}"
    length = round((b - a) / 1000000, 1)
    if xx[7] in mosMdely.index[3:8]:
        if xx[4] > -1:
            pr = "+"
            rf = ",x2,mos,~"
        else:
            pr = "-"
            rf = ",x0,mos,~"
        r = f"{pr}{cn}{loc}({loc1}→{loc2},~{length}Mb{rf}{xx[7]})"
    elif xx[7] in mosMdely.index[8:]:
        if xx[4] > -1:
            pr = "+"
            rf = ",x2"
        if xx[4] > -1 and 2 ** (xx[4] + 1) > 2.45:
            pr = "+"
            rf = f",x{round(2 ** (xx[4] + 1), 1)}"
        if xx[4] < -1:
            pr = "-"
            rf = ",x0"
        r = f"{pr}{cn}{loc}({loc1}→{loc2},~{length}Mb{rf})"
    return r


def parV_gz_double(xx: pd.Series, v: int) -> str:
    cyto = cyto38 if v == 38 else cyto19
    pr = ""
    rf = ""
    r = ""
    cn = xx[0].replace("chr", "")
    a = xx[1]
    b = xx[2]
    i = cyto.loc[cyto.iloc[:, 0] == xx[0], :]
    i.columns = ["chr", "start", "end", "site", "flag"]
    loc1 = i.loc[i["start"] < a, "site"].iloc[-1]
    if b >= i.iloc[-1, 2]:
        b = i.iloc[-1, 2]
        loc2 = i.iloc[-1, 3]
    else:
        loc2 = i.loc[i["end"] > b, "site"].iloc[0]
    if loc1[0] == loc2[0]:
        loc = loc1[0]
    else:
        loc = f"{loc1[0]}{loc2[0]}"
    length = round((b - a) / 1000000, 1)
    if xx[5] in mosMdel.index[3:8]:
        if xx[4] > 0:
            pr = "+"
            rf = ",x3,mos,~"
        else:
            pr = "-"
            rf = ",x1,mos,~"
        r = f"{pr}{cn}{loc}({loc1}→{loc2},~{length}Mb{rf}{xx[5]})"

    elif xx[5] in mosMdel.index[8:]:
        if xx[4] > 0:
            pr = "+"
            rf = ",x3"
        if xx[4] > 0 and 2 ** (xx[4] + 1) > 3.45:
            pr = "+"
            rf = f",x{round(2 ** (xx[4] + 1), 1)}"
        if xx[4] < 0:
            pr = "-"
            rf = ",x1"
        r = f"{pr}{cn}{loc}({loc1}→{loc2},~{length}Mb{rf})"
    return r


def get_parV_gz(xx: pd.Series, g: str, v: int) -> str:

    r = ""
    cn = xx[0].replace("chr", "")

    if g == "XO":

        if cn == "X" or cn == "Y":
            pass

        else:
            r = parV_gz_double(xx, v)

    # 如果出现Y染色体，则按照XY性别处理
    if "Y" in g:

        if cn == "X" or cn == "Y":
            r = parV_gz_single(xx, v)

        else:
            r = parV_gz_double(xx, v)

    # 如果只有X染色体，且X重复再两次及以上，按照XX处理
    if re.fullmatch(r"^X{2,}$", g):

        if cn == "Y":
            pass

        else:
            r = parV_gz_double(xx, v)

    return r


def chrV_gz(d: pd.DataFrame, v: int):
    gen = gender(d)
    dp = None
    duplicate_values = d.loc[d.iloc[:, 0].duplicated(), "Chr"]
    if len(duplicate_values.values) > 0:
        dw = d.loc[~d["Chr"].isin(duplicate_values.values), :]
        dp = d.loc[d["Chr"].isin(duplicate_values.values), :]
        # 获取分段的常染色体
        dn = dp.loc[~dp["Chr"].isin(["chrX", "chrY"]), :]
        # 如果chrX 分段 且 chrY分段
        if "chrX" in duplicate_values.values and "chrY" in duplicate_values.values:
            chrntotal = len(np.unique(dn.iloc[:, 0])) * 2 + np.sum(np.round(dw.iloc[:, 6])) \
                        + gen.count("X") + gen.count("Y")
        # 如果只有X染色体分段：
        elif "chrX" in duplicate_values.values:
            chrntotal = len(np.unique(dn.iloc[:, 0])) * 2 + np.sum(np.round(dw.iloc[:, 6])) + gen.count("X")
        # 如果只有Y染色体分段
        elif "chrY" in duplicate_values.values:
            chrntotal = len(np.unique(dn.iloc[:, 0])) * 2 + np.sum(np.round(dw.iloc[:, 6])) + gen.count("Y")
        # 如果只有常染色体
        else:
            chrntotal = len(np.unique(dn.iloc[:, 0])) * 2 + np.sum(np.round(dw.iloc[:, 6]))
    else:
        dw = d
        chrntotal = np.sum(np.round(dw.iloc[:, 6]))
    chrntotal = int(chrntotal)
    whoV = dw.apply(get_whoV_gz, args=(gen, ), axis=1)
    whoV = ", ".join(i for i in whoV.tolist() if i)
    if dp is None:
        res = ", ".join(str(i) for i in [chrntotal, gen, whoV] if i)
    else:
        parV = dp.apply(get_parV_gz, args=(gen, v, ), axis=1)
        parV = ", ".join(i for i in parV.tolist() if i)
        res = ", ".join(str(i) for i in [chrntotal, gen, whoV, parV] if i)
    res = re.sub(r"[\s,]+$", "", res)
    return res


def whoVN_double(xx: pd.Series) -> str:
    pr = ""
    r = ""
    rf = ""
    cn = xx[0].replace("chr", "")
    if xx[5] in mosMdel.index[3:8]:
        if xx[4] > 0:
            pr = "+"
            rf = " (mos,x3,~"
        if xx[4] < 0:
            pr = "-"
            rf = " (mos,x1,~"
        r = pr + cn + rf + xx[5] + ")"
    if xx[5] in mosMdel.index[8:]:
        if xx[4] > 0:
            pr = "+"
        if xx[4] < 0:
            pr = "-"
        r = pr + cn
    return r


def whoVN_single(xx: pd.Series) -> str:
    pr = ""
    r = ""
    rf = ""
    cn = xx[0].replace("chr", "")
    if xx[7] in mosMdely.index[3:8]:
        if xx[4] > -1:
            pr = "+"
            rf = " (mos,x3,~"
        if xx[4] < -1:
            pr = "-"
            rf = " (mos,x1,~"
        r = pr + cn + rf + xx[7] + ")"
    if xx[7] in mosMdely.index[8:]:
        if xx[4] > -1:
            pr = "+"
        if xx[4] < -1:
            pr = "-"
        r = pr + cn
    return r


def get_whoVN(xx: pd.Series, g: str) -> str:
    r = ""
    cn = xx[0].replace("chr", "")
    if g == "XO":
        if cn == "X" or cn == "Y":
            pass
        else:
            r = whoVN_double(xx)
    if "Y" in g:
        if cn == "X" or cn == "Y":
            r = whoVN_single(xx)
        else:
            r = whoVN_double(xx)
    if re.fullmatch(r"^X{2,}$", g):
        if cn == "Y":
            pass
        else:
            r = whoVN_double(xx)
    return r


def parVN_double(xx: pd.Series, v: int) -> str:
    cyto = cyto38 if v == 38 else cyto19
    chrL = chrL38 if v == 38 else chrL19
    pr, rf, r = "", "", ""
    cn = xx[0].replace("chr", "")
    a = xx[1]
    b = xx[2]
    ou = True
    ccc = b-a
    if b < chrL[xx[0], 2] and ccc < chrL[xx[0], 5]:
        ou = False
    if a > chrL[xx[0], 2] and ccc < chrL[xx[0], 6]:
        ou = False
    if a < chrL[xx[0], 2] < b and ccc < chrL[xx[0], 4]:
        ou = False
    i = cyto.loc[cyto.iloc[:, 0] == xx[0], :]
    loc1 = i.loc[i["start"] < a, "site"].iloc[-1]
    if b >= i.iloc[-1, 2]:
        b = i.iloc[-1, 2]
        loc2 = i.iloc[-1, 3]
    else:
        loc2 = i.loc[i["end"] > b, "site"].iloc[0]
    startloc = insert_com(a)
    endloc = insert_com(b)
    if xx[5] in mosMdel.index[3:8] and ou:
        if xx[4] > 0:
            pr = "dup"
            rf = " (mos,x3,~"
        if xx[4] < 0:
            pr = "del"
            rf = " (mos,x1,~"
        r = f"{pr}({cn})({loc1}{loc2}), chr{cn}:{startloc}-{endloc}{pr}{rf}{xx[5]})"
    elif xx[5] in mosMdel.index[8:] and ou:
        if xx[4] > 0:
            pr = "dup"
        if xx[4] < 0:
            pr = "del"
        r = f"{pr}({cn})({loc1}{loc2}), chr{cn}:{startloc}-{endloc}{pr}"
    return r


def parVN_single(xx: pd.Series, v: int) -> str:
    cyto = cyto38 if v == 38 else cyto19
    chrL = chrL38 if v == 38 else chrL19
    pr, rf, r = "", "", ""
    cn = xx[0].replace("chr", "")
    a = xx[1]
    b = xx[2]
    ou = True
    ccc = b-a
    if b < chrL[xx[0], 2] and ccc < chrL[xx[0], 5]:
        ou = False
    if a > chrL[xx[0], 2] and ccc < chrL[xx[0], 6]:
        ou = False
    if a < chrL[xx[0], 2] < b and ccc < chrL[xx[0], 4]:
        ou = False
    i = cyto.loc[cyto.iloc[:, 0] == xx[0], :]
    loc1 = i.loc[i["start"] < a, "site"].iloc[-1]
    if b >= i.iloc[-1, 2]:
        b = i.iloc[-1, 2]
        loc2 = i.iloc[-1, 3]
    else:
        loc2 = i.loc[i["end"] > b, "site"].iloc[0]
    startloc = insert_com(a)
    endloc = insert_com(b)
    if xx[7] in mosMdely.index[3:8] and ou:
        if xx[4] > -1:
            pr = "dup"
            rf = " (mos,x2,~"
        if xx[4] < -1:
            pr = "del"
            rf = " (mos,x0,~"
        r = f"{pr}({cn})({loc1}{loc2}), chr{cn}:{startloc}-{endloc}{pr}{rf}{xx[7]})"
    elif xx[7] in mosMdely.index[8:] and ou:
        if xx[4] > -1:
            pr = "dup"
        if xx[4] < -1:
            pr = "del"
        r = f"{pr}({cn})({loc1}{loc2}), chr{cn}:{startloc}-{endloc}{pr}"
    return r


def get_parVN(xx: pd.Series, g: str, v: int) -> str:
    cn = xx[0].replace("chr", "")
    r = ""
    if g == "XO":
        if cn == "X" or cn == "Y":
            pass
        else:
            r = parVN_double(xx, v)
    if "Y" in g:
        if cn == "X" or cn == "Y":
            r = parVN_single(xx, v)
        else:
            r = parVN_double(xx, v)
    if re.fullmatch(r"^X{2,}$", g):
        if cn == "Y":
            pass
        else:
            r = parVN_double(xx, v)
    return r


def chrVN(d: pd.DataFrame, v: int) -> str:
    gen = gender(d)
    dp = None
    duplicate_values = d.loc[d.iloc[:, 0].duplicated(), "Chr"]
    if len(duplicate_values.values) > 0:
        dw = d.loc[~d["Chr"].isin(duplicate_values.values), :]
        dp = d.loc[d["Chr"].isin(duplicate_values.values), :]
        # 获取分段的常染色体
        dn = dp.loc[~dp["Chr"].isin(["chrX", "chrY"]), :]
        # 如果chrX 分段 且 chrY分段
        if "chrX" in duplicate_values.values and "chrY" in duplicate_values.values:
            chrntotal = len(np.unique(dn.iloc[:, 0])) * 2 + np.sum(np.round(dw.iloc[:, 6])) + gen.count(
                "X") + gen.count("Y")
        # 如果只有X染色体分段：
        elif "chrX" in duplicate_values.values:
            chrntotal = len(np.unique(dn.iloc[:, 0])) * 2 + np.sum(np.round(dw.iloc[:, 6])) + gen.count("X")
        # 如果只有Y染色体分段
        elif "chrY" in duplicate_values.values:
            chrntotal = len(np.unique(dn.iloc[:, 0])) * 2 + np.sum(np.round(dw.iloc[:, 6])) + gen.count("Y")
        # 如果只有常染色体
        else:
            chrntotal = len(np.unique(dn.iloc[:, 0])) * 2 + np.sum(np.round(dw.iloc[:, 6]))
    else:
        dw = d
        chrntotal = np.sum(np.round(dw.iloc[:, 6]))
    chrntotal = int(chrntotal)
    whoV = dw.apply(get_whoVN, args=(gen, ), axis=1)
    whoV = ", ".join(i for i in whoV.tolist() if i)
    if dp is None:
        res = ", ".join(str(i) for i in [chrntotal, gen, whoV] if i)
    else:
        parV = dp.apply(get_parVN, args=(gen, v, ), axis=1)
        parV = ", ".join(i for i in parV.tolist() if i)
        res = ", ".join(str(i) for i in [chrntotal, gen, whoV, parV] if i)
    res = re.sub(r"[\s,]+$", "", res)
    return res

def func1(s):
    c = s["Chr"]
    if int(chrL38_dict[c]) - int(s["End"]) < 500000:
        return int(chrL38_dict[c])
    else:
        return int(s["End"])

def func(x, sample, version):
    data = outpd(x)
    data["End"] = data.apply(func1, axis=1) 
    r1 = chrV(data, version)
    r2 = chrV_gz(data, version)
    r = pd.Series([sample, r1, r2])
    return r
    

def check_chr_number(x: str) -> str:
    """
    根据结果对染色体总数进行修正
    :param x: 判读结果
    :return: 修正后的染色体总数, 删除 +Y/X 和 -X/Y 的形式
    """
    print(f"#### 修正前 #### {x}")
    items = x.split(", ")
    gen = items[1]
    num = gen.count("X") + gen.count("Y")
    new_items = list()
    n = 0
    if len(items) > 2:
        for i in items[2:]:
            if re.search(r"^\+\d+(?:\s|$)", i):
                n += 1
                new_items.append(i)
            elif re.search(r"^-\d+(?:\s|$)", i):
                n -= 1
                new_items.append(i)
            # X 或 Y 染色体 整体缺失 舍弃
            elif re.search(r"^[+-][XY](?:\s\(x\d+\))?$", i):
                continue
            else:
                new_items.append(i)
        if new_items:
            print(f"#### 修正后 #### {44 + n + num}, {gen}, {', '.join(new_items)}")
            return f"{44 + n + num}, {gen}, {', '.join(new_items)}"
        else:
            print(f"#### 修正后 #### {44 + num}, {gen}")
            return f"{44 + num}, {gen}"
    else:
        print(f"#### 修正后 #### {44 + num}, {gen}")
        return f"{44 + num}, {gen}"
        
def main():

    ap = ArgumentParser()
    ap.add_argument("-i", "--filepath", type=str, default="./inputs", help="input a dir [required]")
    ap.add_argument("-o", "--output", type=str, default="./outputs_py", help="output path")
    ap.add_argument("-c", "--info", type=str, required=True)
    args = ap.parse_args()

    in_dir = Path(args.filepath).absolute()
    out_dir = Path(args.output).absolute()
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(args.info, sheet_name=0)

    # 样本名对应的要使用的基因组版本
    df["version"] = df["test_project"].apply(lambda x: 19 if re.search(r"^NBS", x) else 38)
    sample_version = df.set_index('sampleid')['version'].to_dict()

    # 所有的500K文件
    in_500K_files = list(in_dir.glob("*.500k.wig.seg"))
    in_500K_samples = [
        re.sub(r"\.500k\.wig\.seg$", "", i.name)
        for i in in_500K_files
    ]
    r = {
        s: func(x=f, sample=s, version=sample_version[s])
        for s, f
        in zip(in_500K_samples, in_500K_files)
    }
    df_500 = pd.DataFrame(r).T
    df_500.columns = ["sampleId", "result", "resultgz"]
    df_500["result"] = df_500["result"].apply(check_chr_number)
    df_500["resultgz"] = df_500["resultgz"].apply(check_chr_number)
    
    # 所有的300K文件
    in_300K_files = list(in_dir.glob("*.300k.wig.seg"))
    in_300K_samples = [
        re.sub(r"\.300k\.wig\.seg$", "", i.name)
        for i in in_300K_files
    ]
    r = {
        s: func(x=f, sample=s, version=sample_version[s])
        for s, f
        in zip(in_300K_samples, in_300K_files)
    }
    df_300 = pd.DataFrame(r).T
    df_300.columns = ["sampleId", "result-300k", "resultgz-300k"]
    df_300["result-300k"] = df_300["result-300k"].apply(check_chr_number)
    df_300["resultgz-300k"] = df_300["resultgz-300k"].apply(check_chr_number)
    
    df = pd.merge(left=df_500, right=df_300, on="sampleId", how="inner")
    for i in ["remark", "explanation", "suggest", "gender"]:
        df[i] = ""
    
    df.to_csv(out_dir.joinpath("pre-result-merge.txt"), sep="\t", header=True, index=False)
    
    with open(out_dir.joinpath("pre-result.txt"), "w") as fw:
        fw.write("sampleId\tresult\tresultgz\n")
        for _, row in df_500.iterrows():
            fw.write(f"{row['sampleId']}\t{row['result']}\t{row['resultgz']}\n")
    with open(out_dir.joinpath("pre-result-300k.txt"), "w") as fw:
        fw.write("sampleId\tresult-300k\tresultgz-300k\n")
        for _, row in df_300.iterrows():
            fw.write(f"{row['sampleId']}\t{row['result-300k']}\t{row['resultgz-300k']}\n")


if __name__ == "__main__":
    main()
