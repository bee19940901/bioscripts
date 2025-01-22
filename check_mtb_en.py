from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
import sys
import pymysql


# 数据库连接配置
def create_connection():
    return pymysql.connect(
        host="127.0.0.1",  # 或 'localhost'
        port=3306,
        user="tngs",
        password='Tngs@123654',  # 确保密码正确
        database="emngs",
        charset='utf8mb4',
        cursorclass=pymysql.cursors.DictCursor
    )


# 常量定义

# 上机情况表 流程类型
SAMPLE_PROJECTS = [
    "General_tNGS",
    "mNGS",
    "Ophth_tNGS",
    "Myco_tNGS",
    "UTI_tNGS",
    "STI_tNGS",
    "Pediatric_tNGS"
]

# 上机情况表 样本信息表 送检实验室
LABS = ["Local Lab"]

# 上机情况表 样本呢类型
SAMPLE_TYPES = [
    "Clinical_Sample",
    "No_Template_Control"
]

# 样本信息表 模板类型
TEMPLATE_TYPES = [
    "Standard Template",
    "Ophth Template"
]

# 样本信息表 检测项目
TEST_PROJECTS = [
    "T-NGS-200",
    "Ophth_tNGS",
    "Myco_tNGS"
]


class CheckTngs:

    def __init__(self, fc_id: str):
        self.fc_id = fc_id
        self.out_dir = Path("/data/en_NGS/mtb_project/runmtb/result")
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.out_file = self.out_dir.joinpath(f"{self.fc_id}.err")
        self.out_info = str()
        self.csv_df = pd.DataFrame()
        self.lims_df = pd.DataFrame()

    def write_output(self, message: str):
        """写入输出文件并退出程序"""
        with open(self.out_file, "w", encoding="utf-8") as fw:
            fw.write(message)
        sys.exit(1)

    def query_csv(self):
        """从数据库查询上机情况表"""
        with create_connection() as connection:
            with connection.cursor() as cursor:
                sql = (
                    "SELECT laboratory, date, sampleid, indexid, indexseq, sample_project, description "
                    "FROM tngs_run_job_sample "
                    "WHERE fcid=%s "
                )
                cursor.execute(sql, (self.fc_id,))
                csv_records = cursor.fetchall()
        if not csv_records:
            self.write_output(f"Can't access run information of chip {self.fc_id}!")
        self.csv_df = pd.DataFrame(csv_records)
        print(self.csv_df)
        return self

    def query_lims(self):
        """从数据库查询样本信息表"""
        sampleids = self.csv_df['sampleid'].tolist()
        if not sampleids:
            self.write_output(f"Can't access sample information table of chip {self.fc_id}.")
        placeholders = ', '.join(['%s'] * len(sampleids))
        sql = (
            f"SELECT lab_name, project_type, reportid, sampleid, hospitalid, name, "
            f"sex, age, hospitalseq, bed, phone, hospital, department, doctor, "
            f"collectdate, receivedate, testdate, reportdate, reporttemp, sampletype, "
            f"clinical_diagnosis, target_pathogen, note "
            f"FROM tngs_tpseq_sample_info "
            f"WHERE sampleid IN ({placeholders});"
        )
        with create_connection() as connection:
            with connection.cursor() as cursor:
                cursor.execute(sql, sampleids)
                lims_records = cursor.fetchall()
        if not lims_records:
            self.write_output(f"Can't access sample information table of chip {self.fc_id}.")
        self.lims_df = pd.DataFrame(lims_records)
        print(self.lims_df)
        return self

    def check_lims(self):
        """检查样本信息表中的信息是否正确"""
        for _, row in self.lims_df.iterrows():
            # 检查实验室
            if not row["lab_name"]:
                self.out_info += f"The Lab of sample {row['sampleid']} is missing, please add it.<br>"
            elif row["lab_name"] not in LABS:
                self.out_info += f"The Lab of sample {row['sampleid']} is incorrect, please modify it.<br>"
            # 检查检测项目
            if not row["project_type"]:
                self.out_info += f"The Test Type of sample {row['sampleid']} is missing, please add it.<br>"
            elif row["project_type"] not in TEST_PROJECTS:
                self.out_info += f"The Test Type of sample {row['sample']} is incorrect, please modify it.<br>"
            # 检查报告编号
            if not row["reportid"]:
                self.out_info += f"The report ID of sample {row['sampleid']} is missing, please provide it.<br>"
            # 检查报告日期
            if not row["reportdate"]:
                self.out_info += f"The report date of sample {row['sampleid']} is missing, please provide it.<br>"
            # 检查模板类型
            if not row["reporttemp"]:
                self.out_info += f"The report template of sample {row['sampleid']} is missing, please provide it.<br>"
            elif row["reporttemp"] not in TEMPLATE_TYPES:
                self.out_info += f"The report template of sample {row['sampleid']} if incorrect, please modify it.<br>"
            # 检查样类型
            if not row["sampletype"]:
                self.out_info += f"The type of sample {row['sampleid']} id missing, please provide it.<br>"
        return self

    def check_csv(self):
        """检查上机情况表中的信息是否正确"""
        for _, row in self.csv_df.iterrows():
            # 检查实验室
            if not row["laboratory"]:
                self.out_info += f"The Laboratory of sample {row['sampleid']} " \
                                 f"in run information table is missing, please provide it.<br>"
            elif row["laboratory"] not in LABS:
                self.out_info += f"The laboratory of sample {row['sampleid']} " \
                                 f"in run information table is incorrect, please modify it.<br>"
            # 检查上机日期
            if not row["date"]:
                self.out_info += f"The test date of sample {row['sampleid']} is missing, please provide it.<br>"
            # 检查样本编号
            if not row["sampleid"]:
                self.out_info += f"The chip {self.fc_id} has an unknown sample id，please modify it.<br>"
            # 检查index编号
            if not row["indexid"]:
                self.out_info += f"The index id of sample {row['sampleid']} is missing, please provide it.<br>"
            # 检查流程
            if not row["sample_project"]:
                self.out_info += f"The project type of sample {row['sampleid']} is missing, please provide it.<br>"
            elif row["sample_project"] not in SAMPLE_PROJECTS:
                self.out_info += f"The project type of sample {row['sampleid']} is incorrect, please modify it.<br>"
            # 检查样本类型
            if not row["description"]:
                self.out_info += f"The type of sample {row['sampleid']} is missing, please provide it.<br>"
            elif row["description"] not in SAMPLE_TYPES:
                self.out_info += f"The type of sample {row['sampleid']} is incorrect, please modify it.<br>"
        return self

    def check_csv_lims(self):
        """检查上机情况表和样本信息表对应关系"""
        samples = self.lims_df['sampleid'].to_list()
        for _, row in self.csv_df.iterrows():
            if row['description'] == "Clinical_Sample" and row["sampleid"] not in samples:
                self.out_info += f"The information of sample {row['sampleid']} is missing, please provide it.<br>"
        return self

    def __call__(self):
        """执行检查流程"""
        try:
            (
                self
                .query_csv()
                .query_lims()
                .check_csv()
                .check_lims()
                .check_csv_lims()
            )
            if self.out_info:
                self.write_output(self.out_info)
                sys.exit(1)
            else:
                with open(self.out_file, "w", encoding="utf-8") as fw:
                    fw.write("")
                sys.exit(0)
        except Exception as e:
            self.write_output(f"Error: {str(e)}")
            sys.exit(1)


if __name__ == "__main__":
    ap = ArgumentParser()
    ap.add_argument("-f", "--fc_id", type=str, required=True)
    args = ap.parse_args()
    CheckTngs(fc_id=args.fc_id)()
