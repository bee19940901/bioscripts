using DataFrames
using ArgParse
using Random

#=
从所有的甲流病毒参考序列中获取部分NA序列用于建树
=#

# 定义 TreeRef 结构体
mutable struct TreeRef
    infile::String
    outfile::String
    infas::Dict
    outdf::DataFrame
    function TreeRef(infile::String, outfile::String; infas=Dict(), outdf=DataFrame())
        return new(infile, outfile, infas, outdf)
    end
end

# 创建输出目录
function mkdirs(treeRef::TreeRef)
    outdir = dirname(treeRef.outfile)
    if !isdir(outdir)
        mkpath(outdir)
        println("目录已创建： $outdir")
    end
    return treeRef
end

# 读取fasta文件到字典中
function readFasta(treeRef::TreeRef)
    head = ""
    for line in eachline(treeRef.infile)
        line = strip(line)
        if line != ""
            if startswith(line, ">")
                head = replace(line, r"^>" => "")
                treeRef.infas[head] = ""
            else
                treeRef.infas[head] *= line
            end
        end
    end
    return treeRef
end

# 将序列字典转换为数据框格式，筛选出需要的序列，并将数据框输出到fasta文件中
function filterFas(treeRef::TreeRef)
    # 获取基因组中的第六个片段
    fas = treeRef.infas  # 使用 treeRef.infas 字段
    treeRef.infas = Dict(k => split(v, r"N+")[6] for (k, v) in fas if occursin(r"\(H1N1\)", k))

    # 将字典转换为数据框
    records = []
    for (k, v) in treeRef.infas
        m = match(r"(\d{4})\(H1N1\)", k)
        if m !== nothing
            year = m.captures[1]
            acce = split(k, r"\s+")[1]
            push!(records, (year=year, acce=acce, name="$year/$acce", seq=v))
        end
    end
    treeRef.outdf = DataFrame(records)

    # 筛选需要的序列
    grouped = groupby(treeRef.outdf, :year)
    sampled_groups = vcat([g[shuffle(1:nrow(g))[1:min(5, nrow(g))], :] for g in grouped]...)

    # 将序列输出到fasta文件
    open(treeRef.outfile, "w") do fw
        for row in eachrow(sampled_groups)
            id = row[:name]
            seq = row[:seq]
            write(fw, ">$id\n$seq\n")
        end
    end

    return treeRef
end

# 主函数
function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--infile", "-i"
            help = "Input FASTA file"
            arg_type = String
            required = true
        "--outfile", "-o"
            help = "Output FASTA file"
            arg_type = String
            required = true
    end
    args = parse_args(s)
    tree_ref = TreeRef(args["infile"], args["outfile"])
    tree_ref |> mkdirs |> readFasta |> filterFas
end

main()
