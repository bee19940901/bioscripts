using ArgParse

mutable struct GetNa
    infile::String
    outfile::String
    infas::Dict
    outfas::Dict
    function GetNa(infile::String, outfile::String, infas=Dict(), outfas=Dict())
        # 创建输出目录
        mkpath(dirname(outfile))
        return new(infile, outfile, infas, outfas)
    end
end

# 读取fasta文件保存到字典中
function readFasta(getNa::GetNa)
    head::String = ""
    open(getNa.infile, "r") do fr
        for line in readlines(fr)
            line = strip(line)
            if line !== ""
                if startswith(line, ">")
                    head = replace(line, r"^>"=>"")
                    getNa.infas[head] = ""
                else
                    getNa.infas[head] *= line
                end
            end
        end
    end
    return getNa
end

# 挑选序列中的NA序列
function filterFasta(getNa::GetNa)
    getNa.outfas = Dict(k => v for (k, v) in getNa.infas if occursin(r"\s+segment\s+6,?\s+", k))
    return getNa
end

# 将序列写入fasta文件中
function writeFasta(getNa::GetNa)
    open(getNa.outfile, "w") do fw
        for (k, v) in getNa.outfas
            name, _ = splitext(basename(getNa.infile))
            write(fw, ">$name\n$v\n")
        end
    end
    return getNa
end

# 主函数
function main()
    aps = ArgParseSettings()
    @add_arg_table! aps begin
        "--infile", "-i"
            arg_type = String
            required = true
            help = "Input fasta file."
        "--outfile", "-o"
            arg_type = String
            required = true
            help = "Output fasta file."
    end
    args = parse_args(aps)
    getNa = GetNa(abspath(args["infile"]), abspath(args["outfile"]))
    getNa |> readFasta |> filterFasta |> writeFasta
end

# 运行
main()