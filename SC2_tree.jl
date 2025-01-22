using ArgParse

mutable struct SC2Tree
    infas::String
    outfas::String
    function SC2Tree(infas::String, outfas::String)
        mkpath(dirname(outfas))
        return new(infas, outfas)
    end
end

# 将fasta文件转换为字典
function readFasta(infile::String) :: Dict{String, String}
    header::String = ""
    fasDict::Dict{String, String} = Dict()
    for line::String in eachline(infile)
        line = strip(line)
        if line != ""
            if startswith(line, ">")
                header = replace(line, r"^>"=>"")
                fasDict[header] = ""
            else
                fasDict[header] *= line
            end
        end
    end
    return fasDict
end

# 获取文件的stem
function getSamplename(infile::String) :: String
    return replace(basename(infile), r"(?:\.\S+)+$"=>"")
end

# 将fasta字典写入文件中
function writefasta(fasDict:: Dict{String, String}, outFas:: String)
    sampleName::String = getSamplename(outFas)
    initNum::Int64 = 1
    open(outFas, "w") do fw
        for (_, v) in fasDict
            write(fw, ">$(sampleName)-$(initNum)\n$(v)\n")
            initNum += 1
        end
    end
end

# CALL
function run(sc2Tree::SC2Tree)
    fasDict::Dict{String, String} = readFasta(sc2Tree.infas)
    writefasta(fasDict, sc2Tree.outfas)
end

function main()
    ap = ArgParseSettings()
    @add_arg_table! ap begin
        "--infas", "-i"
            arg_type = String
            required = true
            help = "Input Fasta File"
        "--outfas", "-o"
            arg_type = String
            required = true
            help = "Output Fasta File"
    end
    args = parse_args(ap)
    sc2Tree = SC2Tree(abspath(args["infas"]), abspath(args["outfas"]))
    sc2Tree |> run
end

main()