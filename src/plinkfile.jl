

export Marker,
    Sample,
    PlinkFile,
    samples,
    markers,
    nsamples,
    nmarkers,
    sample_index,
    marker_index,
    gt_hom1,
    gt_het,
    gt_hom2,
    gt_missing

using Mmap

"""
PLINK marker info
"""
struct Marker
    chrom::String
    id::String
    cm::Float64
    pos::Int64
    a1::String
    a2::String
end


"""
load marker info from PLINK .bim file
"""
function load_bim(path)
    
    markers = Marker[]
    bimpath = endswith(path, ".bim") ? path : path * ".bim"
    fh = open(bimpath, "r")
    for line in eachline(fh)
        cols = split(chomp(line))
        push!(markers, Marker(cols[1],
                              cols[2],
                              parse(Float64, cols[3]),
                              parse(Int64, cols[4]),
                              cols[5],
                              cols[6]))
    end
    close(fh)

    markers
end


"""
PLINK sample info
"""
struct Sample
    fid::String
    iid::String
    father::String
    mother::String
    sex::Int8
    phenotype::Int8
end


"""
load sample info from PLINK .fam file

columns are 
- fid
- iid
- iid of father
- iid of mother
- sex (1 = male, 2 = female, 0 = unknown)
- phenotype (1 = control, 2 = case, 0/-9 = missing)
"""
function load_fam(path)

    samples = Sample[]
    fampath = endswith(path, ".fam") ? path : path * ".fam"
    fh = open(fampath, "r")
    for line in eachline(fh)
        cols = split(chomp(line))
        push!(samples, Sample(cols[1],
                              cols[2],
                              cols[3],
                              cols[4],
                              parse(UInt8, cols[5]),
                              parse(Int8, cols[6])))
    end
    close(fh)

    samples
end


"""
PLINK files (.bim, .fam, .bed)
"""
struct PlinkFile
    path::String
    nmarkers::Int64
    nsamples::Int64
    markers::Vector{Marker}
    samples::Vector{Sample}
    A1::AbstractArray{Bool, 2}  # [sample, marker]; 
    A2::AbstractArray{Bool, 2}  # A2A1 = 00 (hom a1), 11 (hom a2), 10 (het), 01 (missing)
end

function Base.show(io::IO, p::PlinkFile)
    print(io, "<PLINK file ($(p.nsamples) samples x $(p.nmarkers) markers) at $(p.path)>")
end


samples(p::PlinkFile) = p.samples
markers(p::PlinkFile) = p.markers
nsamples(p::PlinkFile) = p.nsamples
nmarkers(p::PlinkFile) = p.nmarkers

"""
return index of sample in PlinkFile (linear search)
"""
function sample_index(p::PlinkFile, fid::String, iid::String)
    s = samples(p)
    for i = 1:length(s)
        if s[i].iid == iid && s[i].fid == fid
            return i
        end
    end
    error("No such sample: $fid $iid")
end


"""
return index of marker in PlinkFile (linear search)
"""
function marker_index(p::PlinkFile, id::String)
    m = markers(p)
    for i = 1:length(m)
        if m[i].id == id
            return i
        end
    end
    error("No such marker: $iid")
end    



function PlinkFile(path)

    samples = load_fam(path)
    markers = load_bim(path)
    
    nsamples = length(samples)
    nmarkers = length(markers)

    fh = open(path * ".bed", "r")
    bedheader = read(fh, 3)

    @assert bedheader[1] == 0b01101100 && bedheader[2] == 0b00011011 && bedheader[3] == 0b01 "Only snp major BED v1.0 supported"

    # the region of the bed file containing the genotypes 
    # might not be aligned correctly (to 8 bytes)
    #
    # we first need to map the file as UInt8, and 
    # then reinterpret the array. Need to make sure that
    # the UInt8 array  has size as multiple of 8 to 
    # allow reinterpretation as UInt64 

    n_uint8_per_marker = ceil(Int, nsamples/4)
    nsamples_ = 4 * n_uint8_per_marker

    n_uint64 = ceil(Int, n_uint8_per_marker * nmarkers / 8)
  
    chunks = Mmap.mmap(fh, Vector{UInt8}, 8 * n_uint64)
  
    data = BitArray{3}(undef, 0, 0, 0)
    data.len = 2*nsamples_*nmarkers
    data.dims = (2, nsamples_, nmarkers)
    data.chunks = reinterpret(UInt64, chunks)
  
    # use genotype A2A1 as of Plink convention
    A1 = view(data, 1, 1:nsamples, :)
    A2 = view(data, 2, 1:nsamples, :)

    PlinkFile(path,
              nmarkers,
              nsamples,
              markers,
              samples,
              A1,
              A2)
end

# genotype accessors (A2, A1)
@inline Base.getindex(p::PlinkFile, i::Int, j::Int) = (p.A2[i,j], p.A1[i,j])
@inline Base.getindex(p::PlinkFile, i, j) = (p.A2[i,j], p.A1[i,j])


# want to create 3bit presentation for fast computations of kinship
#
# hom1[marker, sample] to allow fast access to all markers of given sample
#


gt_hom1(A1, A2) = .~(A1 .| A2) 
gt_het(A1, A2) = (.~A1) .& A2
gt_hom2(A1, A2) = A1 .& A2
gt_missing(A1, A2) = A1 .& (.~A2)


# generate .ped and .map files
# this is only used for testing

# .ped format
# one line per sample
# first 6 columns from .fam file 
# then 2 columns for each genotype, alleles as strings

# .map file are the first 4 columns from .bim file

"""
generate .map and .ped files (for testing)
"""
function make_ped(p::PlinkFile, outfile)

    # write .map file
    fh = open(outfile * ".map", "w")
    for m in markers(p)
        write(fh, "$(m.chrom) $(m.id) $(m.cm) $(m.pos)\n")
    end 
    close(fh)

    # write .ped file
    fh = open(outfile * ".ped", "w")
    line = IOBuffer()
    for si in 1:nsamples(p)
        s = samples(p)[si]
        print(line, "$(s.fid)\t$(s.iid)\t$(s.father)\t$(s.mother)\t$(s.sex)\t$(s.phenotype)")
        for mi in 1:nmarkers(p)
            m = markers(p)[mi]
            gt = p[si, mi] # (A2, A1)
            if gt == (false, false)      # A2A1 = 00
                print(line, "\t$(m.a1)\t$(m.a1)")
            elseif gt == (true, false)   # A2A1 = 10
                print(line, "\t$(m.a1)\t$(m.a2)")
            elseif gt == (false, true)   # A2A1 = 01
                print(line, "\t0\t0")
            elseif gt == (true, true)    # A2A1 = 11
                print(line, "\t$(m.a2)\t$(m.a2)")
            end 
        end
        print(line, "\n")
        write(fh, String(take!(line)))
    end
    close(fh)

end