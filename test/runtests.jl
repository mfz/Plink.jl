using Plink
using Test

@testset "Plink.jl" begin
    p = PlinkFile("../data/plink")
    @test nsamples(p) == 1056
    @test nmarkers(p) == 4

    Plink.make_ped(p, "/tmp/plink")
    @test success(`cmp -s ../data/plink.ped /tmp/plink.ped`)
    @test success(`cmp -s ../data/plink.map /tmp/plink.map`)
    rm("/tmp/plink.ped")
    rm("/tmp/plink.map")
end

