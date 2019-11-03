using PiecewiseIncreasingRanges, Test, Interpolations

@testset "PiecewiseIncreasingRange"  begin

    @testset "constructor" begin
        function test(rgs, divisor...)
            vcrgs = vcat(rgs...)
            !isempty(divisor) && (vcrgs /= divisor[1])
            rg = PiecewiseIncreasingRange(rgs, divisor...)
            @test length(rg.ranges) == 4
            @test vcrgs == rg

            yi = LinearInterpolation(1:length(vcrgs),convert(Vector{Float64}, vcrgs))
            @test resample(rg, 3//7) ≈ yi(1:7//3:length(yi))
            @test resample(rg, 7//3) ≈ yi(1:3//7:length(yi))
            @test resample(rg, 5//9) ≈ yi(1:9//5:length(yi))
            @test resample(rg, 61//3) ≈ yi(1:3//61:length(yi))

            for i = 1:length(rg)
                @test searchsortedfirst(rg, rg[i]) == i
                @test searchsortedfirst(rg, rg[i]-1//17) == i
                @test searchsortedfirst(rg, rg[i]+1//17) == i+1

                @test searchsortedlast(rg, rg[i]) == i
                @test searchsortedlast(rg, rg[i]-1//17) == i-1
                @test searchsortedlast(rg, rg[i]+1//17) == i

                for within_half_sample in (false, true)
                    @test findnearest(rg, rg[i], within_half_sample) == i
                    @test findnearest(rg, rg[i]+1//16, within_half_sample) == i
                    @test findnearest(rg, rg[i]-1//17, within_half_sample) == i
                    @test findnearest(rg, rg[i]+1//17, within_half_sample) == i
                end

                for j = i:length(rg)
                    @test vcrgs[i:j] == rg[i:j]
                    for k = 1:i-j+1
                        @test vcrgs[i:k:j] == rg[i:k:j]
                    end
                end
            end

            @test findnearest(rg, -1) == 1
            @test_throws NoNearestSampleError findnearest(rg, -1, true)
            @test findnearest(rg, 15.9) == searchsortedfirst(rg, 15)
            @test findnearest(rg, 16.1) == searchsortedfirst(rg, 17)
            @test_throws NoNearestSampleError findnearest(rg, 15.9, true)
            @test_throws NoNearestSampleError findnearest(rg, 16.1, true)
            @test findnearest(rg, 19) == length(rg)
            @test_throws NoNearestSampleError findnearest(rg, 19, true)
        end

        test(StepRange{Int,Int}[0:4:40, 41:1:80, 82:2:112, 114:2:120, 136:4:144], 8)
        test(StepRange{Rational{Int},Rational{Int}}[0:1//2:5, 5+1//8:1//8:10, 10+1//4:1//4:14, 14+1//4:1//4:15, 17:1//2:18])
        test(StepRangeLen{Float64}[0:1//2:5., 5+1//8:1//8:10., 10+1//4:1//4:14., 14+1//4:1//4:15., 17:1//2:18.])

        # Empty test
        rg = PiecewiseIncreasingRange(StepRange{Rational{Int},Rational{Int}}[])
        @test searchsortedfirst(rg, 1) == 1
        @test searchsortedlast(rg, 1) == 1
        @test_throws NoNearestSampleError findnearest(rg, 1)
        rg = PiecewiseIncreasingRange(UnitRange{Int}[0:-1])
        @test isempty(rg)

        # Test with non-monotonic range
        @test_throws ArgumentError PiecewiseIncreasingRange(StepRange{Rational{Int},Rational{Int}}[0:1//2:1, 3//4:1//2:4])
        @test_throws ArgumentError PiecewiseIncreasingRange(StepRange{Rational{Int},Rational{Int}}[0:-1//2:-1, 10:1//2:8])
        @test_throws ArgumentError PiecewiseIncreasingRange(StepRange{Rational{Int},Rational{Int}}[0:1//2:1, 10:-1//2:8])
        @test_throws ArgumentError PiecewiseIncreasingRange(StepRangeLen{Float64}[0:0.5:1, 0.75:0.5:4])
        @test_throws ArgumentError PiecewiseIncreasingRange(StepRangeLen{Float64}[0:-0.5:-1, 10:0.5:8])
        @test_throws ArgumentError PiecewiseIncreasingRange(StepRangeLen{Float64}[0:0.5:1, 10:-0.5:8])

        # Test with empty ranges interspersed with non-empty ranges
        rg = PiecewiseIncreasingRange(UnitRange{Int}[0:-1, 1:0, 0:-1, 1:3, 0:-1, 5:10])
        @test rg.ranges == UnitRange{Int}[1:3, 5:10]
    end

    @testset "Slicing using StepRange" begin
        vurgs = [1:3,6:7,10:12]
        vcrgs = vcat(vurgs...)
        p = PiecewiseIncreasingRange(vurgs)
        inds = 1:2:5
        ps = p[inds]
        @test ps == vcrgs[inds]
        @test isa(ps,PiecewiseIncreasingRange)
        @test ps.ranges[1] == 1:2:3
        @test ps.ranges[2] == 7:2:7
        @test isa(p[reverse(inds)],Vector)
        @test reverse(p[inds]) == p[reverse(inds)]

        inds_incr = IncreasingStepRange(inds)
        @test p[inds_incr] == ps
        ps = p[inds_incr]
        @test ps == vcrgs[inds_incr]
        @test isa(ps,PiecewiseIncreasingRange)
        @test ps.ranges[1] == 1:2:3
        @test ps.ranges[2] == 7:2:7
        @test isa(p[reverse(inds_incr)],Vector)
        @test reverse(p[inds_incr]) == p[reverse(inds_incr)]
    end

    @testset "checkindex" begin
        p = PiecewiseIncreasingRange([0:4:20, 61:15:80], 8)
        @test checkindex(Bool,p,1)
        @test checkindex(Bool,p,1.5)
        @test !checkindex(Bool,p,4)
        @test checkindex(Bool,p,:)
        @test checkindex(Bool,p,Base.Slice(1:3))

        v=Vector{Bool}(undef,1); v.=false;
        @test !checkindex(Bool,p,v)
        v=Array{Bool,2}(undef,1,1); v.=false;
        @test !checkindex(Bool,p,v)

        @test_throws ArgumentError checkindex(Bool,p,"a")
    end
end

@testset "PiecewiseIncreasingOrdinalRange" begin

    @testset "PiecewiseUnitRange" begin

        @testset "constructor" begin
            function test(rgs)
                vcrgs = vcat(rgs...)
                rg = PiecewiseUnitRange(rgs);
                @test length(rg.ranges)==3
                @test vcrgs == rg

                for i = 1:length(rg)
                    @test searchsortedfirst(rg, rg[i]) == i

                    @test searchsortedlast(rg, rg[i]) == i

                    for within_half_sample in (false, true)
                        @test findnearest(rg, rg[i], within_half_sample) == i
                    end

                    for j = i:length(rg)
                        @test vcrgs[i:j] == rg[i:j]
                        for k = 1:i-j+1
                            @test vcrgs[i:k:j] == rg[i:k:j]
                        end
                    end
                end

                @test findnearest(rg, -1) == 1
                @test_throws NoNearestSampleError findnearest(rg, -1, true)
                
                val = rg.ranges[1][end]
                @test findnearest(rg, val - 0.1) == searchsortedfirst(rg, val)
                @test findnearest(rg, val + 0.1) == searchsortedfirst(rg, val)
                
                val = (rg.ranges[1][end] + rg.ranges[2][1])/2
                @test_throws NoNearestSampleError findnearest(rg, val, true)
                @test findnearest(rg, last(rg)+1) == length(rg)
                @test_throws NoNearestSampleError findnearest(rg,last(rg)+1, true)
            end
            vurgs = [1:3,6:7,10:12]
            test(vurgs)
        end

        @testset "Slicing using StepRange" begin
            vurgs = [1:3,6:7,10:12]
            vcrgs = vcat(vurgs...)
            p = PiecewiseUnitRange(vurgs)
            inds = 1:2:5
            ps = p[inds]
            @test ps == vcrgs[inds]
            @test isa(ps,PiecewiseIncreasingStepRange)
            @test ps.ranges[1] == 1:2:3
            @test ps.ranges[2] == 7:2:7
            @test isa(p[reverse(inds)],Vector)
            @test reverse(p[inds]) == p[reverse(inds)]

            inds_incr = IncreasingStepRange(inds)
            @test p[inds_incr] == p[inds]
            ps = p[inds_incr]
            @test ps == vcrgs[inds_incr]
            @test isa(ps,PiecewiseIncreasingStepRange)
            @test ps.ranges[1] == 1:2:3
            @test ps.ranges[2] == 7:2:7
            @test isa(p[reverse(inds_incr)],Vector)
            @test reverse(p[inds_incr]) == p[reverse(inds_incr)]
        end

        @testset "checkindex" begin
            p = PiecewiseUnitRange([1:3,6:7,10:12])
            @test checkindex(Bool,p,1)
            @test checkindex(Bool,p,1.5)
            @test !checkindex(Bool,p,4)
            @test checkindex(Bool,p,:)
            @test checkindex(Bool,p,2:3)
            @test checkindex(Bool,p,2:1:3)
            @test !checkindex(Bool,p,2:6)
            @test !checkindex(Bool,p,2:2:6)
            @test checkindex(Bool,p,Base.Slice(1:3))

            v=Vector{Bool}(undef,last(p)); v.=false; @. v[p] = true;
            @test !checkindex(Bool,p,v)
            v=Array{Bool,2}(undef,last(p),1); v.=false;
            @test !checkindex(Bool,p,v)

            @test_throws ArgumentError checkindex(Bool,p,"a")
        end
    end

    @testset "PiecewiseIncreasingStepRange" begin
        @testset "constructor" begin
            function test(rgs)
                vcrgs = vcat(rgs...)
                rg = PiecewiseIncreasingStepRange(rgs);
                @test length(rg.ranges)==3
                @test vcrgs == rg

                for i = 1:length(rg)
                    @test searchsortedfirst(rg, rg[i]) == i

                    @test searchsortedlast(rg, rg[i]) == i

                    for within_half_sample in (false, true)
                        @test findnearest(rg, rg[i], within_half_sample) == i
                    end

                    for j = i:length(rg)
                        @test vcrgs[i:j] == rg[i:j]
                        for k = 1:i-j+1
                            @test vcrgs[i:k:j] == rg[i:k:j]
                        end
                    end
                end

                @test findnearest(rg, -1) == 1
                @test_throws NoNearestSampleError findnearest(rg, -1, true)
                
                val = rg.ranges[1][end]
                Δ = min(step(rg.ranges[1])/10,(rg.ranges[2][1]-rg.ranges[1][end])/10)
                @test findnearest(rg, val - Δ) == searchsortedfirst(rg, val)
                @test findnearest(rg, val + Δ) == searchsortedfirst(rg, val)
                
                val = (rg.ranges[1][end] + rg.ranges[2][1])/2
                @test_throws NoNearestSampleError findnearest(rg, val, true)

                Δ = 2step(rg.ranges[end])
                @test findnearest(rg, last(rg) + Δ) == length(rg)
                @test_throws NoNearestSampleError findnearest(rg,last(rg) + Δ, true)
            end

            vsrgs = [1:2:5,16:3:22,25:1:28]
            test(vsrgs)
        end

        @testset "Slicing using StepRange" begin
            vrgs = [1:2:5,10:3:16]
            vcrgs = vcat(vrgs...)
            p = PiecewiseIncreasingStepRange(vrgs);
            inds = 1:2:5
            ps = p[inds]
            @test ps == vcrgs[inds]
            @test isa(ps,PiecewiseIncreasingStepRange)
            @test isa(p[reverse(inds)],Vector)
            @test reverse(p[inds]) == p[reverse(inds)]

            inds_incr = IncreasingStepRange(inds)
            @test p[inds_incr] == ps
            ps = p[inds_incr]
            @test ps == vcrgs[inds_incr]
            @test isa(ps,PiecewiseIncreasingStepRange)
            @test ps.ranges[1] == 1:4:5
            @test ps.ranges[2] == 13:6:13
            @test isa(p[reverse(inds_incr)],Vector)
            @test reverse(p[inds_incr]) == p[reverse(inds_incr)]
        end

        @testset "checkindex" begin
            p = PiecewiseIncreasingStepRange([1:2:5,10:3:16])
            @test checkindex(Bool,p,1)
            @test checkindex(Bool,p,1.5)
            @test checkindex(Bool,p,2)
            @test !checkindex(Bool,p,6)
            @test checkindex(Bool,p,:)
            @test checkindex(Bool,p,1:2:5)
            @test checkindex(Bool,p,1:1:5)
            @test checkindex(Bool,p,1:5)
            @test checkindex(Bool,p,Base.Slice(1:5))

            v=Vector{Bool}(undef,last(p)); v.=false; @. v[p] = true;
            @test !checkindex(Bool,p,v)
            v=Array{Bool,2}(undef,last(p),1); v.=false;
            @test !checkindex(Bool,p,v)

            @test_throws ArgumentError checkindex(Bool,p,"a")
        end
    end
end