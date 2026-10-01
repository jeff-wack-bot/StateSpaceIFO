using StateSpaceIFO
using LinearAlgebra
using Test

band(lo, hi, n) = exp10.(range(log10(lo), log10(hi); length = n))

@testset "StateSpaceIFO.jl" begin

    @testset "Kimura worked examples" begin
        # Ex 4.4: Θ = diag((s-2)/(s+2), (s+1)/(s-1)) is J_{1,1}-lossless, P = diag(1/4, 1/2)
        Θ44 = SS([-2.0 0; 0 1], [-4.0 0; 0 2], [1.0 0; 0 1], [1.0 0; 0 1])
        r = jjlossless_residual(Θ44, 1, 1)
        @test r.dD < 1.0e-12
        @test r.dLyap < 1.0e-12
        @test r.dFreq < 1.0e-12
        @test lyap(transpose(Θ44.A), transpose(Θ44.C) * Jsig(1, 1) * Θ44.C) ≈ [0.25 0; 0 0.5] atol = 1.0e-10

        # Ex 4.6: degree-2 J-lossless example, P = diag(1, 1/2)
        Θ46 = SS([0.0 -1; -2 -1], [-1.0 -1; -5 / 2 3 / 2], [1.0 5 / 4; -1 3 / 4], [1.0 0; 0 1])
        r = jjlossless_residual(Θ46, 1, 1)
        @test r.dD < 1.0e-12
        @test r.dLyap < 1.0e-12
        @test r.dFreq < 1.0e-12
        @test r.minP ≈ 0.5 atol = 1.0e-10

        # Ex 6.3/6.4: (J,J')-lossless factorization of an unstable G
        G64 = SS([-2.0 0; 0 1], [-2.0 2; 0 2], [2.0 0; 0 2], [1.0 -1; 0 1])
        f = jj_factor(G64, 1, 1)
        @test f.ok
        @test f.X ≈ [1 0; 0 0] atol = 1.0e-9
        @test f.X̄ ≈ [0 0; 0 0.5] atol = 1.0e-9
        Πref(s) = [1 -1; 0 (s + 3) / (s + 1)]
        for Ω in [0.1, 0.7, 1.3, 5.0, 20.0]
            s = im * Ω
            @test norm(evalfr(f.Θ, s) * evalfr(f.Π, s) - evalfr(G64, s)) < 1.0e-9
            # Θ, Π are unique only up to a constant J'-unitary factor: compare Π'J'Π
            @test norm(
                evalfr(f.Π, s)' * Jsig(1, 1) * evalfr(f.Π, s) -
                    Πref(s)' * Jsig(1, 1) * Πref(s)
            ) < 1.0e-8
        end
        r = jjlossless_residual(f.Θ, 1, 1)
        @test r.dD < 1.0e-9 && r.dLyap < 1.0e-9 && r.dFreq < 1.0e-9

        # Ex 6.5: a tall (J_{21}, J_{11})-lossless factorization
        s2 = sqrt(2)
        G65 = SS(
            [3.0 0; 0 -3], [-s2 4 - s2; s2 / 2 s2 / 2], [0.0 -4; 1 -2; 2 -1],
            [s2 / 2 s2 / 2; s2 / 2 s2 / 2; 0 1]
        )
        f5 = jj_factor(G65, 2, 1)
        @test f5.ok
        @test f5.X ≈ [1 / 14 0; 0 2] atol = 1.0e-9
        @test f5.X̄ ≈ [2 0; 0 0] atol = 1.0e-9
        for Ω in [0.05, 0.4, 1.0, 3.0, 11.0]
            @test norm(evalfr(f5.Θ, im * Ω) * evalfr(f5.Π, im * Ω) - evalfr(G65, im * Ω)) < 1.0e-8
        end
    end

    @testset "optical elements" begin
        # a lossless passive cavity is inner
        for (δ, Δ) in [(1.0, 0.0), (0.5156013216, -0.5156013216 * 1.7671310136), (1.3, 2.0)]
            G = filter_cavity(δ, Δ)
            e = maximum(norm(evalfr(G, im * Ω)' * evalfr(G, im * Ω) - I) for Ω in band(0.01, 100, 400))
            @test e < 1.0e-10
        end

        # quadrature matrix = exp(-i α_com) R(α_rot)
        δ, Δ = 0.7, -0.9
        ξ = -Δ / δ
        for Ω in [0.05, 0.3, 1.0, 2.5, 8.0]
            αr, αc = rotation_angle(ξ, δ, Ω), common_phase(ξ, δ, Ω)
            αp, αm = sideband_phases(ξ, δ, Ω)
            @test αr ≈ (αp + αm) / 2
            @test αc ≈ (αp - αm) / 2
            R = [cos(αr) -sin(αr); sin(αr) cos(αr)]
            @test norm(evalfr(filter_cavity(δ, Δ), im * Ω) - exp(-im * αc) * R) < 1.0e-10
        end

        # the rotation phasor is the Blaschke factor in s = Ω²; the difference is not
        Δ, δ = 0.7, 0.4
        ξ = -Δ / δ
        p = (Δ - im * δ)^2
        for Ω in [0.1, 0.35, 0.9, 1.7, 3.3]
            αp, αm = sideband_phases(ξ, δ, Ω)
            @test abs(cis(αp + αm) - (Ω^2 - conj(p)) / (Ω^2 - p)) < 1.0e-12
        end
        @test maximum(
            abs(cis(-(sideband_phases(ξ, δ, Ω)...)) - (Ω^2 - conj(p)) / (Ω^2 - p))
                for Ω in [0.1, 0.35, 0.9, 1.7, 3.3]
        ) > 0.1
        # a tuned cavity does not rotate
        @test maximum(abs(rotation_angle(0.0, 0.4, Ω)) for Ω in [0.1, 1.0, 5.0]) < 1.0e-15

        # a squeezer is symplectic, and J-unitary with J = diag(1,-1) in the (a, a†) basis
        Sq = [exp(0.8) 0.0; 0.0 exp(-0.8)]
        W = [1 im; 1 -im] / sqrt(2)
        Bg = W * Sq * inv(W)
        Jb = Diagonal([1.0, -1.0])
        Σ = [0.0 1; -1 0]
        @test norm(Sq'Sq - I) > 1
        @test norm(transpose(Sq) * Σ * Sq - Σ) < 1.0e-12
        @test norm(Bg' * Jb * Bg - Jb) < 1.0e-12
    end

    @testset "KLMTV filter cavities (Eq. 89)" begin
        table = [
            (0.25, 1.4041853665, 0.3741199545, -0.168117389, 1.0812267357),
            (1.0, 1.7671310136, 0.5156013216, -0.2772297408, 1.3017526994),
            (10.0, 2.1681653525, 0.8512981061, -0.3687198181, 2.06433239),
        ]
        for (I0, ξI, δI, ξII, δII) in table
            c = klmtv_cavities(I0)
            @test c[1].ξ ≈ ξI atol = 1.0e-9
            @test c[1].δ ≈ δI atol = 1.0e-9
            @test c[2].ξ ≈ ξII atol = 1.0e-9
            @test c[2].δ ≈ δII atol = 1.0e-9
        end

        # the same cavities from the rational K = 2/(x² + x), x = Ω²
        c = cavities_from_rational([2.0], [0.0, 1.0, 1.0])
        @test c[1].ξ ≈ 1.7671310136 atol = 1.0e-9
        @test c[2].ξ ≈ -0.2772297408 atol = 1.0e-9

        # two cavities + homodyne at π/2 reach the variational optimum exactly
        c = klmtv_cavities(1.0)
        p = [c[1].ξ, c[1].δ, c[2].ξ, c[2].δ, π / 2]
        @test maximum(abs(excess_ratio(p, Ω) - 1) for Ω in band(0.02, 50, 400)) < 1.0e-9

        # the cascade is inner, and Potapov splitting recovers cavity I
        F = filter_cavity(c[1].δ, c[1].Δ) * filter_cavity(c[2].δ, c[2].Δ)
        @test nx(F) == 4
        @test jjlossless_residual(F, 2, 2; Ωs = band(0.01, 100, 300)).dFreq < 1.0e-9
        p1 = poles(filter_cavity(c[1].δ, c[1].Δ))
        Θ1, Θ2 = potapov_split(F, 2, 2, λ -> minimum(abs.(λ .- p1)) < 1.0e-6)
        for Ω in band(0.01, 100, 200)
            @test norm(evalfr(Θ1, im * Ω) * evalfr(Θ2, im * Ω) - evalfr(F, im * Ω)) < 1.0e-8
            @test norm(evalfr(Θ1, im * Ω)' * evalfr(Θ1, im * Ω) - I) < 1.0e-9
        end
        @test sort(poles(Θ1); by = imag) ≈ sort(p1; by = imag) atol = 1.0e-8
    end

    @testset "KLMTV plant" begin
        P = klmtv_plant()
        for Ω in [0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0]
            Ga, _, Gh = transfer_matrices(P, Ω)
            @test real(backaction_gain(Ga)) ≈ klmtv_K(Ω) rtol = 1.0e-12
            @test arm_phase(Ga) ≈ atan(Ω) atol = 1.0e-12
            @test abs(Ga[1, 2]) < 1.0e-14
            @test abs(Ga[2, 2] - Ga[1, 1]) < 1.0e-14
            Σ = Ga * Ga'
            @test sh_min(Gh, Σ) ≈ 1 / abs2(Gh[2]) rtol = 1.0e-9
            v = optimal_readout(Gh, Σ)
            @test real(v[1] / v[2]) ≈ klmtv_K(Ω) rtol = 1.0e-9
        end

        # closed form and canonical realization of K, with a suspension
        Ps = klmtv_plant(μ = 0.02, k = 0.02)
        Kss = ponderomotive_gain_ss(μ = 0.02, k = 0.02)
        for Ω in band(0.02, 50, 100)
            Kp = ponderomotive_gain(im * Ω; μ = 0.02, k = 0.02)
            @test backaction_gain(evalfr(Ps.Gn, im * Ω)) ≈ Kp rtol = 1.0e-9
            @test evalfr(Kss, im * Ω)[1] ≈ Kp rtol = 1.0e-9
        end
        @test ponderomotive_gain(im * 0.7) ≈ klmtv_K(0.7)
    end

    @testset "reality defect" begin
        # tuned arm, loss at the input coupler: rho = 0 and no real-readout penalty
        for γl in [1.0e-3, 1.0e-2, 1.0e-1], Ω in [0.2, 1.0, 5.0]
            g = readout_gap(klmtv_plant(γl = γl), Ω)
            @test g.ρ < 1.0e-12
            @test g.χdef < 1.0e-12
            @test g.Sreal ≈ g.Sany rtol = 1.0e-9
        end
        # detuned arm: rho = O(1), and a real readout pays for it
        g = readout_gap(klmtv_plant(γl = 1.0e-2, Δc = 1.0), 1.0)
        @test g.ρ > 0.1
        @test g.χdef > 1.0e-3
        @test 10log10(g.Sreal / g.Sany) ≈ 0.589 atol = 1.0e-3
    end

    @testset "marginal pole (Lemma 5.5)" begin
        for (μ, k, exists) in [(0.0, 0.0, false), (0.1, 0.0, false), (0.0, 0.01, false), (0.01, 0.01, true), (0.2, 1.0, true)]
            P = klmtv_plant(μ = μ, k = k)
            Θc, _ = jlossless_conjugator(P.A, [P.Ba P.Bh], Jsig(1, 2))
            @test (Θc !== nothing) == exists
        end
        for (μ, k, exists) in [(0.0, 0.0, false), (0.1, 0.0, false), (0.0, 0.1, false), (0.02, 0.02, true)]
            H = dual_chain_plant(ponderomotive_gain_ss(μ = μ, k = k), 1.05)
            @test dual_jj_factor(H, 1, 1; strict = false).ok == exists
        end
    end

    @testset "dual factorization of the KLMTV plant" begin
        Kss = ponderomotive_gain_ss(μ = 0.02, k = 0.02)
        J = Jsig(1, 2)
        # det(H J H~) = 1 - lev², and (6.50) reports feasibility iff lev > 1
        for lev in [0.5, 0.9, 0.999, 1.001, 1.05, 1.5, 5.0]
            H = dual_chain_plant(Kss, lev)
            for Ω in band(0.02, 50, 50)
                T = evalfr(H, im * Ω) * J * transpose(evalfr(H, -im * Ω))
                @test abs(det(T) - (1 - lev^2)) < 1.0e-9
            end
            @test (signature_sqrt(H.D * J * transpose(H.D), 1, 1) !== nothing) == (lev > 1)
        end

        d = dual_factor_design(1.05; μ = 0.02, k = 0.02)
        @test d.ok
        @test d.eH < 1.0e-8
        @test d.eJ < 1.0e-8
        r = dual_lossless_residual(d.Ψ, 1, 1; Ωs = band(0.02, 50, 200))
        @test r.dD < 1.0e-9
        @test r.dFreq < 1.0e-7
        @test isstable(minreal(d.Ω).A)
        @test isstable(minreal(d.iΩ).A)
        @test sqrt(d.supR) <= 1.05

        @test !dual_factor_design(0.99; μ = 0.02, k = 0.02).ok

        # q -> K as lev -> 1
        Ωs = band(0.05, 20, 400)
        errs = map([1.0e-2, 1.0e-3, 1.0e-4]) do δ
            dd = dual_factor_design(1 + δ; μ = 0.02, k = 0.02, Ωs = Ωs)
            maximum(abs(dd.q[i] - ponderomotive_gain(im * Ωs[i]; μ = 0.02, k = 0.02)) for i in eachindex(Ωs))
        end
        @test issorted(errs; rev = true)
        @test errs[end] < 1.0e-3
    end

    @testset "quadratic invariance" begin
        Z2 = zeros(2, 2)
        arm(Ω) = (1 - im * Ω) / (1 + im * Ω) * [1.0 0; -2 / (Ω^2 * (1 + Ω^2)) 1]
        srm = 0.9 * exp(0.3im) * Matrix(1.0I, 2, 2)
        G(Ω) = [Z2 srm; arm(Ω) Z2]
        Ωs = [0.05, 0.2, 0.7, 1.0, 2.5, 8.0, 30.0]
        defect(basis, Gf) = maximum(first(qi_defect(basis, Gf(Ω))) for Ω in Ωs)

        @test defect(qi_basis(:full), G) < 1.0e-10
        @test defect(qi_basis(:tied), Ω -> zeros(4, 4)) < 1.0e-10
        @test defect(qi_basis(:untied), Ω -> [arm(Ω) Z2; Z2 srm]) < 1.0e-10
        @test defect(qi_basis(:single), Ω -> srm * arm(Ω)) < 1.0e-10
        @test defect(qi_basis(:antidiag), G) < 1.0e-10
        @test defect(qi_basis(:tied), G) > 1 - 1.0e-10
        @test defect(qi_basis(:untied), G) > 1 - 1.0e-10
        @test_throws ArgumentError qi_basis(:nope)
    end
end
