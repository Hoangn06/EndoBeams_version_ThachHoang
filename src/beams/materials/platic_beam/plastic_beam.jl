function compute_plastic_contributions( 
    K̄ⁱⁿᵗ, ū, Θ̅₁, Θ̅₂, B¹, B¹¹, B²¹, B¹², B²², B¹⁴, B²⁴, exact, Tₛ⁻¹Θ̅₁, Tₛ⁻¹Θ̅₂,
    P¹¹, P²¹, P¹², P²², P¹⁴, P²⁴, η, lₙ, D₃, Gᵀ¹, Gᵀ², Gᵀ⁴, Rₑ, r¹,
    B̄⁺¹¹, B̄⁺²¹, B̄⁺¹², B̄⁺²², B̄⁺¹⁴, B̄⁺²⁴,
    radius,
    E, G, ν, y0, K, l₀, nᴳ, zᴳ, ωᴳ, nˢ, qˢ, rˢ, tˢ, wˢ,
    plastic_state::PlasticStates, beam_ind::Int=0
)
    strain_energy = 0.0
    local_converged = true

    T̄ⁱⁿᵗ = zeros(Vec7{Float64})
    K̄ⁱⁿᵗ = zeros(Mat77{Float64})

    id_G = 1
    for i in 1:nᴳ
        x1 = (l₀/2)*(zᴳ[i]+1)
        for j in 1:nˢ
            for k in 1:nˢ
                x2 = radius * qˢ[j] * rˢ[k]
                x3 = radius * tˢ[j] * rˢ[k]

                sigma, C_const, converged = plasticity_return_mapping(E, G, ν, y0, K,
                    plastic_state,
                    id_G, x1, x2, x3, l₀,
                    ū, Θ̅₁, Θ̅₂, beam_ind)

                if !converged
                    local_converged = false
                end

                Be   = matrix_Be(x1, x2, x3, l₀)
                pref = (l₀/2)*((pi*radius^2)/nˢ)
                wtot = ωᴳ[i]*wˢ[k]

                T̄ⁱⁿᵗ += (pref*wtot) * (Be' * sigma)
                K̄ⁱⁿᵗ += (pref*wtot) * (Be' * C_const * Be)

                id_G += 1
            end
        end
    end

    # Extract components from T̄ⁱⁿᵗ (Vec7): [ū, Θ̅₁[1:3], Θ̅₂[1:3]]
    T̄ⁱⁿᵗū  = T̄ⁱⁿᵗ[1]
    T̄ⁱⁿᵗΘ̅₁ = Vec3(T̄ⁱⁿᵗ[2], T̄ⁱⁿᵗ[3], T̄ⁱⁿᵗ[4])
    T̄ⁱⁿᵗΘ̅₂ = Vec3(T̄ⁱⁿᵗ[5], T̄ⁱⁿᵗ[6], T̄ⁱⁿᵗ[7])

    # Transform forces: Tⁱⁿᵗ = Bᵀ T̄ⁱⁿᵗ
    Tⁱⁿᵗ¹ = B¹'*T̄ⁱⁿᵗū + B¹¹'*T̄ⁱⁿᵗΘ̅₁ + B²¹'*T̄ⁱⁿᵗΘ̅₂
    Tⁱⁿᵗ² =             B¹²'*T̄ⁱⁿᵗΘ̅₁ + B²²'*T̄ⁱⁿᵗΘ̅₂
    Tⁱⁿᵗ³ = -Tⁱⁿᵗ¹
    Tⁱⁿᵗ⁴ =             B¹⁴'*T̄ⁱⁿᵗΘ̅₁ + B²⁴'*T̄ⁱⁿᵗΘ̅₂

    Tⁱⁿᵗ = [Tⁱⁿᵗ¹;
            Tⁱⁿᵗ²;
            Tⁱⁿᵗ³;
            Tⁱⁿᵗ⁴]

    # Compute [N̄ M̄⁺₁ M̄⁺₂] = B̄ᵀ T̄ⁱⁿᵗ (for geometric stiffness)
    N̄   = T̄ⁱⁿᵗū
    M̄⁺₁ = exact ? Tₛ⁻¹Θ̅₁' * T̄ⁱⁿᵗΘ̅₁ : T̄ⁱⁿᵗΘ̅₁
    M̄⁺₂ = exact ? Tₛ⁻¹Θ̅₂' * T̄ⁱⁿᵗΘ̅₂ : T̄ⁱⁿᵗΘ̅₂

    # Qₛ = Pᵀ [M̄⁺₁ M̄⁺₂]
    Qₛ¹ = P¹¹' * M̄⁺₁ + P²¹' * M̄⁺₂
    Qₛ² = P¹²' * M̄⁺₁ + P²²' * M̄⁺₂
    Qₛ⁴ = P¹⁴' * M̄⁺₁ + P²⁴' * M̄⁺₂

    # Q = S(Qₛ)
    Q¹ = skew(Qₛ¹)
    Q² = skew(Qₛ²)
    Q⁴ = skew(Qₛ⁴)

    a = @SVector [0, η*(M̄⁺₁[1] + M̄⁺₂[1])/lₙ + (M̄⁺₁[2] + M̄⁺₂[2])/lₙ, (M̄⁺₁[3] + M̄⁺₂[3])/lₙ]

    # DN̄
    DN̄¹¹ = D₃*N̄

    # QGᵀ
    QGᵀ¹¹ = Q¹*Gᵀ¹
    QGᵀ¹² = Q¹*Gᵀ²
    QGᵀ¹⁴ = Q¹*Gᵀ⁴
    QGᵀ²² = Q²*Gᵀ²
    QGᵀ²⁴ = Q²*Gᵀ⁴
    QGᵀ⁴⁴ = Q⁴*Gᵀ⁴

    # EGa
    EGa¹ = Rₑ*Gᵀ¹'*a

    # EGar
    EGar¹¹ = EGa¹*r¹

    # Kₘ = DN̄ - EQGᵀEᵀ + EGar (geometric stiffness)
    Kₘ¹¹ = DN̄¹¹ - Rₑ*QGᵀ¹¹*Rₑ' + EGar¹¹
    Kₘ¹² =      - Rₑ*QGᵀ¹²*Rₑ'
    Kₘ¹⁴ =      - Rₑ*QGᵀ¹⁴*Rₑ'
    Kₘ²² =      - Rₑ*QGᵀ²²*Rₑ'
    Kₘ²⁴ =      - Rₑ*QGᵀ²⁴*Rₑ'
    Kₘ⁴⁴ =      - Rₑ*QGᵀ⁴⁴*Rₑ'

    # K̃ (geometric stiffness with exact rotation)
    if exact
        η₁, μ₁ = compute_η_μ(Θ̅₁)
        η₂, μ₂ = compute_η_μ(Θ̅₂)

        M̄₁ = T̄ⁱⁿᵗΘ̅₁
        M̄₂ = T̄ⁱⁿᵗΘ̅₂

        K̄ₕ₁ = compute_K̄ₕ(Θ̅₁, M̄₁, Tₛ⁻¹Θ̅₁, η₁, μ₁)
        K̄ₕ₂ = compute_K̄ₕ(Θ̅₂, M̄₂, Tₛ⁻¹Θ̅₂, η₂, μ₂)
    end

    K̃¹¹ = exact ?  B̄⁺¹¹' * K̄ₕ₁ * B̄⁺¹¹  +  B̄⁺²¹' * K̄ₕ₂ * B̄⁺²¹  +  Kₘ¹¹  :  Kₘ¹¹
    K̃¹² = exact ?  B̄⁺¹¹' * K̄ₕ₁ * B̄⁺¹²  +  B̄⁺²¹' * K̄ₕ₂ * B̄⁺²²  +  Kₘ¹²  :  Kₘ¹²
    K̃¹⁴ = exact ?  B̄⁺¹¹' * K̄ₕ₁ * B̄⁺¹⁴  +  B̄⁺²¹' * K̄ₕ₂ * B̄⁺²⁴  +  Kₘ¹⁴  :  Kₘ¹⁴
    K̃²² = exact ?  B̄⁺¹²' * K̄ₕ₁ * B̄⁺¹²  +  B̄⁺²²' * K̄ₕ₂ * B̄⁺²²  +  Kₘ²²  :  Kₘ²²
    K̃²⁴ = exact ?  B̄⁺¹²' * K̄ₕ₁ * B̄⁺¹⁴  +  B̄⁺²²' * K̄ₕ₂ * B̄⁺²⁴  +  Kₘ²⁴  :  Kₘ²⁴
    K̃⁴⁴ = exact ?  B̄⁺¹⁴' * K̄ₕ₁ * B̄⁺¹⁴  +  B̄⁺²⁴' * K̄ₕ₂ * B̄⁺²⁴  +  Kₘ⁴⁴  :  Kₘ⁴⁴

    zero_row = SMatrix{1,3,Float64}(0.0, 0.0, 0.0)
    B = vcat(
        hcat(B¹, zero_row, -B¹, zero_row),
        hcat(B¹¹, B¹², -B¹¹, B¹⁴),
        hcat(B²¹, B²², -B²¹, B²⁴)
    )

    BᵀKⁱⁿᵗB = B' * K̄ⁱⁿᵗ * B

    K̃ = hcat(vcat(K̃¹¹, K̃¹²', -K̃¹¹', K̃¹⁴'),
             vcat(K̃¹², K̃²², -K̃¹², K̃²⁴'),
             vcat(-K̃¹¹, -K̃¹²', K̃¹¹, -K̃¹⁴'),
             vcat(K̃¹⁴, K̃²⁴, -K̃¹⁴, K̃⁴⁴))

    Kⁱⁿᵗ = BᵀKⁱⁿᵗB + K̃

    return Tⁱⁿᵗ, Kⁱⁿᵗ, strain_energy, local_converged
end
