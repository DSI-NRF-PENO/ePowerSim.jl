# (C) 2025 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za

# Sauer: section 6.10:  pg 135 - 136, 6.242 - 6.247

####################################################

#-----------------------------------------------------
#-----------------------------------------------------
# generic algebraic eq i.e network and stator equations
# This is powerflow and netwkork equations
#-----------------------------------------------------
#-----------------------------------------------------

# Millano , page 329

function common_equations(
    dx, x, (θh, vh, id, iq, τm), t; kwd)

    δ, ω = x

    # Ωb = 2 * pi * 60
    
    # ωs = 1.0
    
    Ωb, D, H, ωs = kwd

    M = 2 * H
    
    vd = vh * sin(δ -  θh)
    
    vq = vh * cos(δ -  θh)
    
    τe = ψd * iq - ψq * id

    # 0 = τm0 - τm
    
    # 0 = vf0 - vf
    
    dx[1] = Ωb * (ω -  ωs)

    dx[2] = (1/M) * (τm -  τe - D * (ω -  ωs))

    return nothing

end


function stator_equations_by_flux_dynamics(
    dx, x, (ω, vd, vq, id, iq), t; kwd)

    ψd, ψq = x
    
    (ra, Ωb) = kwd
    
    dx[1] = Ωb * (ra * id + ω * ψq + vd)

    dx[2] = Ωb * (ra * iq - ω * ψd + vq)

    return nothing

end


function stator_equations_by_zero_flux_dynamics(
    dx, x, (ω, vd, vq, id, iq), t; kwd)

    ψd, ψq = x
    
    (ra, ) = kwd
    
    dx[1] = ra * id + ω * ψq + vd

    dx[2] = ra * iq - ω * ψd + vq

    return nothing

end


function stator_equations_by_small_ω_approximation(
    dx, x, (vd, vq, id, iq), t; kwd)

    ψd, ψq = x
    
    (ra, ) = kwd
    
    dx[1] = ra * id + ψq + vd

    dx[2] = ra * iq - ψd + vq

    return nothing

end


# "T_d_dash": 8.96,
# "X_d_dash": 0.0608,
# "D": 0.01254,
# "X_d": 0.146,
# "Ωb": 314.159,
# "T_q_dash": 0.31,
# "X_q_dash": 0.0969,
# "ra": 0.000125,
# "X_q": 0.1000,
# "X_q_2dash": 0.06,
# "ωs": 376.99111843077515,
# "xℓ": 0.0146,
# "T_d_2dash": 0.01,
# "T_q_2dash": 0.01,
# "X_d_2dash": 0.06,



function magnetic_equations_by_sauer_pai_model(
    dx, x, (vf, vd, vq, id, iq), t; kwd)

    ed_dash, eq_dash, ψd_2dash, ψq_2dash = x
    
    (;X_d,
     X_q,
     X_d_dash,
     X_q_dash,
     X_d_2dash,
     X_q_2dash,
     xℓ,
     T_d_dash,
     T_q_dash,
     T_d_2dash,
     T_q_2dash) =
         NamedTupleTools.select(
             kwd,
             (:X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash,
              :X_d_2dash,
              :X_q_2dash,
              :xℓ,
              :T_d_dash,
              :T_q_dash,
              :T_d_2dash,
              :T_q_2dash))

    γd1 = (X_d_2dash - xℓ)/(X_d_dash - xℓ )

    γq1 = (X_q_2dash - xℓ)/(X_q_dash - xℓ )

    γd2 = ( 1 - γd1 )/(X_d_dash - xℓ )

    γq2 = ( 1 - γq1 )/(X_q_dash - xℓ )

    # ed_dash, eq_dash, ψd_2dash, ψq_2dash
    # X_d, X_q,
    
    dx[1] = (-eq_dash - (X_d - X_d_dash)*(
        id - γd2 * ψd_2dash- (1 - γd1) * id +
            γd2 * eq_dash) + vf )/T_d_dash

    dx[2] = (-ed_dash + (X_q - X_q_dash)*(
        iq - γq2 * ψq_2dash- (1 - γq1) * iq +
            γd2 * ed_dash ))/T_q_dash

    dx[3] = (-ψd_2dash + eq_dash - (
        X_d_dash - xℓ) * id)/T_d_2dash

    dx[4] = (-ψq_2dash - ed_dash - (
        X_q_dash - xℓ) * iq)/T_q_2dash

    
    return nothing

end


function algebraic_magnetic_equations_by_sauer_pai_model(
    dy, y, (x, id, iq), t; kwd)

    ψd, ψq = y
    
    ed_dash, eq_dash, ψd_2dash, ψq_2dash = x    
    
    (;X_d,
     X_q,
     X_d_dash,
     X_q_dash,
     X_d_2dash,
     X_q_2dash,
     xℓ,
     T_d_dash,
     T_q_dash,
     T_d_2dash,
     T_q_2dash) =
         NamedTupleTools.select(
             kwd,
             (:X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash,
              :X_d_2dash,
              :X_q_2dash,
              :xℓ,
              :T_d_dash,
              :T_q_dash,
              :T_d_2dash,
              :T_q_2dash))

    γd1 = (X_d_2dash - xℓ)/(X_d_dash - xℓ )

    γq1 = (X_q_2dash - xℓ)/(X_q_dash - xℓ )

    # ed_dash, eq_dash, ψd_2dash, ψq_2dash
    
    # X_d, X_q,

    dy[1] = ψd + X_d_2dash * id - γd1 * eq_dash - (
        1 - γd1) * ψd_2dash

    dy[2] = ψq + X_q_2dash * iq + γq1 * ed_dash - (
        1 - γq1) * ψq_2dash
    
    return nothing

end



function magnetic_equations_by_marconato_model(
    dx, x, (vf, vd, vq, id, iq), t; kwd)

    ed_dash, eq_dash, ed_2dash, eq_2dash = x
    
    (;X_d,
     X_q,
     X_d_dash,
     X_q_dash,
     X_d_2dash,
     X_q_2dash,
     xℓ,
     T_d_dash,
     T_q_dash,
     T_d_2dash,
     T_q_2dash) =
         NamedTupleTools.select(
             kwd,
             (:X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash,
              :X_d_2dash,
              :X_q_2dash,
              :xℓ,
              :T_d_dash,
              :T_q_dash,
              :T_d_2dash,
              :T_q_2dash))

    γd = (T_d_2dash / T_d_dash) * (
        X_d_2dash / X_d_dash) * (X_d - X_d_dash)

    γq = (T_q_2dash / T_q_dash) * (
        X_q_2dash / X_q_dash) * (X_q - X_q_dash)
    
    dx[1] = (-eq_dash - (X_d - X_d_dash - γd)* id + (
         1 - TAA/T_d_dash)vf )/T_d_dash

    dx[2] = (-ed_dash + (X_q - X_q_dash - γq)* iq)/T_q_dash

    dx[3] = (-eq_2dash + eq_dash - (X_d_dash - X_d_2dash +
        γd)* id + (TAA/T_d_dash)vf )/T_d_2dash

    dx[4] = (-ed_2dash + ed_dash + (X_q_dash - X_q_2dash +
        γq)* iq )/T_q_2dash
    
    return nothing

end



function algebraic_magnetic_equations_by_marconato_model(
    dy, y, (x, id, iq), t; kwd)

    ψd, ψq = y
    
    ed_dash, eq_dash, ed_2dash, eq_2dash = x    
    
    (;X_d,
     X_q,
     X_d_dash,
     X_q_dash,
     X_d_2dash,
     X_q_2dash,
     xℓ,
     T_d_dash,
     T_q_dash,
     T_d_2dash,
     T_q_2dash) =
         NamedTupleTools.select(
             kwd,
             (:X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash,
              :X_d_2dash,
              :X_q_2dash,
              :xℓ,
              :T_d_dash,
              :T_q_dash,
              :T_d_2dash,
              :T_q_2dash))

    dy[1] = ψd + X_d_2dash * id - eq_2dash

    dy[2] = ψq + X_q_2dash * iq + ed_2dash
    
    return nothing

end



function magnetic_equations_by_anderson_fouad_model(
    dx, x, (vf, vd, vq, id, iq), t; kwd)

    ed_dash, eq_dash, ed_2dash, eq_2dash = x
    
    (;X_d,
     X_q,
     X_d_dash,
     X_q_dash,
     X_d_2dash,
     X_q_2dash,
     xℓ,
     T_d_dash,
     T_q_dash,
     T_d_2dash,
     T_q_2dash) =
         NamedTupleTools.select(
             kwd,
             (:X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash,
              :X_d_2dash,
              :X_q_2dash,
              :xℓ,
              :T_d_dash,
              :T_q_dash,
              :T_d_2dash,
              :T_q_2dash))
    
    dx[1] = (-eq_dash - (X_d - X_d_dash)* id + vf )/T_d_dash

    dx[2] = (-ed_dash + (X_q - X_q_dash)* iq)/T_q_dash

    dx[3] = (-eq_2dash + eq_dash -
        (X_d_dash - X_d_2dash) * id )/T_d_2dash

    dx[4] = (-ed_2dash + ed_dash + (X_q_dash - X_q_2dash)* iq )/T_q_2dash
    
    return nothing

end


function algebraic_magnetic_equations_by_anderson_fouad_model(
    dy, y, (x, id, iq), t; kwd)

    ψd, ψq = y
    
    ed_dash, eq_dash, ed_2dash, eq_2dash = x    
    
    (;X_d,
     X_q,
     X_d_dash,
     X_q_dash,
     X_d_2dash,
     X_q_2dash,
     xℓ,
     T_d_dash,
     T_q_dash,
     T_d_2dash,
     T_q_2dash) =
         NamedTupleTools.select(
             kwd,
             (:X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash,
              :X_d_2dash,
              :X_q_2dash,
              :xℓ,
              :T_d_dash,
              :T_q_dash,
              :T_d_2dash,
              :T_q_2dash))

    γd1 = (X_d_2dash - xℓ)/(X_d_dash - xℓ )

    γq1 = (X_q_2dash - xℓ)/(X_q_dash - xℓ )

    dy[1] = ψd + X_d_2dash * id - eq_2dash

    dy[2] = ψq + X_q_2dash * iq - ed_2dash
    
    return nothing

end


#----------------------------------------
# Synchronous generator
#----------------------------------------


function gen_fun__SM_2axis_cb_v6__init(
    gen_vh,
    gen_θh,
    P_g,
    Q_g;
    kwd_para =
        gen_para)
    
    (ra, X_d, X_q, X_d_dash, X_q_dash) =
        NamedTupleTools.select(
            kwd_para,
            (:ra,
             :X_d,
             :X_q,
             :X_d_dash,
             :X_q_dash))

    
    Ig = (P_g - im * Q_g) / (
        gen_vh * exp(-im * gen_θh))

    E_gen = gen_vh * exp(im * gen_θh) +
        (ra + im * X_q) * Ig

    δ = angle( E_gen )

    E = abs(E_gen)

    i_dq =  abs(Ig) * exp(im * (angle(Ig) - δ + π/2) )

    v_dq =  gen_vh * exp(im * (gen_θh - δ + π/2) )

    i_d = real(i_dq )

    i_q = imag(i_dq )

    v_d = real(v_dq )

    v_q = imag(v_dq )
    
    ed_dash = v_d  + ra * i_d - X_q_dash * i_q
    
    eq_dash = v_q  + ra * i_q + X_d_dash * i_d

    vf_tilade = eq_dash + (X_d - X_d_dash) * i_d

    τm_tilade = ed_dash * i_d + eq_dash * i_q +
        (X_q_dash - X_d_dash) * i_d * i_q
    
    return  (;δ,  E,
             i_d, i_q,
             v_d, v_q,
             ed_dash, eq_dash,
             τm_tilade,
             vf_tilade) 

end

#---------------------------------------------

function ode_gen_fun__SM_2axis_cb_v6__func!(
    dx, x,
    gen_id_iq_τmtilade_vftilade_para,
    t;
    kwd_para =
        gen_kwd_para )
    
    (i_d,
     i_q,
     τm_tilade,
     vf_tilade) =
        gen_id_iq_τmtilade_vftilade_para
    
    (gen_para,
     ωs) =
         kwd_para
    
    (;H,
     D,
     
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gen_para,
             (:H,
              :D,
              
              :X_d,
              :X_q,
              
              :X_d_dash,
              :X_q_dash,
              
              :T_d_dash,
              :T_q_dash))

    # Inputs:  i_d, i_q, Tm, vf_tilade    
    # Output:  gen_ω === x[4]
    
    δ, ω, ed_dash, eq_dash = x

    dx[1] = ω - ωs
    
    dx[2] = ( ωs/(2*H)) * (τm_tilade - ed_dash *
        i_d - eq_dash *
        i_q - (X_q_dash - X_d_dash) * i_d * i_q  -
        D * (ω - ωs) ) 
    
    dx[3] = (1.0 / T_q_dash) *
        ( (X_q - X_q_dash) * i_q  - ed_dash)

    dx[4] = (1/T_d_dash) *
        (vf_tilade - (X_d - X_d_dash) * i_d  -
        eq_dash)
    
    return nothing    

end



function dae_gen_fun__SM_2axis_cb_v6__func!(
    res,
    dx, x,
    gen_id_iq_τmtilade_vftilade_para,
    t;
    kwd_para =
        gen_kwd_para )
    
    (i_d,
     i_q,
     τm_tilade,
     vf_tilade) =
        gen_id_iq_τmtilade_vftilade_para
    
    (gen_para,
     ωs) =
         kwd_para
    
    (;H,
     D,
     
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gen_para,
             (:H,
              :D,
              
              :X_d,
              :X_q,
              
              :X_d_dash,
              :X_q_dash,
              
              :T_d_dash,
              :T_q_dash))

    # Inputs:  i_d, i_q, Tm, vf_tilade    
    # Output:  gen_ω === x[4]
    
    δ, ω, ed_dash, eq_dash = x
    
    # if abs( dx[2] ) < Δω_toleraance
        
    #     Δω = ωs - ω
        
    #     res[1] = ω - ωs - dx[1] +  Δω

    # else
        
    #     res[1] = ω - ωs - dx[1]
        
    # end

    res[1] = ω - ωs - dx[1]
    
    res[2] = ( ωs/(2*H)) * (τm_tilade - ed_dash *
        i_d - eq_dash *
        i_q - (X_q_dash - X_d_dash) * i_d * i_q -
        D * (ω - ωs) ) - dx[2]
    
    res[3] = (1/T_q_dash) *
        ( (X_q - X_q_dash) * i_q - ed_dash) - dx[3]

    res[4] = (1/T_d_dash) *
        (vf_tilade - (X_d - X_d_dash) * i_d  -
        eq_dash) - dx[4]    
    
    return nothing    

end



function gen_fun__SM_2axis_cb_v6__func!(
    res,
    dx, x,
    gen_id_iq_τmtilade_vftilade_para,
    t;
    kwd_para =
        gen_kwd_para )
    
    (i_d,
     i_q,
     τm_tilade,
     vf_tilade) =
        gen_id_iq_τmtilade_vftilade_para
    
    (gen_para,
     ωs) =
         kwd_para
    
    (;H,
     D,
     
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gen_para,
             (:H,
              :D,
              
              :X_d,
              :X_q,
              
              :X_d_dash,
              :X_q_dash,
              
              :T_d_dash,
              :T_q_dash))

    # Inputs:  i_d, i_q, Tm, vf_tilade    
    # Output:  gen_ω === x[4]
    
    δ, ω, ed_dash, eq_dash = x

    res[1] = ω - ωs - dx[1]
    
    res[2] = (ωs/(2*H)) * (τm_tilade - ed_dash *
        i_d - eq_dash *
        i_q - (X_q_dash - X_d_dash) * i_d * i_q -
        D * (ω - ωs) ) - dx[2]
    
    res[3] = (1/T_q_dash) *
        ( (X_q - X_q_dash) * i_q - ed_dash) - dx[3]

    res[4] = (1/T_d_dash) *
        (vf_tilade - (X_d - X_d_dash) * i_d  -
        eq_dash) - dx[4]    
    
    return nothing    

end



function gen_fun__SM_2axis_cb_v6__output(
    x,
    ωref0_vref0_porder0_id_iq_vh;
    kwd_para =
        gen_kwd_para )
    
    δ, ω, ed_dash, eq_dash = x
    
    return  ω 

end


function gen_fun__SM_2axis_cb_v6__callback_paras(
    plant_state_idx_in_states,
    plant_states_syms ;
    comp_paras =
        gen_kwd_para )


    return [make_a_component_state_no_callback_paras(
        :nothing)]

end

#----------------------------------------
# Synchronous condenser
#----------------------------------------


function gen_fun__SC_2axis_cb_v6__init(
    gen_vh,
    gen_θh,
    P_g,
    Q_g;
    kwd_para =
        gen_para)
    
    (ra, X_d, X_q, X_d_dash, X_q_dash) =
        NamedTupleTools.select(
            kwd_para,
            (:ra,
             :X_d,
             :X_q,
             :X_d_dash,
             :X_q_dash))

    
    Ig = (P_g - im * Q_g) / (
        gen_vh * exp(-im * gen_θh))

    E_gen = gen_vh * exp(im * gen_θh) +
        (ra + im * X_q) * Ig

    δ = angle( E_gen )

    E = abs(E_gen)

    i_dq =  abs(Ig) * exp(im * (angle(Ig) - δ + π/2) )

    v_dq =  gen_vh * exp(im * (gen_θh - δ + π/2) )

    i_d = real(i_dq )

    i_q = imag(i_dq )

    v_d = real(v_dq )

    v_q = imag(v_dq )

    #   ed_dash = (X_q - X_q_dash) * i_q
    
    ed_dash = v_d  + ra * i_d - X_q_dash * i_q
    
    eq_dash = v_q  + ra * i_q + X_d_dash * i_d

    vf_tilade = eq_dash + (X_d - X_d_dash) * i_d

    τm_tilade = ed_dash * i_d + eq_dash * i_q +
        (X_q_dash - X_d_dash) * i_d * i_q
    
    return  (;δ,  E,
             i_d, i_q,
             v_d, v_q,
             ed_dash, eq_dash,
             τm_tilade,
             vf_tilade) 

end

#---------------------------------------------


function ode_gen_fun__SC_2axis_cb_v6__func!(
    dx, x,
    gen_id_iq_τmtilade_vftilade_para,
    t;
    kwd_para =
        gen_kwd_para )
    
    i_d, i_q, τm_tilade, vf_tilade =
        gen_id_iq_τmtilade_vftilade_para
    
    (gen_para,
     ωs) =
         kwd_para
    
    (;H,
     D,
     
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gen_para,
             (:H,
              :D,
              :X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash,
              :T_d_dash,
              :T_q_dash))

    # Inputs:  i_d, i_q, Tm, vf_tilade    
    # Output:  gen_ω === x[4]
    
    δ, ω, ed_dash, eq_dash = x

    dx[1] = ω - ωs
    
    dx[2] = (ωs/(2*H)) * (τm_tilade - ed_dash *
        i_d - eq_dash * i_q -
        (X_q_dash - X_d_dash) * i_d * i_q - D * (ω - ωs))
        
    dx[3] = (1/T_q_dash) *
        ( (X_q - X_q_dash) * i_q  - ed_dash)

    dx[4] = (1/T_d_dash) *
        (vf_tilade - (X_d - X_d_dash) * i_d  - eq_dash)
    
    return nothing    

end


function dae_gen_fun__SC_2axis_cb_v6__func!(
    res,
    dx, x,
    gen_id_iq_τmtilade_vftilade_para,
    t;
    kwd_para =
        gen_kwd_para )
    
    i_d, i_q, τm_tilade, vf_tilade =
        gen_id_iq_τmtilade_vftilade_para
    
    (gen_para,
     ωs) =
         kwd_para
    
    (;H,
     D,
     
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gen_para,
             (:H,
              :D,
              :X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash,
              :T_d_dash,
              :T_q_dash))

    # Inputs:  i_d, i_q, Tm, vf_tilade    
    # Output:  gen_ω === x[4]
    
    δ, ω, ed_dash, eq_dash = x
    
    # if abs( dx[2] ) < Δω_toleraance
        
    #     Δω = ωs - ω
        
    #     res[1] = ω - ωs - dx[1] +  Δω

    # else
        
    #     res[1] = ω - ωs - dx[1]
        
    # end

    res[1] = ω - ωs - dx[1]
    
    res[2] = ( ωs/(2*H)) * (τm_tilade - ed_dash *
        i_d - eq_dash *
        i_q - (X_q_dash - X_d_dash) * i_d * i_q -
        D * (ω - ωs))  - dx[2]
        
    res[3] = (1/T_q_dash) *
        ( (X_q - X_q_dash) * i_q  - ed_dash) - dx[3]

    res[4] = (1/T_d_dash) *
        (vf_tilade - (X_d - X_d_dash) * i_d  -
        eq_dash) - dx[4]
    
    
    return nothing    

end


function gen_fun__SC_2axis_cb_v6__func!(
    res,
    dx, x,
    gen_id_iq_τmtilade_vftilade_para,
    t;
    kwd_para =
        gen_kwd_para )
    
    i_d, i_q, τm_tilade, vf_tilade =
        gen_id_iq_τmtilade_vftilade_para
    
    (gen_para,
     ωs) =
         kwd_para
    
    (;H,
     D,
     
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gen_para,
             (:H,
              :D,
              :X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash,
              :T_d_dash,
              :T_q_dash))

    # Inputs:  i_d, i_q, Tm, vf_tilade    
    # Output:  gen_ω === x[4]
    
    δ, ω, ed_dash, eq_dash = x

    res[1] = ω - ωs - dx[1]
    
    res[2] = ( ωs/(2*H)) * (τm_tilade - ed_dash *
        i_d - eq_dash *
        i_q - (X_q_dash - X_d_dash) * i_d * i_q -
        D * (ω - ωs) )  - dx[2]

    
    # res[2] = ( ωs/(2*H)) * ( - D * (ω - ωs) )  - dx[2]
    
    res[3] = (1/T_q_dash) *
        ( (X_q - X_q_dash) * i_q  - ed_dash) - dx[3]

    res[4] = (1/T_d_dash) *
        (vf_tilade - (X_d - X_d_dash) * i_d  -
        eq_dash) - dx[4]
    
    
    return nothing    

end



function gen_fun__SC_2axis_cb_v6__output(
    x,
    ωref0_vref0_porder0_id_iq_vh;
    kwd_para =
        gen_kwd_para )
    
    δ, ω, ed_dash, eq_dash = x
    
    return  ω 

end


function gen_fun__SC_2axis_cb_v6__callback_paras(
    plant_state_idx_in_states,
    plant_states_syms ;
    comp_paras =
        gen_kwd_para )


    return [make_a_component_state_no_callback_paras(
        :nothing)]

end

#----------------------------------------
# Gov
#----------------------------------------

function gov_fun__gov_ieee_tgov1_cb__init(
    ωs, τm_tilade;
    kwd_para =
        gov_kwd_para )
    
    (T1,    
     T2,    
     T3,    
     Dt,    
     p_max, 
     p_min, 
     R) =
         NamedTupleTools.select(
             kwd_para,
             (:T1,    
              :T2,    
              :T3,    
              :Dt,    
              :p_max, 
              :p_min, 
              :R))


    ω_ref0  = ωs
    

    b = [0.0, τm_tilade]

    A = [(1 - T2/T3)  -1.0;
         (T2/T3)       1.0
         ]
    
    x = A \ b

    x1     = x[1]
    x2     = x[2]

    p_order = x1

    return (; state_var = Float64[x1, x2],
            ref = (; ω_ref0, p_order ) )

end

#---------------------------------------------

function ode_gov_fun__gov_ieee_tgov1_cb__func!(
    dx, x,
    gov_gen_ω_ωref0_porder0_para,
    t;
    kwd_para =
        gov_kwd_para_wt_gov_cb_sw )


    (gov_kwd_para,
     gov_cb_sw) = kwd_para

    sw = gov_cb_sw[1][1]

    # Inputs :  gen_ω, ω_ref0, p_order0    
    # Output : xg2 === Tm === x[2]    
         
    (gen_ω,
     ω_ref0,
     p_order0)   =
        gov_gen_ω_ωref0_porder0_para
    
    (T1,    
     T2,    
     T3,    
     Dt,    
     p_max, 
     p_min, 
     R) =
         NamedTupleTools.select(
             gov_kwd_para,
             (:T1,    
              :T2,    
              :T3,    
              :Dt,    
              :p_max, 
              :p_min, 
              :R))

    p_ref    = p_order0 * R

    xg1, xg2 = x

    # dx[1] = (1/T1) * (( p_order0 - (1/R) *
    #     ( gen_ω/ω_ref0 - 1)) - xg1 )
    
    if sw == 0

        dx[1] = (1/T1) * (( p_order0 - (1/R) *
            ( gen_ω/ω_ref0 - 1)) - xg1 )
        
    elseif sw == 1
        
        dx[1] = 0.0
        
    elseif sw == 2
        
        dx[1] = 0.0
    end    
    
    dx[2] = (1/T3) * (( 1 - T2/T3 ) *
        xg1 - xg2 )
   
    return nothing    

end



function dae_gov_fun__gov_ieee_tgov1_cb__func!(
    res,
    dx, x,
    gov_gen_ω_ωref0_porder0_para, 
    t;
    kwd_para =
        gov_kwd_para_wt_gov_cb_sw )


    (gov_kwd_para,
     gov_cb_sw) = kwd_para

    sw = gov_cb_sw[1][1]

    # Inputs :  gen_ω, ω_ref0, p_order0    
    # Output : xg2 === Tm === x[2]    
         
    (gen_ω,
     ω_ref0,
     p_order0)   =
        gov_gen_ω_ωref0_porder0_para

    
    (T1,    
     T2,    
     T3,    
     Dt,    
     p_max, 
     p_min, 
     R) =
         NamedTupleTools.select(
             gov_kwd_para,
             (:T1,    
              :T2,    
              :T3,    
              :Dt,    
              :p_max, 
              :p_min, 
              :R))

    p_ref    = p_order0 * R

    xg1, xg2 = x


    # res[1] = (1/T1) * (( p_order0 - (1/R) *
    #     ( gen_ω/ω_ref0 - 1)) - xg1 ) - dx[1]
    
    if sw == 0

        res[1] = (1/T1) * (( p_order0 - (1/R) *
            ( gen_ω/ω_ref0 - 1)) - xg1 ) - dx[1]
        
    elseif sw == 1

        dx[1] = 0.0
        res[1] = 0.0 - dx[1]
        
    elseif sw == 2

        dx[1] = 0.0
        res[1] = 0.0 - dx[1]
    end    
    
    res[2] = (1/T3) * (( 1 - T2/T3 ) *
        xg1 - xg2 ) - dx[1]
   
    return nothing    

end


function gov_fun__gov_ieee_tgov1_cb__func!(
    res,
    dx, x,
    gov_gen_ω_ωref0_porder0_para, 
    t;
    kwd_para =
        gov_kwd_para_wt_gov_cb_sw )


    (gov_kwd_para,
     gov_cb_sw) = kwd_para

    sw = gov_cb_sw[1][1]

    # Inputs :  gen_ω, ω_ref0, p_order0    
    # Output : xg2 === Tm === x[2]    
         
    (gen_ω,
     ω_ref0,
     p_order0)   =
        gov_gen_ω_ωref0_porder0_para

    
    (T1,    
     T2,    
     T3,    
     Dt,    
     p_max, 
     p_min, 
     R) =
         NamedTupleTools.select(
             gov_kwd_para,
             (:T1,    
              :T2,    
              :T3,    
              :Dt,    
              :p_max, 
              :p_min, 
              :R))

    p_ref    = p_order0 * R

    xg1, xg2 = x


    # res[1] = (1/T1) * (( p_order0 - (1/R) *
    #     ( gen_ω/ω_ref0 - 1)) - xg1 ) - dx[1]
    
    if sw == 0

        res[1] = (1/T1) * (( p_order0 - (1/R) *
            ( gen_ω/ω_ref0 - 1)) - xg1 ) - dx[1]
        
    elseif sw == 1

        dx[1] = 0.0
        res[1] = 0.0 - dx[1]
        
        
    elseif sw == 2

        dx[1] = 0.0
        res[1] = 0.0 - dx[1]
        
    end    
    
    res[2] = (1/T3) * (( 1 - T2/T3 ) *
        xg1 - xg2 ) - dx[1]
   
    return nothing    

end



function gov_fun__gov_ieee_tgov1_cb__output(
    x,
    gov_gen_ω_ωref0_porder0_para;
    kwd_para =
        gov_kwd_para )

    # Inputs :  gen_ω, ω_ref0, p_order0    
    # Output : xg2 === Tm === x[2]    
         
    (gen_ω,
     ω_ref0,
     p_order0)   =
        gov_gen_ω_ωref0_porder0_para

    
    (T1,    
     T2,    
     T3,    
     Dt,    
     p_max, 
     p_min, 
     R) =
         NamedTupleTools.select(
             kwd_para,
             (:T1,    
              :T2,    
              :T3,    
              :Dt,    
              :p_max, 
              :p_min, 
              :R))

    xg1, xg2 = x
    
    return xg2 + (T2/T3) * xg1 - Dt * (
        gen_ω/ω_ref0 - 1)
   

end


function gov_fun__gov_ieee_tgov1_cb__callback_paras(
    plant_state_idx_in_states,
    plant_states_syms;
    comp_paras =
        gov_kwd_para )

    component_parameter_switch = [0]
    
    (p_max, 
     p_min) =
         NamedTupleTools.select(
             comp_paras,
             (:p_max, 
              :p_min))
    
    lower_lim_anti_windup_cb_para =
        make_a_component_state_callback_paras(    
            :lower_lim_anti_windup,

            plant_state_idx_in_states,
            plant_states_syms,

            p_min,
            p_min,

            :xg1,
            :xg1,

            cb_fun_condition_lower_lim,
            cb_fun_affect_anti_windup_lower_lim!;

            component_parameter_switch =
                component_parameter_switch )
    
    upper_lim_anti_windup_cb_para =
        make_a_component_state_callback_paras(    
            :upper_lim_anti_windup,

            plant_state_idx_in_states,
            plant_states_syms,

            p_max,
            p_max,

            :xg1,
            :xg1,

            cb_fun_condition_upper_lim,
            cb_fun_affect_anti_windup_upper_lim!;

            component_parameter_switch =
                component_parameter_switch)

    
    return [ lower_lim_anti_windup_cb_para,
             upper_lim_anti_windup_cb_para]
   

end

#---------------------------------------------------

function gov_fun__gov_t1_cb__init(
    ωs, τm_tilade;
    kwd_para =
        gov_kwd_para )

    
    (T3,    
    T4,    
    T5,    
    Tc,    
    Ts,    
    p_max, 
    p_min, 
     R) =
         NamedTupleTools.select(
             kwd_para,
             (:T3,    
              :T4,    
              :T5,    
              :Tc,    
              :Ts,    
              :p_max, 
              :p_min, 
              :R))
     

    ω_ref0 = ωs
    

    b = [0.0, 0.0, τm_tilade]

    A = [(1 - T3/Tc)          (-1.0)        0.0;
         (1 - T4/T5)*(T3/Tc)  (1 - T4/T5)   -1.0;
         (T4/T5)*(T3/Tc)      (T4/T5)       1.0
         ]
    
    x = A \ b

    xg1     = x[1]
    xg2     = x[2]
    xg3     = x[3]

    p_order =  xg1 - 1/R * (ω_ref0/ωs - 1.0 )
    
    return (; state_var = Float64[xg1, xg2, xg3],
            ref = (; ω_ref0, p_order ))


end


function ode_gov_fun__gov_t1_cb__func!(
    dx, x,
    gov_gen_ω_ωref0_porder0_para,
    t;
    kwd_para =
        gov_kwd_para_wt_gov_cb_sw )


    (gov_kwd_para,
     gov_cb_sw) = kwd_para

    # Inputs :  gen_ω, ω_ref0, p_order0    
    # Output : xg2 === Tm === x[2]    
         
    (gen_ω,
     ω_ref0,
     p_order0)   =
        gov_gen_ω_ωref0_porder0_para

    
    (T3,    
    T4,    
    T5,    
    Tc,    
    Ts,    
    p_max, 
    p_min, 
     R) =
         NamedTupleTools.select(
             gov_kwd_para,
             (:T3,    
              :T4,    
              :T5,    
              :Tc,    
              :Ts,    
              :p_max, 
              :p_min, 
              :R))
     
 
    xg1, xg2, xg3 = x

    # phat_in = p_order0 - (1/R) * (u_gen_ω/ω_ref0-1)
    
    # res[1] = (1/Ts) * (phat_in - xg1) - dx[1]

    dx[1] = (1/Ts) * (p_order0 - (1/R) *
        (gen_ω/ω_ref0 - 1) - xg1)

    dx[2] = (1/Tc) * ((1 - T3/Tc) * xg1 - xg2)

    dx[3] = (1/T5) * ((1 - T4/T5) * (xg2 +
        (T3/Tc) * xg1) - xg3 )
    
    return nothing    

end


function dae_gov_fun__gov_t1_cb__func!(
    res,
    dx, x,
    gov_gen_ω_ωref0_porder0_para,
    t;
    kwd_para =
        gov_kwd_para_wt_gov_cb_sw )


    (gov_kwd_para,
     gov_cb_sw) = kwd_para

    # Inputs :  gen_ω, ω_ref0, p_order0    
    # Output : xg2 === Tm === x[2]    
         
    (gen_ω,
     ω_ref0,
     p_order0)   =
        gov_gen_ω_ωref0_porder0_para

    
    (T3,    
    T4,    
    T5,    
    Tc,    
    Ts,    
    p_max, 
    p_min, 
     R) =
         NamedTupleTools.select(
             gov_kwd_para,
             (:T3,    
              :T4,    
              :T5,    
              :Tc,    
              :Ts,    
              :p_max, 
              :p_min, 
              :R))
     
 
    xg1, xg2, xg3 = x

    # phat_in = p_order0 - (1/R) * (u_gen_ω/ω_ref0-1)
    
    # res[1] = (1/Ts) * (phat_in - xg1) - dx[1]

    res[1] = (1/Ts) * (p_order0 - (1/R) *
        (gen_ω/ω_ref0 - 1) - xg1) - dx[1] 

    res[2] = (1/Tc) * ((1 - T3/Tc) * xg1 - xg2) - dx[2]

    res[3] = (1/T5) * ((1 - T4/T5) * (xg2 +
        (T3/Tc) * xg1) - xg3 ) - dx[3]
    
    return nothing    

end



function gov_fun__gov_t1_cb__func!(
    res,
    dx, x,
    gov_gen_ω_ωref0_porder0_para,
    t;
    kwd_para =
        gov_kwd_para_wt_gov_cb_sw )


    (gov_kwd_para,
     gov_cb_sw) = kwd_para

    # Inputs :  gen_ω, ω_ref0, p_order0    
    # Output : xg2 === Tm === x[2]    
         
    (gen_ω,
     ω_ref0,
     p_order0)   =
        gov_gen_ω_ωref0_porder0_para

    
    (T3,    
    T4,    
    T5,    
    Tc,    
    Ts,    
    p_max, 
    p_min, 
     R) =
         NamedTupleTools.select(
             gov_kwd_para,
             (:T3,    
              :T4,    
              :T5,    
              :Tc,    
              :Ts,    
              :p_max, 
              :p_min, 
              :R))
     
 
    xg1, xg2, xg3 = x

    # phat_in = p_order0 - (1/R) * (u_gen_ω/ω_ref0-1)
    
    # res[1] = (1/Ts) * (phat_in - xg1) - dx[1]

    res[1] = (1/Ts) * (p_order0 - (1/R) *
        (gen_ω/ω_ref0 - 1) - xg1) - dx[1] 

    res[2] = (1/Tc) * ((1 - T3/Tc) * xg1 - xg2) - dx[2]

    res[3] = (1/T5) * ((1 - T4/T5) * (xg2 +
        (T3/Tc) * xg1) - xg3 ) - dx[3]
    
    return nothing    

end



function gov_fun__gov_t1_cb__output(
    x,
    gov_gen_ω_ωref0_porder0_para;
    kwd_para =
        gov_kwd_para )

    # Inputs :  gen_ω, ω_ref0, p_order0    
    # Output : xg3 + (T4/T5) * (xg2 + (T3/Tc) * xg1)  
         
    (gen_ω,
     ω_ref0,
     p_order0)   =
        gov_gen_ω_ωref0_porder0_para

    
    (T3,    
    T4,    
    T5,    
    Tc,    
    Ts,    
    p_max, 
    p_min, 
     R) =
         NamedTupleTools.select(
             kwd_para,
             (:T3,    
              :T4,    
              :T5,    
              :Tc,    
              :Ts,    
              :p_max, 
              :p_min, 
              :R))
     
 
    xg1, xg2, xg3  = x

    return xg3 + (T4/T5) * (xg2 + (T3/Tc) * xg1) 

end



function gov_fun__gov_t1_cb__callback_paras(
    plant_state_idx_in_states,
    plant_states_syms;
    comp_paras =
        gov_kwd_para )

    
    (p_max, 
     p_min) =
         NamedTupleTools.select(
             comp_paras,
             (:p_max, 
              :p_min))
    
    lower_lim_windup_cb_para =
        make_a_component_state_callback_paras(    
            :lower_lim_windup,

            plant_state_idx_in_states,
            plant_states_syms,

            p_min,
            p_min,

            :xg1,
            :xg1,

            cb_fun_condition_lower_lim,
            cb_fun_affect_windup_lower_lim! )
    
    upper_lim_windup_cb_para =
        make_a_component_state_callback_paras(    
            :upper_lim_windup,

            plant_state_idx_in_states,
            plant_states_syms,

            p_max,
            p_max,

            :xg1,
            :xg1,

            cb_fun_condition_upper_lim,
            cb_fun_affect_windup_upper_lim! )

    
    return [ lower_lim_windup_cb_para,
             upper_lim_windup_cb_para]
   

end

#---------------------------------------------------

function gov_fun__gov_t1_cb_sauer__init(
    ωs, τm_tilade;
    kwd_para =
        gov_kwd_para )

    
    (;Tc,   
     Ts,   
     p_max,
     p_min,
     R) =
         NamedTupleTools.select(
             kwd_para,
             (:Tc,   
              :Ts,   
              :p_max,
              :p_min,
              :R))     

    ω_ref0 =  ωs
    
    xg2 = τm_tilade 

    xg1 = xg2 

    p_order = xg1 

    return (state_var = Float64[xg1, xg2],
            ref = (;ω_ref0, p_order))
end



function gov_fun__gov_t1_cb_sauer__func!(
    res,
    dx, x,
    gov_gen_ω_ωref0_porder0_para, 
    t;
    kwd_para =
        gov_kwd_para_wt_gov_cb_sw )


    (gov_kwd_para,
     gov_cb_sw) = kwd_para
    
    # Inputs :  gen_ω, ω_ref0, p_order0    
    # Output : xg2 === Tm === x[2]    
         
    (gen_ω,
     ω_ref0,
     p_order0)   =
        gov_gen_ω_ωref0_porder0_para

    
    (;Tc,   
     Ts,   
     p_max,
     p_min,
     R) =
         NamedTupleTools.select(
             gov_kwd_para,
             (:Tc,   
              :Ts,   
              :p_max,
              :p_min,
              :R))     
    
    # Psv === x[1]
    # TM  === x[2] 
    
    xg1, xg2 = x

    res[1] = (1/Ts) * (p_order0 - xg1 -
        (1 / R) * (gen_ω/ω_ref0 - 1)  ) -  dx[1]

    res[2] = (1/Tc) * (xg1 - xg2) -  dx[2]
    
    return nothing    

end


function ode_gov_fun__gov_t1_cb_sauer__func!(
    dx, x,
    gov_gen_ω_ωref0_porder0_para,
    t;
    kwd_para =
        gov_kwd_para_wt_gov_cb_sw )


    (gov_kwd_para,
     gov_cb_sw) = kwd_para
    
    # Inputs :  gen_ω, ω_ref0, p_order0    
    # Output : xg2 === Tm === x[2]    
         
    (gen_ω,
     ω_ref0,
     p_order0)   =
        gov_gen_ω_ωref0_porder0_para

    
    (;Tc,   
     Ts,   
     p_max,
     p_min,
     R) =
         NamedTupleTools.select(
             gov_kwd_para,
             (:Tc,   
              :Ts,   
              :p_max,
              :p_min,
              :R))     
    
    # Psv === x[1]
    # TM  === x[2] 
    
    xg1, xg2 = x

    dx[1] = (1/Ts) * (p_order0 - xg1 -
        (1 / R) * (gen_ω/ω_ref0 - 1)) 

    dx[2] = (1/Tc) * (xg1 - xg2)
    
    return nothing    

end


function dae_gov_fun__gov_t1_cb_sauer__func!(
    res,
    dx, x,
    gov_gen_ω_ωref0_porder0_para, 
    t;
    kwd_para =
        gov_kwd_para_wt_gov_cb_sw )

    (gov_kwd_para,
     gov_cb_sw) = kwd_para

    
    # Inputs :  gen_ω, ω_ref0, p_order0    
    # Output : xg2 === Tm === x[2]    
         
    (gen_ω,
     ω_ref0,
     p_order0)   =
        gov_gen_ω_ωref0_porder0_para

    
    (;Tc,   
     Ts,   
     p_max,
     p_min,
     R) =
         NamedTupleTools.select(
             gov_kwd_para,
             (:Tc,   
              :Ts,   
              :p_max,
              :p_min,
              :R))     
    
    # Psv === x[1]
    # TM  === x[2] 
    
    xg1, xg2 = x

    res[1] = (1/Ts) * (p_order0 - xg1 -
        (1 / R) * (gen_ω/ω_ref0 - 1)  ) -  dx[1]

    res[2] = (1/Tc) * (xg1 - xg2) -  dx[2]
    
    return nothing    

end


function gov_fun__gov_t1_cb_sauer__output(
    x,
    gov_gen_ω_ωref0_porder0_para;
    kwd_para =
        gov_kwd_para )

    # Output : xg2 === Tm === x[2]    
         
    (gen_ω,
     ω_ref0,
     p_order0)   =
        gov_gen_ω_ωref0_porder0_para

    
    (;Tc,   
     Ts,   
     p_max,
     p_min,
     R) =
         NamedTupleTools.select(
             kwd_para,
             (:Tc,   
              :Ts,   
              :p_max,
              :p_min,
              :R))     
    
    xg1, xg2 = x
    
    return xg2   

end



function gov_fun__gov_t1_cb_sauer__callback_paras(
    plant_state_idx_in_states,
    plant_states_syms;
    comp_paras =
        gov_kwd_para )

    
    (p_max, 
     p_min) =
         NamedTupleTools.select(
             comp_paras,
             (:p_max, 
              :p_min))
    
    lower_lim_windup_cb_para =
        make_a_component_state_callback_paras(    
            :lower_lim_windup,

            plant_state_idx_in_states,
            plant_states_syms,

            p_min,
            p_min,

            :xg1,
            :xg1,

            cb_fun_condition_lower_lim,
            cb_fun_affect_windup_lower_lim! )
    
    upper_lim_windup_cb_para =
        make_a_component_state_callback_paras(    
            :upper_lim_windup,

            plant_state_idx_in_states,
            plant_states_syms,

            p_max,
            p_max,

            :xg1,
            :xg1,

            cb_fun_condition_upper_lim,
            cb_fun_affect_windup_upper_lim! )

    
    return [ lower_lim_windup_cb_para,
             upper_lim_windup_cb_para]
   

end

#----------------------------------------
# Avr
#----------------------------------------


function avr_fun__avr_t0_cb__init(
    vh, vf_tilade;
    kwd_para =
        avr_kwd_para )

    (;Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max,
     V_R_min, Ae, Be) =
         NamedTupleTools.select(
             kwd_para,
             (:Ta, :Te, :Tf, :Tr, :Ka, :Ke, :Kf,
              :V_R_max, :V_R_min, :Ae, :Be))
    vm = vh
    
    vr1 = vf_tilade * (Ke + Sevf(Ae, Be, vf_tilade))
    
    vr2 = (Kf / Tf) * vf_tilade
    
    v_ref = (vr1 / Ka) + vm - vr2 + (Kf / Tf) * vf_tilade
      
    return (; state_var = Float64[
        vm, vr1, vr2, vf_tilade],
            ref = (;v_ref,))

end



function ode_avr_fun__avr_t0_cb__func!(
    dx, x,
    avr_vh_vref_para, 
    t;
    kwd_para =
        avr_kwd_para_wt_avr_cb_sw )


    (avr_kwd_para,
     avr_cb_sw) =
         kwd_para

    vh, v_ref0 = avr_vh_vref_para

    sw = avr_cb_sw[1][1]
    
    (;Ta,      
     Te,      
     Tf,      
     Tr,      
     Ka,      
     Ke,      
     Kf,      
     V_R_max, 
     V_R_min, 
     Ae,      
     Be) =
         NamedTupleTools.select(
             avr_kwd_para,
             (:Ta,      
              :Te,      
              :Tf,      
              :Tr,      
              :Ka,      
              :Ke,      
              :Kf,      
              :V_R_max, 
              :V_R_min, 
              :Ae,      
              :Be))

    # vm, vr1, vr2, vf_tilade,  vr1_hat = x
    
    vm, vr1, vr2, vf_tilade  = x   
     
    dx[1] = (1/Tr)  * (vh - vm )

    dx[2] = (Ka * ((v_ref0 - vm) + vr2 -
        (Kf/Tf) * vf_tilade ) - vr1) /Ta 
    
    dx[3] = -vr2/Tf + (Kf/Tf^2) * vf_tilade
    
    if sw == 0
        
    dx[4] = (-1/Te) * (vf_tilade * (Ke + Sevf(
        Ae, Be, vf_tilade)) - vr1)
        
    elseif sw == 1
        
        dx[4] = 0.0
        
    elseif sw == 2
        
        dx[4] = 0.0
    end
    
    return nothing    

end



function dae_avr_fun__avr_t0_cb__func!(
    res, dx, x, avr_vh_vref_para, t;
    kwd_para =
        avr_kwd_para_wt_avr_cb_sw )

    (avr_kwd_para,
     avr_cb_sw) =
         kwd_para

    sw = avr_cb_sw[1][1]

    vh, v_ref0 = avr_vh_vref_para

    (;Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max,
     V_R_min, Ae, Be) =
         NamedTupleTools.select(
             avr_kwd_para, (
                 :Ta, :Te, :Tf, :Tr, :Ka, :Ke, :Kf,
                 :V_R_max, :V_R_min, :Ae, :Be))
    
    vm, vr1, vr2, vf_tilade  = x   
     
    res[1] = (1/Tr) * (vh - vm ) - dx[1] 

    res[2] = (Ka * ((v_ref0 - vm) + vr2 -
        (Kf/Tf) * vf_tilade ) - vr1) /Ta -  dx[2]
    
    res[3] = -vr2/Tf + (Kf/Tf^2) * vf_tilade - dx[3]

    
    if sw == 0
        
        res[4] = (-1/Te) * (vf_tilade * (Ke + Sevf(
        Ae, Be, vf_tilade)) - vr1) - dx[4] 
        
    elseif sw == 1

        dx[4]  = 0.0
        res[4] = 0.0 -  dx[4]
        
        
    elseif sw == 2

        dx[4] = 0.0
        res[4] = 0.0 - dx[4]
        
    end
    
    return nothing    

end



function avr_fun__avr_t0_cb__func!(
    res, dx, x, avr_vh_vref_para, t;
    kwd_para =
        avr_kwd_para_wt_avr_cb_sw )

    (avr_kwd_para,
     avr_cb_sw) =
         kwd_para

    sw = avr_cb_sw[1][1]

    vh, v_ref0 = avr_vh_vref_para

    (;Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max,
     V_R_min, Ae, Be) = NamedTupleTools.select(
         avr_kwd_para, (
             :Ta, :Te, :Tf, :Tr, :Ka, :Ke, :Kf,
             :V_R_max, :V_R_min, :Ae, :Be))
    
    vm, vr1, vr2, vf_tilade  = x   
     
    res[1] = (1/Tr)  * (vh - vm ) - dx[1] 

    res[2] = (Ka * ((v_ref0 - vm) + vr2 -
        (Kf/Tf) * vf_tilade ) - vr1) /Ta -  dx[2]
    
    res[3] = -vr2/Tf + (Kf/Tf^2) * vf_tilade - dx[3]

    
    if sw == 0
        
        res[4] = (-1/Te) * (vf_tilade * (Ke + Sevf(
        Ae, Be, vf_tilade)) - vr1) - dx[4] 
        
    elseif sw == 1

        dx[4]  = 0.0
        res[4] = 0.0 -  dx[4]
        
        
    elseif sw == 2

        dx[4] = 0.0
        res[4] = 0.0 - dx[4]
        
    end
    
    return nothing    

end


function avr_fun__avr_t0_cb__output(
    x,
    avr_vh_vref_para;
    kwd_para =
        avr_kwd_para )

    vm, vr1, vr2, vf_tilade  = x   
      
    return vf_tilade  

end


function avr_fun__avr_t0_cb__callback_paras(
    plant_state_idx_in_states,
    plant_states_syms;
    comp_paras =
        avr_para )
    
    (V_R_max, 
     V_R_min) =
         NamedTupleTools.select(
             comp_paras,
             (:V_R_max, 
              :V_R_min))

    component_parameter_switch = [0]
    
    lower_lim_anti_windup_cb_para =
        make_a_component_state_callback_paras(    
            :lower_lim_anti_windup,

            plant_state_idx_in_states,
            plant_states_syms,

            V_R_min,
            V_R_min,

            :vr1,
            :vr1,

            cb_fun_condition_lower_lim,
            cb_fun_affect_anti_windup_lower_lim!;
            
            component_parameter_switch =
                component_parameter_switch)
    
    upper_lim_anti_windup_cb_para =
        make_a_component_state_callback_paras(    
            :upper_lim_anti_windup,

            plant_state_idx_in_states,
            plant_states_syms,

            V_R_max,
            V_R_max,

            :vr1,
            :vr1,

            cb_fun_condition_upper_lim,
            cb_fun_affect_anti_windup_upper_lim!;
            
            component_parameter_switch =
                component_parameter_switch )

    
    return [ lower_lim_anti_windup_cb_para,
             upper_lim_anti_windup_cb_para]
   

end

#---------------------------------------------------

function avr_fun__avr_t1_cb__init(
    vh,
    vf_tilade;
    kwd_para =
        avr_kwd_para )


    # Inputs:  vh, V_ref0    
    # Output:  vf_tilade === x[4]    
    
    (;Ta,      
     Te,      
     Tf,      
     Tr,      
     Ka,      
     Ke,      
     Kf,      
     V_R_max, 
     V_R_min, 
     Ae,      
     Be) =
         NamedTupleTools.select(
             kwd_para,
             (:Ta,      
              :Te,      
              :Tf,      
              :Tr,      
              :Ka,      
              :Ke,      
              :Kf,      
              :V_R_max, 
              :V_R_min, 
              :Ae,      
              :Be))

    
    vm = vh
    
    vr1 = vf_tilade * (Ke + Sevf(Ae, Be, vf_tilade))

    
    vr2 = (Kf / Tf) * vf_tilade
    
    v_ref = (vr1 / Ka) + vm - vr2 + (Kf / Tf) * vf_tilade

    
    return (; state_var =
        Float64[vm, vr1, vr2, vf_tilade],
            ref = (; v_ref,) )    

end



function ode_avr_fun__avr_t1_cb__func!(
    dx, x,
    avr_vh_vref_para, 
    t;
    kwd_para =
        avr_kwd_para_wt_avr_cb_sw )

    (avr_kwd_para,
     avr_cb_sw) =
         kwd_para

    sw = avr_cb_sw[1][1]
    
    # Inputs:  vh, V_ref0    
    # Output:  vf_tilade === x[4]    

    vh, v_ref0 = avr_vh_vref_para
    
    (;Ta,      
     Te,      
     Tf,      
     Tr,      
     Ka,      
     Ke,      
     Kf,      
     V_R_max, 
     V_R_min, 
     Ae,      
     Be) =
         NamedTupleTools.select(
             avr_kwd_para,
             (:Ta,      
              :Te,      
              :Tf,      
              :Tr,      
              :Ka,      
              :Ke,      
              :Kf,      
              :V_R_max, 
              :V_R_min, 
              :Ae,      
              :Be))
    
    vm, vr1, vr2, vf_tilade  = x   


    dx[1] = (1/Tr)  * (vh - vm )

    # dx[2] = (Ka * ((v_ref0 - vm) + vr2 -
    #     (Kf/Tf) * vf_tilade ) - vr1) /Ta
    
    if sw == 0
        
        dx[2] = (Ka * ((v_ref0 - vm) + vr2 -

            (Kf/Tf) * vf_tilade ) - vr1) /Ta
        
    elseif sw == 1
        
        dx[2] = 0.0
        
    elseif sw == 2
        
        dx[2] = 0.0
    end

    dx[3] = -vr2/Tf + (Kf/Tf^2) * vf_tilade

    dx[4] = (1/Te) * (vr1 - vf_tilade *
        (Ke + Sevf(Ae, Be, vf_tilade)) )

    
    return nothing    

end


function dae_avr_fun__avr_t1_cb__func!(
    res,
    dx, x,
    avr_vh_vref_para,
    t;
    kwd_para =
        avr_kwd_para_wt_avr_cb_sw )


    (avr_kwd_para,
     avr_cb_sw) =
         kwd_para

    sw = avr_cb_sw[1][1]

    # Inputs:  vh, V_ref0    
    # Output:  vf_tilade === x[4]    

    vh, v_ref0 = avr_vh_vref_para
    
    (;Ta,      
     Te,      
     Tf,      
     Tr,      
     Ka,      
     Ke,      
     Kf,      
     V_R_max, 
     V_R_min, 
     Ae,      
     Be) =
         NamedTupleTools.select(
             avr_kwd_para,
             (:Ta,      
              :Te,      
              :Tf,      
              :Tr,      
              :Ka,      
              :Ke,      
              :Kf,      
              :V_R_max, 
              :V_R_min, 
              :Ae,      
              :Be))
    
    vm, vr1, vr2, vf_tilade  = x   


    res[1] = (1/Tr)  * (vh - vm ) - dx[1] 
    
    if sw == 0
        
        res[2] = (Ka * ((v_ref0 - vm) + vr2 -
            (Kf/Tf) * vf_tilade ) - vr1) /Ta - dx[2]
        
    elseif sw == 1

        dx[2] = 0.0
        res[2] = 0.0 - dx[2]
        
    elseif sw == 2

        dx[2] = 0.0
        res[2] = 0.0 - dx[2]
    end

    # res[2] = (Ka * ((v_ref0 - vm) + vr2 -
    #     (Kf/Tf) * vf_tilade ) - vr1) /Ta - dx[2]

    res[3] = -vr2/Tf + (Kf/Tf^2) * vf_tilade - dx[3] 

    res[4] = (1/Te) * (vr1 - vf_tilade *
        (Ke + Sevf(Ae, Be, vf_tilade)) )  - dx[4]  

    
    return nothing    

end


function avr_fun__avr_t1_cb__func!(
    res,
    dx, x,
    avr_vh_vref_para, 
    t;
    kwd_para =
        avr_kwd_para_wt_avr_cb_sw )

    (avr_kwd_para,
     avr_cb_sw) =
         kwd_para

    sw = avr_cb_sw[1][1]
    

    # Inputs:  vh, V_ref0    
    # Output:  vf_tilade === x[4]    

    vh, v_ref0 = avr_vh_vref_para
    
    (;Ta,      
     Te,      
     Tf,      
     Tr,      
     Ka,      
     Ke,      
     Kf,      
     V_R_max, 
     V_R_min, 
     Ae,      
     Be) =
         NamedTupleTools.select(
             avr_kwd_para,
             (:Ta,      
              :Te,      
              :Tf,      
              :Tr,      
              :Ka,      
              :Ke,      
              :Kf,      
              :V_R_max, 
              :V_R_min, 
              :Ae,      
              :Be))
    
    vm, vr1, vr2, vf_tilade  = x   


    res[1] = (1/Tr)  * (vh - vm ) - dx[1] 
    
    if sw == 0
        
        res[2] = (Ka * ((v_ref0 - vm) + vr2 -
            (Kf/Tf) * vf_tilade ) - vr1) /Ta - dx[2]
        
    elseif sw == 1

        dx[2] = 0.0
        res[2] = 0.0 - dx[2]
        
    elseif sw == 2

        dx[2] = 0.0
        res[2] = 0.0 - dx[2]
    end

    # res[2] = (Ka * ((v_ref0 - vm) + vr2 -
    #     (Kf/Tf) * vf_tilade ) - vr1) /Ta - dx[2]

    res[3] = -vr2/Tf + (Kf/Tf^2) * vf_tilade - dx[3] 

    res[4] = (1/Te) * (vr1 - vf_tilade *
        (Ke + Sevf(Ae, Be, vf_tilade)) )  - dx[4]  

    
    return nothing    

end


function avr_fun__avr_t1_cb__output(
    x,
    avr_vh_vref_para;
    kwd_para =
        avr_kwd_para )


    # Inputs:  vh, V_ref0    
    # Output:  vf_tilade === x[4]    
    
    vm, vr1, vr2, vf_tilade  = x   
    
    return vf_tilade  

end


function avr_fun__avr_t1_cb__callback_paras(
    plant_state_idx_in_states,
    plant_states_syms ;
    comp_paras =
        avr_para )

    
    (V_R_max, 
     V_R_min) =
         NamedTupleTools.select(
             comp_paras,
             (:V_R_max, 
              :V_R_min))

    component_parameter_switch = [0]
    
    lower_lim_anti_windup_cb_para =
        make_a_component_state_callback_paras(    
            :lower_lim_anti_windup,

            plant_state_idx_in_states,
            plant_states_syms,

            V_R_min,
            V_R_min,

            :vr1,
            :vr1,

            cb_fun_condition_lower_lim,
            cb_fun_affect_anti_windup_lower_lim!;
            
            component_parameter_switch =
                component_parameter_switch  )
    
    upper_lim_anti_windup_cb_para =
        make_a_component_state_callback_paras(    
            :upper_lim_anti_windup,

            plant_state_idx_in_states,
            plant_states_syms,

            V_R_max,
            V_R_max,

            :vr1,
            :vr1,

            cb_fun_condition_upper_lim,
            cb_fun_affect_anti_windup_upper_lim!;
            
            component_parameter_switch =
                component_parameter_switch  )

    
    return [ lower_lim_anti_windup_cb_para,
             upper_lim_anti_windup_cb_para]
   

end

#---------------------------------------------------


function avr_fun__avr_t1_cb_sauer__init(
   vh, vf_tilade ;
    kwd_para =
        avr_kwd_para )

    
    (;Ta,      
     Te,      
     Tf,      
     Tr,      
     Ka,      
     Ke,      
     Kf,      
     V_R_max, 
     V_R_min, 
     Ae,      
     Be) =
         NamedTupleTools.select(
             kwd_para,
             (:Ta,      
              :Te,      
              :Tf,      
              :Tr,      
              :Ka,      
              :Ke,      
              :Kf,      
              :V_R_max, 
              :V_R_min, 
              :Ae,      
              :Be))

        
    vr1 = vf_tilade * (Ke + Sevf(Ae, Be, vf_tilade))
    
    vr2 = (Kf / Tf) * vf_tilade
    
    v_ref = vh + (vr1 / Ka)

    return (; state_var = Float64[vr1, vr2, vf_tilade],
            ref = (;v_ref,))
    
    return nothing    

end


function ode_avr_fun__avr_t1_cb_sauer__func!(
    dx, x,
    avr_vh_vref_para, 
    t;
    kwd_para =
        avr_kwd_para_wt_avr_cb_sw )


    (avr_kwd_para,
     avr_cb_sw) =
         kwd_para
    
    sw = avr_cb_sw[1][1]
    

    # Inputs:  vh, V_ref0    
    # Output:  vf_tilade === x[3]    

    vh, v_ref0 = avr_vh_vref_para

    
    (;Ta,      
     Te,      
     Tf,      
     Tr,      
     Ka,      
     Ke,      
     Kf,      
     V_R_max, 
     V_R_min, 
     Ae,      
     Be) =
         NamedTupleTools.select(
             avr_kwd_para,
             (:Ta,      
              :Te,      
              :Tf,      
              :Tr,      
              :Ka,      
              :Ke,      
              :Kf,      
              :V_R_max, 
              :V_R_min, 
              :Ae,      
              :Be))
    
    vr1,  vr2, vf_tilade  = x   



    # dx[1] = ( Ka * vr2 - (Ka * Kf/Tf) *
    #         vf_tilade  + Ka * (v_ref0 - vh) - vr1 ) /Ta
    
    if sw == 0

        dx[1] = (-vr1 + Ka * vr2 - (Ka * Kf/Tf) *
            vf_tilade  + Ka * (v_ref0 - vh)  ) /Ta
                
    elseif sw == 1
        
        dx[1] = 0.0
        
    elseif sw == 2
        
        dx[1] = 0.0
    end

    dx[2] = (1/Tf) * ( (Kf/Tf) * vf_tilade - vr2  ) 

    dx[3] = (1/Te) * ( vr1 - vf_tilade * (Ke +
        Sevf(Ae, Be, vf_tilade) ) )
        
    return nothing    

end


function dae_avr_fun__avr_t1_cb_sauer__func!(
    res,
    dx, x,
    avr_vh_vref_para,
    t;
    kwd_para =
        avr_kwd_para_wt_avr_cb_sw )


    (avr_kwd_para,
     avr_cb_sw) =
         kwd_para
    
    sw = avr_cb_sw[1][1]
    
    # Inputs:  vh, V_ref0    
    # Output:  vf_tilade === x[3]    

    vh, v_ref0 = avr_vh_vref_para

    
    (;Ta,      
     Te,      
     Tf,      
     Tr,      
     Ka,      
     Ke,      
     Kf,      
     V_R_max, 
     V_R_min, 
     Ae,      
     Be) =
         NamedTupleTools.select(
             avr_kwd_para,
             (:Ta,      
              :Te,      
              :Tf,      
              :Tr,      
              :Ka,      
              :Ke,      
              :Kf,      
              :V_R_max, 
              :V_R_min, 
              :Ae,      
              :Be))
    
    vr1, vr2, vf_tilade  = x   
    
    if sw == 0

        res[1] = (-vr1 + Ka * vr2 - (Ka * Kf/Tf) *
            vf_tilade  + Ka * (v_ref0 - vh)) /Ta - dx[1]
                
    elseif sw == 1

        dx[1] = 0.0
        res[1] = 0.0 - dx[1]
        
    elseif sw == 2

        dx[1] = 0.0
        res[1] = 0.0 - dx[1]
    end


    # res[1] = ( Ka * vr2 - (Ka * Kf/Tf) *
    #  vf_tilade + Ka * (v_ref0 - vh) - vr1) /Ta - dx[1]

    res[2] = (1/Tf) * ((Kf/Tf) * vf_tilade - vr2 ) - dx[2]

    res[3] = (1/Te) * ( vr1 - vf_tilade * (Ke +
        Sevf(Ae, Be, vf_tilade) ) ) - dx[3]   
        
    return nothing    

end


function avr_fun__avr_t1_cb_sauer__func!(
    res,
    dx, x,
    avr_vh_vref_para,
    t;
    kwd_para =
        avr_kwd_para_wt_avr_cb_sw )


    (avr_kwd_para,
     avr_cb_sw) =
         kwd_para
    
    sw = avr_cb_sw[1][1]

    # Inputs:  vh, V_ref0    
    # Output:  vf_tilade === x[3]    

    vh, v_ref0 = avr_vh_vref_para

    
    (;Ta,      
     Te,      
     Tf,      
     Tr,      
     Ka,      
     Ke,      
     Kf,      
     V_R_max, 
     V_R_min, 
     Ae,      
     Be) =
         NamedTupleTools.select(
             avr_kwd_para,
             (:Ta,      
              :Te,      
              :Tf,      
              :Tr,      
              :Ka,      
              :Ke,      
              :Kf,      
              :V_R_max, 
              :V_R_min, 
              :Ae,      
              :Be))
    
    vr1,  vr2, vf_tilade  = x   
    
    if sw == 0

        res[1] = (-vr1 + Ka * vr2 - (Ka * Kf/Tf) *
            vf_tilade  + Ka * (v_ref0 - vh)) /Ta -  dx[1]
                
    elseif sw == 1

        dx[1] = 0.0
        res[1] = 0.0 - dx[1]
        
    elseif sw == 2

        dx[1] = 0.0
        res[1] = 0.0 - dx[1]
    end


    # res[1] = ( Ka * vr2 - (Ka * Kf/Tf) *
    # vf_tilade + Ka * (v_ref0 - vh) - vr1 )/Ta - dx[1]

   res[2] = (1/Tf) * ( (Kf/Tf) * vf_tilade - vr2 ) - dx[2]

   res[3] = (1/Te) * ( vr1 - vf_tilade * (Ke +
        Sevf(Ae, Be, vf_tilade) ) ) - dx[3]   
        
    return nothing    

end



function avr_fun__avr_t1_cb_sauer__output(
    x,
    avr_vh_vref_para;
    kwd_para =
        avr_kwd_para )

    
    vr1,  vr2, vf_tilade  = x   
        
    return vf_tilade    

end


function avr_fun__avr_t1_cb_sauer__callback_paras(
    plant_state_idx_in_states,
    plant_states_syms ;
    comp_paras =
        avr_para )

    
    (V_R_max, 
     V_R_min) =
         NamedTupleTools.select(
             comp_paras,
             (:V_R_max, 
              :V_R_min))

    component_parameter_switch = [0]
    
    lower_lim_anti_windup_cb_para =
        make_a_component_state_callback_paras(    
            :lower_lim_anti_windup,

            plant_state_idx_in_states,
            plant_states_syms,

            V_R_min,
            V_R_min,

            :vr1,
            :vr1,

            cb_fun_condition_lower_lim,
            cb_fun_affect_anti_windup_lower_lim!;
            
            component_parameter_switch =
                component_parameter_switch )
    
    upper_lim_anti_windup_cb_para =
        make_a_component_state_callback_paras(    
            :upper_lim_anti_windup,

            plant_state_idx_in_states,
            plant_states_syms,

            V_R_max,
            V_R_max,

            :vr1,
            :vr1,

            cb_fun_condition_upper_lim,
            cb_fun_affect_anti_windup_upper_lim!;
            
            component_parameter_switch =
                component_parameter_switch )

    
    return [ lower_lim_anti_windup_cb_para,
             upper_lim_anti_windup_cb_para]
   

end


#---------------------------------------------------
#---------------------------------------------------
# Special functions
#---------------------------------------------------
#---------------------------------------------------


"It is assumed that all SM machines have the same
characteristics. The same assumption holds for SC.
Consequently, only one function is created for
each type"


#----------------------------------------
# Generators with local loads
#----------------------------------------


function gen_fun__SM_2axis_wt_loc_load_cb_v6__init(
    gen_vh,
    gen_θh,
    P_g,
    Q_g;
    kwd_para = gen_para)
    
        return gen_fun__SM_2axis_cb_v6__init(
            gen_vh,
            gen_θh,
            P_g,
            Q_g;
            kwd_para = kwd_para)

end


function ode_gen_fun__SM_2axis_wt_loc_load_cb_v6__func!(
    dx,
    x,
    gen_id_iq_τmtilade_vftilade_para,
    t;
    kwd_para =
        gen_kwd_para )
    
       return  ode_gen_fun__SM_2axis_cb_v6__func!(
            dx,
            x,
            gen_id_iq_τmtilade_vftilade_para,
            t;
            kwd_para =
                kwd_para )

end



function dae_gen_fun__SM_2axis_wt_loc_load_cb_v6__func!(
    res,
    dx,
    x,
    gen_id_iq_τmtilade_vftilade_para,
    t;
    kwd_para =
        gen_kwd_para )
    
       return  dae_gen_fun__SM_2axis_cb_v6__func!(
            res,
            dx,
            x,
            gen_id_iq_τmtilade_vftilade_para,
            t;
            kwd_para =
                kwd_para )

end



function gen_fun__SM_2axis_wt_loc_load_cb_v6__func!(
    res,
    dx,
    x,
    gen_id_iq_τmtilade_vftilade_para,
    t;
    kwd_para =
        gen_kwd_para )
    
       return  gen_fun__SM_2axis_cb_v6__func!(
            res,
            dx,
            x,
            gen_id_iq_τmtilade_vftilade_para,
            t;
            kwd_para =
                kwd_para )

end


function gen_fun__SM_2axis_wt_loc_load_cb_v6__output(
    x,
    ωref0_vref0_porder0_id_iq_vh;
    kwd_para =
        gen_kwd_para )
    
        return gen_fun__SM_2axis_cb_v6__output(
            x,
            ωref0_vref0_porder0_id_iq_vh;
            kwd_para =
                kwd_para )

end


function gen_fun__SM_2axis_wt_loc_load_cb_v6__callback_paras(
    plant_state_idx_in_states,
    plant_states_syms ;
    comp_paras =
        gen_kwd_para )
    
        return gen_fun__SM_2axis_cb_v6__callback_paras(
    plant_state_idx_in_states,
            plant_states_syms ;
            comp_paras =
        comp_paras )

end

#---------------------------------------------------

function gen_fun__SC_2axis_wt_loc_load_cb_v6__init(
    gen_vh,
    gen_θh,
    P_g,
    Q_g;
    kwd_para = gen_para)
    
        return gen_fun__SC_2axis_cb_v6__init(
            gen_vh,
            gen_θh,
            P_g,
            Q_g;
            kwd_para =
                kwd_para)

end



function ode_gen_fun__SC_2axis_wt_loc_load_cb_v6__func!(
    dx,
    x,
    gen_id_iq_τmtilade_vftilade_para,
    t;
    kwd_para =
        gen_kwd_para ) 
        return ode_gen_fun__SC_2axis_cb_v6__func!(
            dx,
            x,
            gen_id_iq_τmtilade_vftilade_para,
            t;
            kwd_para =
                kwd_para )

end


function dae_gen_fun__SC_2axis_wt_loc_load_cb_v6__func!(
    res,
    dx,
    x,
    gen_id_iq_τmtilade_vftilade_para,
    t;
    kwd_para =
        gen_kwd_para ) 
        return dae_gen_fun__SC_2axis_cb_v6__func!(
            res,
            dx,
            x,
            gen_id_iq_τmtilade_vftilade_para,
            t;
            kwd_para =
                kwd_para )

end

function gen_fun__SC_2axis_wt_loc_load_cb_v6__func!(
    res,
    dx,
    x,
    gen_id_iq_τmtilade_vftilade_para,
    t;
    kwd_para =
        gen_kwd_para ) 
        return gen_fun__SC_2axis_cb_v6__func!(
            res,
            dx,
            x,
            gen_id_iq_τmtilade_vftilade_para,
            t;
            kwd_para =
                kwd_para )

end


function gen_fun__SC_2axis_wt_loc_load_cb_v6__output(
    x,
    ωref0_vref0_porder0_id_iq_vh;
    kwd_para =
        gen_kwd_para )
    
        return gen_fun__SC_2axis_cb_v6__output(
            x,
            ωref0_vref0_porder0_id_iq_vh;
            kwd_para =
                kwd_para )

end


function gen_fun__SC_2axis_wt_loc_load_cb_v6__callback_paras(
    plant_state_idx_in_states,
    plant_states_syms ;
    comp_paras =
        gen_kwd_para )
    
        return gen_fun__SC_2axis_cb_v6__callback_paras(
    plant_state_idx_in_states,
            plant_states_syms ;
            comp_paras =
        comp_paras )

end


#----------------------------------------
# Gov special functions
#----------------------------------------


"These functions are meant for sycnronous machines.
They do not have a governor hence an empty anonymus
function is created for them.

The name of the function is used to test if a gen plant is a
synchronous macheine

`nameof(comp_dyn_fun) == :gov_fun__nothing__func!`
"


function gov_fun__nothing__callback_paras()

    return [make_a_component_state_no_callback_paras(
        :nothing)]
end

gov_fun__nothing__init() = anonymus_func 

gov_fun__nothing__output() = anonymus_func 

gov_fun__nothing__func!() = anonymus_func 

ode_gov_fun__nothing__func!() = anonymus_func 

dae_gov_fun__nothing__func!() = anonymus_func 


#---------------------------------------------------
#---------------------------------------------------
# dictionaries of comp dyn of output functions 
#---------------------------------------------------
#---------------------------------------------------

"""
Dictionaries are created to hold init, dynamic and output
functons. The keys of the dictionaries are symbols of
components names

The design pattern for naming functions is:

component_fun__componentsubtype__functiontype.

The rationale for the design pattern is to split the name
of each function by "__". This will give three substrings
 "component_fun",  "componentsubtype", and  "functiontype".

The second substring is the name of a specific component, hence is is used as a key in a dictionary.

`Dict( componentsubtype =>
       component_fun__componentsubtype__functiontype )`


Addind a new type of device requires, creating an
init function, dynamic function (state equatons) and
output function. onces these functions are created, they should be added to the list in `get_dict_gens_type_init_func`,
`get_dict_gens_type_dyn_func` and
`get_dict_gens_type_output_func` respectively.


"""


function get_dict_gens_type_callback_paras_func()
    
    types_func =
        Function[gen_fun__SM_2axis_cb_v6__callback_paras,
                 gen_fun__SC_2axis_cb_v6__callback_paras,
                 gen_fun__SM_2axis_wt_loc_load_cb_v6__callback_paras,
                 gen_fun__SC_2axis_wt_loc_load_cb_v6__callback_paras ]

    return get_dict_types_dyn_or_output_func(
        types_func)
        
end


function get_dict_gens_type_init_func()


    # :SM_2axis_cb_v6
    # :SM_2axis_wt_loc_load_cb_v6
    
    # :SC_2axis_cb_v6
    # :SC_2axis_wt_loc_load_cb_v6
    
    types_func =
        Function[gen_fun__SM_2axis_cb_v6__init,
                 gen_fun__SC_2axis_cb_v6__init,
                gen_fun__SM_2axis_wt_loc_load_cb_v6__init,
                gen_fun__SC_2axis_wt_loc_load_cb_v6__init]

    return get_dict_types_dyn_or_output_func(
        types_func)
        
end


function get_dict_gens_type_output_func()


    # :SM_2axis_cb_v6
    # :SM_2axis_wt_loc_load_cb_v6
    
    # :SC_2axis_cb_v6
    # :SC_2axis_wt_loc_load_cb_v6
    
    types_func =
        Function[gen_fun__SM_2axis_cb_v6__output,
                 gen_fun__SC_2axis_cb_v6__output,
                 gen_fun__SM_2axis_wt_loc_load_cb_v6__output,
                 gen_fun__SC_2axis_wt_loc_load_cb_v6__output]

    return get_dict_types_dyn_or_output_func(
    types_func)
        
end


function get_namedtuple_dict_gens_type_dyn_func()
    
    ode_types_func =
        Function[ode_gen_fun__SM_2axis_cb_v6__func!,
                 ode_gen_fun__SC_2axis_cb_v6__func!,
                 ode_gen_fun__SM_2axis_wt_loc_load_cb_v6__func!,
                 ode_gen_fun__SC_2axis_wt_loc_load_cb_v6__func!]
    
    dae_types_func =
        Function[dae_gen_fun__SM_2axis_cb_v6__func!,
                 dae_gen_fun__SC_2axis_cb_v6__func!,
                 dae_gen_fun__SM_2axis_wt_loc_load_cb_v6__func!,
                 dae_gen_fun__SC_2axis_wt_loc_load_cb_v6__func!]
    
    return get_namedtuple_dict_types_dyn_or_output_func(
        ode_types_func,
        dae_types_func;
        ode_type = :ode,
        dae_type = :dae )
        
end


function get_dict_gens_type_dyn_func()


    # :SM_2axis_cb_v6
    # :SM_2axis_wt_loc_load_cb_v6
    
    # :SC_2axis_cb_v6
    # :SC_2axis_wt_loc_load_cb_v6
    
    types_func =
        Function[gen_fun__SM_2axis_cb_v6__func!,
                 gen_fun__SC_2axis_cb_v6__func!,
                 gen_fun__SM_2axis_wt_loc_load_cb_v6__func!,
                 gen_fun__SC_2axis_wt_loc_load_cb_v6__func!]
    
    return get_dict_types_dyn_or_output_func(
        types_func)
        
end


function get_dict_gens_type_dyn_func(
    dae_or_ode_type)

    namedtuple_dict_type =
        get_namedtuple_dict_gens_type_dyn_func()

    if dae_or_ode_type == :ode
        
        return namedtuple_dict_type.ode
        
    elseif dae_or_ode_type == :dae

        return namedtuple_dict_type.dae
    else

        throw("The model function type $(dae_or_ode_type) is not known. It should be :ode or :dae")
        
    end
    
        
end


#-----------------------------------------------------


function get_dict_govs_type_callback_paras_func()

    types_func =
        Function[gov_fun__gov_ieee_tgov1_cb__callback_paras,
                 gov_fun__gov_t1_cb__output,
                 gov_fun__gov_t1_cb_sauer__callback_paras,
                 gov_fun__nothing__callback_paras]

    return get_dict_types_dyn_or_output_func(
        types_func)
        
end


function get_dict_govs_type_init_func()

    types_func =
        Function[gov_fun__gov_ieee_tgov1_cb__init,
                 gov_fun__gov_t1_cb__init,
                 gov_fun__gov_t1_cb_sauer__init,
                 gov_fun__nothing__init]

    return get_dict_types_dyn_or_output_func(
        types_func)
        
end



function get_dict_govs_type_output_func()

    types_func =
        Function[gov_fun__gov_ieee_tgov1_cb__output,
                 gov_fun__gov_t1_cb__output,
                 gov_fun__gov_t1_cb_sauer__output,
                 gov_fun__nothing__output]

    return get_dict_types_dyn_or_output_func(
        types_func)
        
end



function get_namedtuple_dict_govs_type_dyn_func()
    
    ode_types_func =
        Function[ode_gov_fun__gov_ieee_tgov1_cb__func!,
                 ode_gov_fun__gov_t1_cb__func!,
                 ode_gov_fun__gov_t1_cb_sauer__func!,
                 ode_gov_fun__nothing__func!]
    
    dae_types_func =
        Function[dae_gov_fun__gov_ieee_tgov1_cb__func!,
                 dae_gov_fun__gov_t1_cb__func!,
                 dae_gov_fun__gov_t1_cb_sauer__func!,
                 dae_gov_fun__nothing__func!]
    
    return get_namedtuple_dict_types_dyn_or_output_func(
        ode_types_func,
        dae_types_func;
        ode_type = :ode,
        dae_type = :dae )
        
end


function get_dict_govs_type_dyn_func()

    types_func =
        Function[gov_fun__gov_ieee_tgov1_cb__func!,
                 gov_fun__gov_t1_cb__func!,
                 gov_fun__gov_t1_cb_sauer__func!,
                 gov_fun__nothing__func!]

    return get_dict_types_dyn_or_output_func(
        types_func)
        
end


function get_dict_govs_type_dyn_func(
    dae_or_ode_type)

    namedtuple_dict_type =
        get_namedtuple_dict_govs_type_dyn_func()

    if dae_or_ode_type == :ode
        
        return namedtuple_dict_type.ode
        
    elseif dae_or_ode_type == :dae

        return namedtuple_dict_type.dae
    else

        throw("The model function type $(dae_or_ode_type) is not known. It should be :ode or :dae")
        
    end
    
        
end


#-----------------------------------------------------


function get_dict_avrs_type_callback_paras_func()

    types_func =
        Function[avr_fun__avr_t0_cb__callback_paras,
                 avr_fun__avr_t1_cb__callback_paras,
                 avr_fun__avr_t1_cb_sauer__callback_paras]

    return get_dict_types_dyn_or_output_func(
        types_func)
        
end


function get_dict_avrs_type_init_func()

    types_func =
        Function[avr_fun__avr_t0_cb__init,
                 avr_fun__avr_t1_cb__init,
                 avr_fun__avr_t1_cb_sauer__init]

    return get_dict_types_dyn_or_output_func(
        types_func)
        
end


function get_dict_avrs_type_output_func()

    types_func =
        Function[avr_fun__avr_t0_cb__output,
                 avr_fun__avr_t1_cb__output,
                 avr_fun__avr_t1_cb_sauer__output]

    return get_dict_types_dyn_or_output_func(
        types_func)
        
end



function get_namedtuple_dict_avrs_type_dyn_func()
    
    ode_types_func =
        Function[ode_avr_fun__avr_t0_cb__func!,
                 ode_avr_fun__avr_t1_cb__func!,
                 ode_avr_fun__avr_t1_cb_sauer__func!]
    
    dae_types_func =
        Function[dae_avr_fun__avr_t0_cb__func!,
                 dae_avr_fun__avr_t1_cb__func!,
                 dae_avr_fun__avr_t1_cb_sauer__func!]
    
    return get_namedtuple_dict_types_dyn_or_output_func(
        ode_types_func,
        dae_types_func;
        ode_type = :ode,
        dae_type = :dae )
        
end


function get_dict_avrs_type_dyn_func()

    types_func =
        Function[avr_fun__avr_t0_cb__func!,
                 avr_fun__avr_t1_cb__func!,
                 avr_fun__avr_t1_cb_sauer__func!]

    return get_dict_types_dyn_or_output_func(
        types_func)
        
end


function get_dict_avrs_type_dyn_func(
    dae_or_ode_type)

    namedtuple_dict_type =
        get_namedtuple_dict_avrs_type_dyn_func()

    if dae_or_ode_type == :ode
        
        return namedtuple_dict_type.ode
        
    elseif dae_or_ode_type == :dae

        return namedtuple_dict_type.dae
    else

        throw("The model function type $(dae_or_ode_type) is not known. It should be :ode or :dae")
        
    end
    
        
end

#----------------------------------------
# Initialisation: Synchronous condenser
#----------------------------------------

function a_SC_plant_generic_model_init_func(
    gen_vh,
    gen_θh,
    P_g,
    Q_g,
    ωs;
    kwd_para =
        plant_init_kwd_para )

    (gen_init_para,
     avr_init_para,
     comps_init_fun,
     ) =
         kwd_para
             
    
    (gen_init_fun,
     avr_init_fun) =
         comps_init_fun
    
    # #----------------------------------------

    (;δ,    
     ed_dash,
     eq_dash,
     τm_tilade,
     vf_tilade,
     i_d,
     i_q) =
         NamedTupleTools.select(
             gen_init_fun(
                 gen_vh,
                 gen_θh,
                 P_g,
                 Q_g;
                 kwd_para =
                     gen_init_para ),
             (:δ,    
              :ed_dash,
              :eq_dash,
              :τm_tilade,
              :vf_tilade,
              :i_d,
              :i_q))

    gen_states_init =
        [δ, ωs, ed_dash, eq_dash]
    
    #----------------------------------------

    (avr_states_init,
     avr_ref ) =
         NamedTupleTools.select(
             avr_init_fun(
                 gen_vh,
                 vf_tilade;
                 kwd_para =
                     avr_init_para ),
             (:state_var,
              :ref))
    
    #----------------------------------------

    plant_states_init = Float64[
        gen_states_init;
        avr_states_init ]

    #----------------------------------------

    p_order = 0
    
    plant_ref =
        Float64[ωs,
         ωs,
         avr_ref.v_ref,
         p_order,
         i_d,
         i_q,
         gen_vh]
    
    #----------------------------------------

    return (;plant_states_init,
            plant_ref )
    
end

#-----------------------------------------------------
# Initialisation: Two-axis Synchronous generator
#-----------------------------------------------------

function a_SM_plant_generic_model_init_func(
    gen_vh,
    gen_θh,
    P_g,
    Q_g,
    ωs;
    kwd_para =
        plant_init_kwd_para )

    (gen_init_para,
     avr_init_para,
     gov_init_para,
     comps_init_fun,
     ) =
         kwd_para
             
    
    (gen_init_fun,
     avr_init_fun,
     gov_init_fun) =
         comps_init_fun
    
    #----------------------------------------

    (;δ,    
     ed_dash,
     eq_dash,
     τm_tilade,
     vf_tilade,
     i_d,
     i_q) =
         NamedTupleTools.select(
             gen_init_fun(
                 gen_vh,
                 gen_θh,
                 P_g,
                 Q_g;
                 kwd_para =
                     gen_init_para ),
             (:δ,    
              :ed_dash,
              :eq_dash,
              :τm_tilade,
              :vf_tilade,
              :i_d,
              :i_q))

    gen_states_init =
        [δ, ωs, ed_dash, eq_dash]
    
    #----------------------------------------

    (avr_states_init,
     avr_ref ) =
         NamedTupleTools.select(
             avr_init_fun(
                 gen_vh,
                 vf_tilade;
                 kwd_para =
                     avr_init_para ),
             (:state_var,
              :ref))
    
    #----------------------------------------

    (gov_states_init,
     gov_ref ) =
         NamedTupleTools.select(
             gov_init_fun(
                 ωs,
                 τm_tilade;
                 kwd_para =
                     gov_init_para),
             (:state_var,
              :ref))

    #----------------------------------------

    plant_states_init = Float64[
        gen_states_init;
        avr_states_init;
        gov_states_init ]

    #----------------------------------------

    plant_ref =
        Float64[ωs,
         gov_ref.ω_ref0,
         avr_ref.v_ref,
         gov_ref.p_order,
         i_d,
         i_q,
         gen_vh]
    
    #----------------------------------------

    return (;plant_states_init, plant_ref )
end

#----------------------------------------
# Dynamics: Two-axis Synchronous generator
#----------------------------------------

function ode_test_a_SM_plant_generic_model_func!(
    dx,
    x,
    a_ωref0_vref0_porder0_id_iq_vh,
    t;
    kwd_para =
        plant_kwd_para )

    (gen_para,
     avr_para,
     gov_para,
     
     comps_Idx,
     
     # gov_avr_cb_sw,
     avr_gov_cb_sw,
     avr_gov_cb_sw_Idx,     
     
     comps_dyn_fun,
     comps_output_fun,
     
     ωs
     ) =
         kwd_para
             
    

    (gen_Idx,
     avr_Idx,
     gov_Idx
     ) =
         comps_Idx
    
    (gen_fun!,
     avr_fun!,
     gov_fun! ) =
         comps_dyn_fun

    
    (gen_output_fun,
     avr_output_fun,
     gov_output_fun) =
         comps_output_fun


    (avr_cb_sw_Idx,
     gov_cb_sw_Idx) =
         NamedTupleTools.select(
             avr_gov_cb_sw_Idx,
             (:avr_cb_para_sw_idx_in_plant,
              :gov_cb_para_sw_idx_in_plant))

    avr_cb_sw = avr_gov_cb_sw[avr_cb_sw_Idx ]
    
    gov_cb_sw = avr_gov_cb_sw[gov_cb_sw_Idx]

    avr_kwd_para_wt_avr_cb_sw =
        (avr_para, avr_cb_sw )
    
    gov_kwd_para_wt_gov_cb_sw =
        (gov_para, gov_cb_sw)
    
    #----------------------------------------

    states_dims  = length.(comps_Idx)
    state_offset = create_offsets(states_dims)
    state_Idx    = create_idxs(state_offset, states_dims)
    
    u_gen  = @view x[state_Idx[1]]
    u_avr  = @view x[state_Idx[2]]    
    u_gov  = @view x[state_Idx[3]]    
    
    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]
    du_gov = @view dx[state_Idx[3]]
    
    
    #----------------------------------------

    ωref0, vref0, porder0, id, iq, vh =
        a_ωref0_vref0_porder0_id_iq_vh

    #----------------------------------------

    avr_vh_Vref0_para =
        [vh, vref0]
    
    avr_fun!(
        du_avr,
        u_avr,
        avr_vh_Vref0_para,
        t;
        kwd_para =
            avr_kwd_para_wt_avr_cb_sw )
    
    avr_vf_tilade = avr_output_fun(
        u_avr,
        avr_vh_Vref0_para;
        kwd_para =
            avr_para )

    #----------------------------------------
    
    gen_ω = gen_output_fun(
        u_gen,
        a_ωref0_vref0_porder0_id_iq_vh;
        kwd_para =
            gen_para  )

    gov_gen_ω_ωref0_porder0_para =
        [gen_ω, ωref0, porder0, ]

    
    gov_fun!(
        du_gov,
        u_gov,
        gov_gen_ω_ωref0_porder0_para, 
        t;
        kwd_para =
            gov_kwd_para_wt_gov_cb_sw )

    
    gov_τm_tilade = gov_output_fun(
        u_gov,
        gov_gen_ω_ωref0_porder0_para;
        kwd_para =
            gov_para )
    
    #----------------------------------------
    
    gen_id_iq_τmtilade_vftilade_para =
        [id, iq, gov_τm_tilade, avr_vf_tilade]
    
    gen_fun!(
        du_gen,
        u_gen,
        gen_id_iq_τmtilade_vftilade_para,
        t;
        kwd_para =
            (gen_para, ωs) )

    #----------------------------------------

    return nothing
end



function ode_a_SM_plant_generic_model_func!(
    dx,
    x,
    a_ωref0_vref0_porder0_id_iq_vh,
    t;
    kwd_para =
        plant_kwd_para )

    (gen_para,
     avr_para,
     gov_para,
     
     comps_Idx,
     
     avr_gov_cb_sw,
     avr_gov_cb_sw_Idx,     
     
     comps_dyn_fun,
     comps_output_fun,
     
     ωs
     ) =
         kwd_para
             
    (gen_Idx,
     avr_Idx,
     gov_Idx
     ) =
         comps_Idx
    
    (gen_fun!,
     avr_fun!,
     gov_fun! ) =
         comps_dyn_fun

    
    (gen_output_fun,
     avr_output_fun,
     gov_output_fun) =
         comps_output_fun


    (avr_cb_sw_Idx,
     gov_cb_sw_Idx) =
         NamedTupleTools.select(
             avr_gov_cb_sw_Idx,
             (:avr_cb_para_sw_idx_in_plant,
              :gov_cb_para_sw_idx_in_plant))

    avr_cb_sw = avr_gov_cb_sw[avr_cb_sw_Idx ]
    
    gov_cb_sw = avr_gov_cb_sw[gov_cb_sw_Idx]

    avr_kwd_para_wt_avr_cb_sw =
        (avr_para, avr_cb_sw )
    
    gov_kwd_para_wt_gov_cb_sw =
        (gov_para, gov_cb_sw)
    
    #----------------------------------------
    #----------------------------------------

    states_dims  = length.(comps_Idx)
    state_offset = create_offsets(states_dims)
    state_Idx    = create_idxs(state_offset, states_dims)
    
    u_gen  = @view x[state_Idx[1]]
    u_avr  = @view x[state_Idx[2]]    
    u_gov  = @view x[state_Idx[3]]    
    
    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]
    du_gov = @view dx[state_Idx[3]]
    
    #----------------------------------------
    #----------------------------------------

    ωref0, vref0, porder0, id, iq, vh =
        a_ωref0_vref0_porder0_id_iq_vh

    #----------------------------------------

    avr_vh_Vref0_para =
        [vh, vref0]
    
    avr_fun!(
        du_avr,
        u_avr,
        avr_vh_Vref0_para,
        t;
        kwd_para =
            avr_kwd_para_wt_avr_cb_sw )
    
    avr_vf_tilade = avr_output_fun(
        u_avr,
        avr_vh_Vref0_para;
        kwd_para =
            avr_para )
    
    #----------------------------------------
    
    gen_ω = gen_output_fun(
        u_gen,
        a_ωref0_vref0_porder0_id_iq_vh;
        kwd_para =
            gen_para  )

    gov_gen_ω_ωref0_porder0_para =
        [gen_ω, ωref0, porder0, ]
    
    gov_fun!(
        du_gov,
        u_gov,
        gov_gen_ω_ωref0_porder0_para, 
        t;
        kwd_para =
            gov_kwd_para_wt_gov_cb_sw )

    
    gov_τm_tilade = gov_output_fun(
        u_gov,
        gov_gen_ω_ωref0_porder0_para;
        kwd_para =
            gov_para )
    
    #----------------------------------------
    
    gen_id_iq_τmtilade_vftilade_para =
        [id, iq, gov_τm_tilade, avr_vf_tilade]

    #  dae_a_SM_gen_generic_model_func!
    
    gen_fun!(
        du_gen,
        u_gen,
        gen_id_iq_τmtilade_vftilade_para,
        t;
        kwd_para =
            (gen_para, ωs) )
    
    #----------------------------------------

    return nothing
end


function dae_a_SM_plant_generic_model_func!(
    res,
    dx, x,
    a_ωref0_vref0_porder0_id_iq_vh,
    t;
    kwd_para =
        plant_kwd_para )

    (gen_para,
     avr_para,
     gov_para,
     
     comps_Idx,
     
     # gov_avr_cb_sw,
     avr_gov_cb_sw,
     avr_gov_cb_sw_Idx,
     
     comps_dyn_fun,
     comps_output_fun,
     
     ωs
     ) =
         kwd_para

    (gen_Idx,
     avr_Idx,
     gov_Idx
     ) =
         comps_Idx
    
    (gen_fun!,
     avr_fun!,
     gov_fun! ) =
         comps_dyn_fun

    
    (gen_output_fun,
     avr_output_fun,
     gov_output_fun) =
         comps_output_fun

    
    (gen_output_fun,
     avr_output_fun) =
         comps_output_fun

    (avr_cb_sw_Idx,
     gov_cb_sw_Idx) =
         NamedTupleTools.select(
             avr_gov_cb_sw_Idx,
             (:avr_cb_para_sw_idx_in_plant,
              :gov_cb_para_sw_idx_in_plant))

    avr_cb_sw = avr_gov_cb_sw[avr_cb_sw_Idx ]
    
    gov_cb_sw = avr_gov_cb_sw[gov_cb_sw_Idx]

    avr_kwd_para_wt_avr_cb_sw =
        (avr_para, avr_cb_sw )
    
    gov_kwd_para_wt_gov_cb_sw =
        (gov_para, gov_cb_sw)
    
    #----------------------------------------
    #----------------------------------------

    states_dims  = length.(comps_Idx)
    state_offset = create_offsets(states_dims)
    state_Idx    = create_idxs(state_offset, states_dims)
    
    res_gen = @view res[state_Idx[1]] 
    res_avr = @view res[state_Idx[2]]
    res_gov = @view res[state_Idx[3]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]
    du_gov = @view dx[state_Idx[3]]
    
    u_gen  = @view x[state_Idx[1]]
    u_avr  = @view x[state_Idx[2]]    
    u_gov  = @view x[state_Idx[3]]    
        
    #----------------------------------------
    #----------------------------------------
    
    ωref0, vref0, porder0, id, iq, vh =
        a_ωref0_vref0_porder0_id_iq_vh

    #----------------------------------------

    avr_vh_Vref0_para =
        [vh, vref0]
    
    avr_fun!(
        res_avr,
        du_avr,
        u_avr,
        avr_vh_Vref0_para,
        t;
        kwd_para =
            avr_kwd_para_wt_avr_cb_sw )
    
    avr_vf_tilade = avr_output_fun(
        u_avr,
        avr_vh_Vref0_para;
        kwd_para =
            avr_para )

    #----------------------------------------
    
    gen_ω = gen_output_fun(
        u_gen,
        a_ωref0_vref0_porder0_id_iq_vh;
        kwd_para =
            gen_para  )

    gov_gen_ω_ωref0_porder0_para =
        [gen_ω, ωref0, porder0, ]

    # dae_gov_sauer_func!
    
    gov_fun!(
        res_gov,
        du_gov,
        u_gov,
        gov_gen_ω_ωref0_porder0_para,
        t;
        kwd_para =
            gov_kwd_para_wt_gov_cb_sw )
    
    gov_τm_tilade = gov_output_fun(
        u_gov,
        gov_gen_ω_ωref0_porder0_para;
        kwd_para =
            gov_para )
    
    #----------------------------------------
    
    gen_id_iq_τmtilade_vftilade_para =
        [id, iq, gov_τm_tilade, avr_vf_tilade]

    #  dae_a_SM_gen_generic_model_func!
    
    gen_fun!(
        res_gen,
        du_gen,
        u_gen,
        gen_id_iq_τmtilade_vftilade_para,
        t;
        kwd_para =
            (gen_para, ωs) )

    #----------------------------------------

    return nothing
end

#----------------------------------------
# Dynamics: Synchronous condenser
#----------------------------------------

function ode_a_SC_plant_generic_model_func!(
    dx,
    x,
    a_ωref0_vref0_porder0_id_iq_vh,
    t;
    kwd_para =
        plant_kwd_para )

    (gen_para,
     avr_para,
     
     comps_Idx,
     
     # gov_avr_cb_sw,
     avr_gov_cb_sw,
     avr_gov_cb_sw_Idx,
          
     comps_dyn_fun,
     comps_output_fun,

     ωs ) =
         kwd_para
                 
    (gen_Idx,
     avr_Idx ) =
         comps_Idx
    
    (gen_fun!,
     avr_fun! ) =
         comps_dyn_fun

    
    (gen_output_fun,
     avr_output_fun) =
         comps_output_fun


    (avr_cb_sw_Idx,
     gov_cb_sw_Idx) =
         NamedTupleTools.select(
             avr_gov_cb_sw_Idx,
             (:avr_cb_para_sw_idx_in_plant,
              :gov_cb_para_sw_idx_in_plant))

    avr_cb_sw = avr_gov_cb_sw[avr_cb_sw_Idx ]
    
    gov_cb_sw = avr_gov_cb_sw[gov_cb_sw_Idx]

    avr_kwd_para_wt_avr_cb_sw =
        (avr_para, avr_cb_sw )
    
    # gov_kwd_para_wt_gov_cb_sw =
    #     (gov_para, gov_cb_sw)
    
    
    #----------------------------------------
    #----------------------------------------

    states_dims  = length.(comps_Idx)
    state_offset = create_offsets( states_dims)
    state_Idx    = create_idxs( state_offset, states_dims)
    
    u_gen  = @view x[state_Idx[1]]
    u_avr  = @view x[state_Idx[2]]    
    # u_gov  = @view x[state_Idx[3]]    
    
    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]
    # du_gov = @view dx[state_Idx[3]]
    
    #----------------------------------------
    #----------------------------------------

    ωref0, vref0, porder0, id, iq, vh =
        a_ωref0_vref0_porder0_id_iq_vh

    #----------------------------------------

    avr_vh_Vref0_para =
        [vh, vref0]
    
    avr_fun!(
        du_avr,
        u_avr,
        avr_vh_Vref0_para,
        t;
        kwd_para =
            avr_kwd_para_wt_avr_cb_sw )
    
    avr_vf_tilade = avr_output_fun(
        u_avr,
        avr_vh_Vref0_para;
        kwd_para =
            avr_para )
    
    #----------------------------------------
    
    gen_ω = gen_output_fun(
        u_gen,
        a_ωref0_vref0_porder0_id_iq_vh;
        kwd_para =
            gen_para  )
    
    #----------------------------------------
    
    gov_τm_tilade = porder0
    
    gen_id_iq_τmtilade_vftilade_para =
        [id, iq, gov_τm_tilade, avr_vf_tilade]

    #  ode_a_SM_gen_generic_model_func!
    
    gen_fun!(
        du_gen,
        u_gen,
        gen_id_iq_τmtilade_vftilade_para,
        t;
        kwd_para =
            (gen_para, ωs) )
    
    #----------------------------------------

    return nothing
end



function dae_a_SC_plant_generic_model_func!(
    res,
    dx, x,
    a_ωref0_vref0_porder0_id_iq_vh,
    t;
    kwd_para =
        plant_kwd_para )

    (gen_para,
     avr_para,
     
     comps_Idx,
     
     avr_gov_cb_sw,
     avr_gov_cb_sw_Idx,

     
     comps_dyn_fun,
     comps_output_fun,

     ωs ) =
         kwd_para
             
    
    (gen_Idx,
     avr_Idx ) =
         comps_Idx
    
    (gen_fun!,
     avr_fun! ) =
         comps_dyn_fun

    
    (gen_output_fun,
     avr_output_fun) =
         comps_output_fun

    (avr_cb_sw_Idx,
     gov_cb_sw_Idx) =
         NamedTupleTools.select(
             avr_gov_cb_sw_Idx,
             (:avr_cb_para_sw_idx_in_plant,
              :gov_cb_para_sw_idx_in_plant))

    avr_cb_sw = avr_gov_cb_sw[avr_cb_sw_Idx ]
    
    gov_cb_sw = avr_gov_cb_sw[gov_cb_sw_Idx]

    avr_kwd_para_wt_avr_cb_sw =
        (avr_para, avr_cb_sw )
    
    # gov_kwd_para_wt_gov_cb_sw =
    #     (gov_para, gov_cb_sw)
    
    #----------------------------------------
    #----------------------------------------

    states_dims  = length.(comps_Idx)
    
    state_offset = create_offsets(states_dims)
    
    state_Idx    = create_idxs(
        state_offset, states_dims)
    
    res_gen = @view res[state_Idx[1]] 
    res_avr = @view res[state_Idx[2]]
    # res_gov = @view res[state_Idx[3]]
    
    du_gen = @view dx[state_Idx[1]] 
    du_avr = @view dx[state_Idx[2]]
    # du_gov = @view dx[state_Idx[3]]
    
    u_gen  = @view x[state_Idx[1]]
    u_avr  = @view x[state_Idx[2]]    
    # u_gov  = @view x[state_Idx[3]]    
    
    #----------------------------------------
    #----------------------------------------
    
    ωref0, vref0, porder0, id, iq, vh =
        a_ωref0_vref0_porder0_id_iq_vh


    #----------------------------------------

    # avr_vh_Vref0_para =
    #     [vh, vref0]
    # # dae_avr_sauer_func!
    
    # avr_fun!(
    #     res[avr_Idx],
    #     dx[avr_Idx],
    #     x[avr_Idx],
    #     avr_vh_Vref0_para,
    #     t;
    #     kwd_para =
    #         avr_kwd_para_wt_avr_cb_sw )
    
    # avr_vf_tilade = avr_output_fun(
    #     x[avr_Idx],
    #     avr_vh_Vref0_para;
    #     kwd_para =
    #         avr_kwd_para )

    avr_vh_Vref0_para =
        [vh, vref0]
    # dae_avr_sauer_func!
    
    avr_fun!(
        res_avr,
        du_avr,
        u_avr,
        avr_vh_Vref0_para,
        t;
        kwd_para =
            avr_kwd_para_wt_avr_cb_sw )
    
    avr_vf_tilade = avr_output_fun(
        u_avr,
        avr_vh_Vref0_para;
        kwd_para =
            avr_para )
    
    #----------------------------------------
    
    gen_ω = gen_output_fun(
        u_gen,
        a_ωref0_vref0_porder0_id_iq_vh;
        kwd_para =
            gen_para  )

    #----------------------------------------

    # # gov_τm_tilade = 0

    # gov_τm_tilade = porder0
    
    # gen_id_iq_τmtilade_vftilade_para =
    #     [id, iq, gov_τm_tilade, avr_vf_tilade]

    # #  dae_a_SM_gen_generic_model_func!
    
    # gen_fun!(
    #     res[gen_Idx],
    #     dx[gen_Idx],
    #     x[gen_Idx],
    #     gen_id_iq_τmtilade_vftilade_para,
    #     t;
    #     kwd_para =
    #         (gen_kwd_para, ωs) )


    # gov_τm_tilade = 0

    gov_τm_tilade = porder0
    
    gen_id_iq_τmtilade_vftilade_para =
        [id, iq, gov_τm_tilade, avr_vf_tilade]

    #  dae_a_SM_gen_generic_model_func!
    
    gen_fun!(
        res_gen,
        du_gen,
        u_gen,
        gen_id_iq_τmtilade_vftilade_para,
        t;
        kwd_para =
            (gen_para, ωs) )
    
    #----------------------------------------

    return nothing
end

#----------------------------------------
# Plants
#----------------------------------------

function ode_gens_plants_generic_model_wt_cb_state_func!(
    dx,
    x,
    model_dynamics_para,
    t;
    kwd_para =
        plants_kwd_para )

    (ωref0_vref0_porder0_id_iq_vh,
     plants_cb_paras_switches) =
        model_dynamics_para

    (state_vars_idx,

     ωref0_vref0_porder0_id_iq_vh_Idx,

     generic_gens_para,
     generic_avrs_para,
     generic_govs_para,

     vec_comp_states_Idx,

     # govs_avrs_cb_sw,
     avrs_govs_cb_sw,
     avrs_govs_cb_sw_Idx,
    
     comps_dyn_funs,
     comps_output_funs,

     ωs) =
        NamedTupleTools.select(
            kwd_para,
            (:state_vars_idx,

             :ωref0_vref0_porder0_id_iq_vh_Idx,

             :generic_gens_para,
             :generic_avrs_para,
             :generic_govs_para,

             :vec_comp_states_Idx,

             # :govs_avrs_cb_sw,
             :avrs_govs_cb_sw,
             :avrs_govs_cb_sw_Idx,

             # :comps_dyn_funs,
             :ode_comps_dyn_funs,
             :comps_output_funs,

             :ωs))
    
    #----------------------------------------

    du_views = [view(dx, idx)
                for idx in state_vars_idx]
    
    u_views  = [view(x, idx)
                for idx in state_vars_idx]

    #----------------------------------------
    
    # (ωref0_Idx, 
    #  vref0_Idx,    
    #  porder0_Idx,
    #  id_Idx,    
    #  iq_Idx, 
    #  vh_Idx) = ωref0_vref0_porder0_id_iq_vh_Idx

    ωref0_Idx =
        ωref0_vref0_porder0_id_iq_vh_Idx.ωref0_Idx
    
    vref0_Idx =
        ωref0_vref0_porder0_id_iq_vh_Idx.vref0_Idx
    
    porder0_Idx =
        ωref0_vref0_porder0_id_iq_vh_Idx.porder0_Idx
    
    id_Idx =
        ωref0_vref0_porder0_id_iq_vh_Idx.id_Idx
    
    iq_Idx =
        ωref0_vref0_porder0_id_iq_vh_Idx.iq_Idx
    
    vh_Idx =
        ωref0_vref0_porder0_id_iq_vh_Idx.vh_Idx

    #----------------------------------------
    
    ωref0 =
        ωref0_vref0_porder0_id_iq_vh[
            ωref0_Idx]
    
    vref0 = ωref0_vref0_porder0_id_iq_vh[
        vref0_Idx]
    
    porder0 =
        ωref0_vref0_porder0_id_iq_vh[
            porder0_Idx]
    
    id =
        ωref0_vref0_porder0_id_iq_vh[
            id_Idx]
    
    iq =
        ωref0_vref0_porder0_id_iq_vh[
            iq_Idx]
    
    vh = ωref0_vref0_porder0_id_iq_vh[
        vh_Idx]

    
    #----------------------------------------
             
    gens_dyn_funs! =
        [an_item.gen_dyn_fun for an_item in
             comps_dyn_funs]

    avrs_dyn_funs! = 
        [an_item.avr_dyn_fun for an_item in
             comps_dyn_funs]

    govs_dyn_funs! = 
        [an_item.gov_dyn_fun for an_item in
             comps_dyn_funs]

    
    #----------------------------------------
             
    gens_output_funs =
        [an_item.gen_output_fun for an_item in
             comps_output_funs]

    avrs_output_funs = 
        [an_item.avr_output_fun for an_item in
             comps_output_funs]

    govs_output_funs = 
        [an_item.gov_output_fun for an_item in
             comps_output_funs]
    
    #----------------------------------------

    
    for (state_var_idx,
         du_view,
         u_view,

         a_ωref0,
         a_vref0,
         a_porder0,
         a_id,
         a_iq,
         a_vh,
             
         gen_para,
         avr_para,
         gov_para,

         comps_Idx,

         avr_gov_cb_sw,
         avr_gov_cb_sw_Idx,

         gen_dyn_fun!,
         avr_dyn_fun!,
         gov_dyn_fun!,
         
         gen_output_fun,
         avr_output_fun,
         gov_output_fun) in
        zip( state_vars_idx,
             du_views,
             u_views,

             ωref0,
             vref0,
             porder0,
             id,
             iq,
             vh,
             
             generic_gens_para,
             generic_avrs_para,
             generic_govs_para,

             vec_comp_states_Idx,

             # govs_avrs_cb_sw,
             avrs_govs_cb_sw,
             avrs_govs_cb_sw_Idx,

             gens_dyn_funs!,
             avrs_dyn_funs!,
             govs_dyn_funs!,
         
             gens_output_funs,
             avrs_output_funs,
             govs_output_funs)

            a_ωref0_vref0_porder0_id_iq_vh =
                [a_ωref0,
                 a_vref0,
                 a_porder0,
                 a_id,
                 a_iq,
                 a_vh]
        
        if Symbol(split(String(
            nameof(gov_dyn_fun!)),
                        "__")[2]) == :nothing
            
            comps_dyn_fun =
                (gen_dyn_fun!,
                 avr_dyn_fun! )
                    
            comps_output_fun =
                (gen_output_fun,
                 avr_output_fun )
            
            plant_kwd_para =
                (gen_para,
                 avr_para,

                 comps_Idx,

                 #gov_avr_cb_sw,
                 avr_gov_cb_sw,
                 avr_gov_cb_sw_Idx,

                 comps_dyn_fun,
                 comps_output_fun,
                 
                  ωs) 
                     
            ode_a_SC_plant_generic_model_func!(
                # dx[state_var_idx],
                # x[state_var_idx],
                du_view,
                u_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    plant_kwd_para )
            
        else

            comps_dyn_fun =
                (gen_dyn_fun!,
                 avr_dyn_fun!,
                 gov_dyn_fun! )
         
                    
            comps_output_fun =
                (gen_output_fun,
                 avr_output_fun,
                 gov_output_fun )
            
            plant_kwd_para =
                (gen_para,
                 avr_para,
                 gov_para,

                 comps_Idx,

                 # gov_avr_cb_sw,
                 avr_gov_cb_sw,
                 avr_gov_cb_sw_Idx,

                 comps_dyn_fun,
                 comps_output_fun,
                 
                  ωs) 
            
            ode_a_SM_plant_generic_model_func!(
                # dx[state_var_idx],
                # x[state_var_idx],
                du_view,
                u_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    plant_kwd_para )
                
        end        
    end
    
    #----------------------------------------

    return nothing
end



function ode_cb_test_gens_plants_generic_model_func!(
    dx, x,
    ωref0_vref0_porder0_id_iq_vh,
    # gens_plants_para,
    t;
    kwd_para =
        plants_kwd_para )

    # (ωref0_vref0_porder0_id_iq_vh,
    #  plants_cb_paras_switches) =
    #     gens_plants_para

    (state_vars_idx,

     ωref0_vref0_porder0_id_iq_vh_Idx,

     generic_gens_para,
     generic_avrs_para,
     generic_govs_para,

     vec_comp_states_Idx,

     # govs_avrs_cb_sw,
     avrs_govs_cb_sw,
     avrs_govs_cb_sw_Idx,
    
     comps_dyn_funs,
     comps_output_funs,

     ωs) =
        NamedTupleTools.select(
            kwd_para,
            (:state_vars_idx,

             :ωref0_vref0_porder0_id_iq_vh_Idx,

             :generic_gens_para,
             :generic_avrs_para,
             :generic_govs_para,

             :vec_comp_states_Idx,

             # :govs_avrs_cb_sw,
             :avrs_govs_cb_sw,
             :avrs_govs_cb_sw_Idx,

             # :comps_dyn_funs,
             :ode_comps_dyn_funs,
             :comps_output_funs,

             :ωs))
    
    #----------------------------------------

    du_views = [view(dx, idx)
                for idx in state_vars_idx]
    
    u_views  = [view(x, idx)
                for idx in state_vars_idx]

    #----------------------------------------
    
    # (ωref0_Idx, 
    #  vref0_Idx,    
    #  porder0_Idx,
    #  id_Idx,    
    #  iq_Idx, 
    #  vh_Idx) = ωref0_vref0_porder0_id_iq_vh_Idx

    ωref0_Idx =
        ωref0_vref0_porder0_id_iq_vh_Idx.ωref0_Idx
    
    vref0_Idx =
        ωref0_vref0_porder0_id_iq_vh_Idx.vref0_Idx
    
    porder0_Idx =
        ωref0_vref0_porder0_id_iq_vh_Idx.porder0_Idx
    
    id_Idx =
        ωref0_vref0_porder0_id_iq_vh_Idx.id_Idx
    
    iq_Idx =
        ωref0_vref0_porder0_id_iq_vh_Idx.iq_Idx
    
    vh_Idx =
        ωref0_vref0_porder0_id_iq_vh_Idx.vh_Idx

    #----------------------------------------
    
    ωref0 =
        ωref0_vref0_porder0_id_iq_vh[
            ωref0_Idx]
    
    vref0 = ωref0_vref0_porder0_id_iq_vh[
        vref0_Idx]
    
    porder0 =
        ωref0_vref0_porder0_id_iq_vh[
            porder0_Idx]
    
    id =
        ωref0_vref0_porder0_id_iq_vh[
            id_Idx]
    
    iq =
        ωref0_vref0_porder0_id_iq_vh[
            iq_Idx]
    
    vh = ωref0_vref0_porder0_id_iq_vh[
        vh_Idx]

    
    #----------------------------------------
             
    gens_dyn_funs! =
        [an_item.gen_dyn_fun for an_item in
             comps_dyn_funs]

    avrs_dyn_funs! = 
        [an_item.avr_dyn_fun for an_item in
             comps_dyn_funs]

    govs_dyn_funs! = 
        [an_item.gov_dyn_fun for an_item in
             comps_dyn_funs]

    
    #----------------------------------------
             
    gens_output_funs =
        [an_item.gen_output_fun for an_item in
             comps_output_funs]

    avrs_output_funs = 
        [an_item.avr_output_fun for an_item in
             comps_output_funs]

    govs_output_funs = 
        [an_item.gov_output_fun for an_item in
             comps_output_funs]
    
    #----------------------------------------

    
    for (state_var_idx,
         du_view,
         u_view,

         a_ωref0,
         a_vref0,
         a_porder0,
         a_id,
         a_iq,
         a_vh,
             
         gen_para,
         avr_para,
         gov_para,

         comps_Idx,

         avr_gov_cb_sw,
         avr_gov_cb_sw_Idx,

         gen_dyn_fun!,
         avr_dyn_fun!,
         gov_dyn_fun!,
         
         gen_output_fun,
         avr_output_fun,
         gov_output_fun) in
        zip( state_vars_idx,
             du_views,
             u_views,

             ωref0,
             vref0,
             porder0,
             id,
             iq,
             vh,
             
             generic_gens_para,
             generic_avrs_para,
             generic_govs_para,

             vec_comp_states_Idx,

             # govs_avrs_cb_sw,
             avrs_govs_cb_sw,
             avrs_govs_cb_sw_Idx,

             gens_dyn_funs!,
             avrs_dyn_funs!,
             govs_dyn_funs!,
         
             gens_output_funs,
             avrs_output_funs,
             govs_output_funs)

            a_ωref0_vref0_porder0_id_iq_vh =
                [a_ωref0,
                 a_vref0,
                 a_porder0,
                 a_id,
                 a_iq,
                 a_vh]
        
        if Symbol(split(String(
            nameof(gov_dyn_fun!)),
                        "__")[2]) == :nothing
            
            comps_dyn_fun =
                (gen_dyn_fun!,
                 avr_dyn_fun! )
                    
            comps_output_fun =
                (gen_output_fun,
                 avr_output_fun )
            
            plant_kwd_para =
                (gen_para,
                 avr_para,

                 comps_Idx,

                 #gov_avr_cb_sw,
                 avr_gov_cb_sw,
                 avr_gov_cb_sw_Idx,

                 comps_dyn_fun,
                 comps_output_fun,
                 
                  ωs) 
                     
            ode_a_SC_plant_generic_model_func!(
                # dx[state_var_idx],
                # x[state_var_idx],
                du_view,
                u_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    plant_kwd_para )
            
        else

            comps_dyn_fun =
                (gen_dyn_fun!,
                 avr_dyn_fun!,
                 gov_dyn_fun! )
         
                    
            comps_output_fun =
                (gen_output_fun,
                 avr_output_fun,
                 gov_output_fun )
            
            plant_kwd_para =
                (gen_para,
                 avr_para,
                 gov_para,

                 comps_Idx,

                 # gov_avr_cb_sw,
                 avr_gov_cb_sw,
                 avr_gov_cb_sw_Idx,

                 comps_dyn_fun,
                 comps_output_fun,
                 
                  ωs) 
            
            ode_a_SM_plant_generic_model_func!(
                # dx[state_var_idx],
                # x[state_var_idx],
                du_view,
                u_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    plant_kwd_para )
                
        end        
    end
    
    #----------------------------------------

    return nothing
end


function ode_gens_plants_generic_model_func!(
    dx, x,
    ωref0_vref0_porder0_id_iq_vh,
    t;
    kwd_para =
        plants_kwd_para )

    (state_vars_idx,

     ωref0_vref0_porder0_id_iq_vh_Idx,

     generic_gens_para,
     generic_avrs_para,
     generic_govs_para,

     vec_comp_states_Idx,

     avrs_govs_cb_sw,
     avrs_govs_cb_sw_Idx,
    
     comps_dyn_funs,
     comps_output_funs,

     ωs) =
        NamedTupleTools.select(
            kwd_para,
            (:state_vars_idx,

             :ωref0_vref0_porder0_id_iq_vh_Idx,

             :generic_gens_para,
             :generic_avrs_para,
             :generic_govs_para,

             :vec_comp_states_Idx,

             :avrs_govs_cb_sw,
             :avrs_govs_cb_sw_Idx,

             # :comps_dyn_funs,
             :ode_comps_dyn_funs,
             :comps_output_funs,

             :ωs))
    
    #----------------------------------------

    du_views = [view(dx, idx)
                for idx in state_vars_idx]
    
    u_views  = [view(x, idx)
                for idx in state_vars_idx]
    
    #----------------------------------------

    (ωref0_Idx, 
     vref0_Idx,    
     porder0_Idx,
     id_Idx,    
     iq_Idx, 
     vh_Idx) =
         NamedTupleTools.select(
             ωref0_vref0_porder0_id_iq_vh_Idx,
             (:ωref0_Idx, 
              :vref0_Idx,    
              :porder0_Idx,
              :id_Idx,    
              :iq_Idx, 
              :vh_Idx))

    #----------------------------------------
    
    ωref0 =
        ωref0_vref0_porder0_id_iq_vh[
            ωref0_Idx]
    
    vref0 = ωref0_vref0_porder0_id_iq_vh[
        vref0_Idx]
    
    porder0 =
        ωref0_vref0_porder0_id_iq_vh[
            porder0_Idx]
    
    id =
        ωref0_vref0_porder0_id_iq_vh[
            id_Idx]
    
    iq =
        ωref0_vref0_porder0_id_iq_vh[
            iq_Idx]
    
    vh = ωref0_vref0_porder0_id_iq_vh[
        vh_Idx]
    
    #----------------------------------------
             
    gens_dyn_funs! =
        [an_item.gen_dyn_fun for an_item in
             comps_dyn_funs]

    avrs_dyn_funs! = 
        [an_item.avr_dyn_fun for an_item in
             comps_dyn_funs]

    govs_dyn_funs! = 
        [an_item.gov_dyn_fun for an_item in
             comps_dyn_funs]

    
    #----------------------------------------
             
    gens_output_funs =
        [an_item.gen_output_fun for an_item in
             comps_output_funs]

    avrs_output_funs = 
        [an_item.avr_output_fun for an_item in
             comps_output_funs]

    govs_output_funs = 
        [an_item.gov_output_fun for an_item in
             comps_output_funs]
    
    #----------------------------------------

    
    for (state_var_idx,
         
         du_view,
         u_view,

         a_ωref0,
         a_vref0,
         a_porder0,
         a_id,
         a_iq,
         a_vh,
             
         gen_para,
         avr_para,
         gov_para,

         comps_Idx,

         #gov_avr_cb_sw,
         avr_gov_cb_sw,
         avr_gov_cb_sw_Idx,

         gen_dyn_fun!,
         avr_dyn_fun!,
         gov_dyn_fun!,
         
         gen_output_fun,
         avr_output_fun,
         gov_output_fun) in
        zip( state_vars_idx,

             du_views,
             u_views,
             
             ωref0,
             vref0,
             porder0,
             id,
             iq,
             vh,
             
             generic_gens_para,
             generic_avrs_para,
             generic_govs_para,

             vec_comp_states_Idx,

             # govs_avrs_cb_sw,
             avrs_govs_cb_sw,
             avrs_govs_cb_sw_Idx,

             gens_dyn_funs!,
             avrs_dyn_funs!,
             govs_dyn_funs!,
         
             gens_output_funs,
             avrs_output_funs,
             govs_output_funs)

            a_ωref0_vref0_porder0_id_iq_vh =
                [a_ωref0,
                 a_vref0,
                 a_porder0,
                 a_id,
                 a_iq,
                 a_vh]
        
        if Symbol(split(String(
            nameof(
                gov_dyn_fun!)),"__")[2]) == :nothing
            
            comps_dyn_fun =
                (gen_dyn_fun!,
                 avr_dyn_fun! )
                    
            comps_output_fun =
                (gen_output_fun,
                 avr_output_fun )
            
            plant_kwd_para =
                (gen_para,
                 avr_para,

                 comps_Idx,

                 #gov_avr_cb_sw,
                 avr_gov_cb_sw,
                 avr_gov_cb_sw_Idx,

                 comps_dyn_fun,
                 comps_output_fun,
                 
                  ωs) 
                     
            ode_a_SC_plant_generic_model_func!(
                # dx[state_var_idx],
                # x[state_var_idx],
                du_view,
                u_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    plant_kwd_para )
            
        else

            comps_dyn_fun =
                (gen_dyn_fun!,
                 avr_dyn_fun!,
                 gov_dyn_fun! )
         
                    
            comps_output_fun =
                (gen_output_fun,
                 avr_output_fun,
                 gov_output_fun )
            
            plant_kwd_para =
                (gen_para,
                 avr_para,
                 gov_para,

                 comps_Idx,

                 avr_gov_cb_sw,
                 avr_gov_cb_sw_Idx,

                 comps_dyn_fun,
                 comps_output_fun,
                 
                  ωs) 
            
            ode_a_SM_plant_generic_model_func!(
                # dx[state_var_idx],
                # x[state_var_idx],
                du_view,
                u_view,                             
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    plant_kwd_para )
                
        end        
    end
    
    #----------------------------------------

    return nothing
end



function dae_gens_plants_generic_model_func!(
    res,
    dx, x,
    ωref0_vref0_porder0_id_iq_vh,
    t;
    kwd_para =
        plants_kwd_para )

    (state_vars_idx,

     ωref0_vref0_porder0_id_iq_vh_Idx,

     generic_gens_para,
     generic_avrs_para,
     generic_govs_para,

     vec_comp_states_Idx,

     # govs_avrs_cb_sw,
     avrs_govs_cb_sw,
     avrs_govs_cb_sw_Idx,
    
     comps_dyn_funs,
     comps_output_funs,

     ωs) =
        NamedTupleTools.select(
            kwd_para,
            (:state_vars_idx,

             :ωref0_vref0_porder0_id_iq_vh_Idx,

             :generic_gens_para,
             :generic_avrs_para,
             :generic_govs_para,

             :vec_comp_states_Idx,

             # :govs_avrs_cb_sw,
             :avrs_govs_cb_sw,
             :avrs_govs_cb_sw_Idx,

             # :comps_dyn_funs,
             :dae_comps_dyn_funs,
             :comps_output_funs,

             :ωs))

    
    #----------------------------------------

    res_views = [view(res, idx)
                for idx in state_vars_idx]

    du_views = [view(dx, idx)
                for idx in state_vars_idx]
    
    u_views  = [view(x, idx)
                for idx in state_vars_idx]

    #----------------------------------------

    (ωref0_Idx,
     vref0_Idx,
     porder0_Idx,
     id_Idx,
     iq_Idx,
     vh_Idx) =
         NamedTupleTools.select(
             ωref0_vref0_porder0_id_iq_vh_Idx,
             (:ωref0_Idx,
              :vref0_Idx,
              :porder0_Idx,
              :id_Idx,
              :iq_Idx,
              :vh_Idx))
    
    #----------------------------------------
    
    ωref0 =
        ωref0_vref0_porder0_id_iq_vh[
            ωref0_Idx]
    
    vref0 = ωref0_vref0_porder0_id_iq_vh[
        vref0_Idx]
    
    porder0 =
        ωref0_vref0_porder0_id_iq_vh[
            porder0_Idx]
    
    id =
        ωref0_vref0_porder0_id_iq_vh[
            id_Idx]
    
    iq =
        ωref0_vref0_porder0_id_iq_vh[
            iq_Idx]
    
    vh = ωref0_vref0_porder0_id_iq_vh[
        vh_Idx]

    
    #----------------------------------------
             
    gens_dyn_funs! =
        [an_item.gen_dyn_fun for an_item in
             comps_dyn_funs]

    avrs_dyn_funs! = 
        [an_item.avr_dyn_fun for an_item in
             comps_dyn_funs]

    govs_dyn_funs! = 
        [an_item.gov_dyn_fun for an_item in
             comps_dyn_funs]

    
    #----------------------------------------
             
    gens_output_funs =
        [an_item.gen_output_fun for an_item in
             comps_output_funs]

    avrs_output_funs = 
        [an_item.avr_output_fun for an_item in
             comps_output_funs]

    govs_output_funs = 
        [an_item.gov_output_fun for an_item in
             comps_output_funs]
    
    #----------------------------------------

    
    for (state_var_idx,

         res_view,
         du_view,
         u_view,
         
         a_ωref0,
         a_vref0,
         a_porder0,
         a_id,
         a_iq,
         a_vh,
             
         gen_para,
         avr_para,
         gov_para,

         comps_Idx,

         #gov_avr_cb_sw,
         avr_gov_cb_sw,
         avr_gov_cb_sw_Idx,


         gen_dyn_fun!,
         avr_dyn_fun!,
         gov_dyn_fun!,
         
         gen_output_fun,
         avr_output_fun,
         gov_output_fun) in
        zip( state_vars_idx,

             res_views,
             du_views,
             u_views,

             ωref0,
             vref0,
             porder0,
             id,
             iq,
             vh,
             
             generic_gens_para,
             generic_avrs_para,
             generic_govs_para,

             vec_comp_states_Idx,

             # govs_avrs_cb_sw,
             avrs_govs_cb_sw,
             avrs_govs_cb_sw_Idx,

             gens_dyn_funs!,
             avrs_dyn_funs!,
             govs_dyn_funs!,
         
             gens_output_funs,
             avrs_output_funs,
             govs_output_funs)

            a_ωref0_vref0_porder0_id_iq_vh =
                [a_ωref0,
                 a_vref0,
                 a_porder0,
                 a_id,
                 a_iq,
                 a_vh]
        
        if Symbol(split(String(
            nameof(
                gov_dyn_fun!)),"__")[2]) == :nothing
            
            comps_dyn_fun =
                (gen_dyn_fun!,
                 avr_dyn_fun! )
                    
            comps_output_fun =
                (gen_output_fun,
                 avr_output_fun )
            
            plant_kwd_para =
                (gen_para,
                 avr_para,

                 comps_Idx,

                 # gov_avr_cb_sw,
                 avr_gov_cb_sw,
                 avr_gov_cb_sw_Idx,

                 comps_dyn_fun,
                 comps_output_fun,
                 
                  ωs) 
                     
            dae_a_SC_plant_generic_model_func!(
                # res[state_var_idx],
                # dx[state_var_idx],
                # x[state_var_idx],
                res_view,
                du_view,
                u_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    plant_kwd_para )
            
        else

            comps_dyn_fun =
                (gen_dyn_fun!,
                 avr_dyn_fun!,
                 gov_dyn_fun! )
         
                    
            comps_output_fun =
                (gen_output_fun,
                 avr_output_fun,
                 gov_output_fun )
            
            plant_kwd_para =
                (gen_para,
                 avr_para,
                 gov_para,

                 comps_Idx,

                 # gov_avr_cb_sw,
                 avr_gov_cb_sw,
                 avr_gov_cb_sw_Idx,
                 

                 comps_dyn_fun,
                 comps_output_fun,
                 
                  ωs) 
            
            dae_a_SM_plant_generic_model_func!(
                # res[state_var_idx],
                # dx[state_var_idx],
                # x[state_var_idx],
                res_view,
                du_view,
                u_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    plant_kwd_para )
                
        end        
    end
    
    #----------------------------------------

    return nothing
end


#-----------------------------------------------------
# generic system model dynamics
#-----------------------------------------------------

function ode_generic_system_model_dynamics!(
    dx,
    x,
    model_dynamics_para,
    t;
    kwd_para =
        generic_model_dynamics_kwd_para  )

    #----------------------------------------

    (;generic_model_dynamics_para,
     plants_cb_paras_switches) =
         model_dynamics_para
    
    #----------------------------------------

    (;ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx, 

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     # pf_generic_gens_para, #

     plants_kwd_para,

     algebraic_generic_model_sol_kwd_para ) =
         NamedTupleTools.select(
             kwd_para,
             (:ωs,
              :loc_load_exist,
              :state_vars_idx,

              :id_iq_pg_vh_Idx,

             :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              # :pf_generic_gens_para,
              
              :plants_kwd_para,

              :algebraic_generic_model_sol_kwd_para))


    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))

    #----------------------------------------

    (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    # -------------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))  

    #----------------------------------------

    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------
    
    du_state_vars_view =
        @view dx[state_var_Idx_in_state]
    
    du_algebraic_vars_view =
        @view dx[algebraic_var_Idx_in_state]
    
    u_state_vars_view =
        @view x[state_var_Idx_in_state]
    
    u_algebraic_vars_view =
        @view x[algebraic_var_Idx_in_state]

    #----------------------------------------
    
    vh  = u_algebraic_vars_view[
        dyn_pf_vh_Idxs]
    
    θh  = u_algebraic_vars_view[
        dyn_pf_θh_Idxs]
    
    gens_i_d = u_algebraic_vars_view[
        dyn_pf_id_Idxs]
    
    gens_i_q = u_algebraic_vars_view[
        dyn_pf_iq_Idxs]
         
    #----------------------------------------
    
    gens_δ = u_state_vars_view[
        δ_idx_in_state]
        
    gens_ed_dash = u_state_vars_view[
        ed_dash_idx_in_state]
    
    gens_eq_dash = u_state_vars_view[
        eq_dash_idx_in_state]
    
    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        vh[ non_gens_nodes_idx ]

    non_gens_nodes_θh =
        θh[ non_gens_nodes_idx ]
            
    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        generic_model_dynamics_para[
            dyn_ω_ref_Idx]
    
    v_ref = generic_model_dynamics_para[
        dyn_v_ref_Idx]
    
    p_order = generic_model_dynamics_para[
        dyn_p_order_Idx]
    
    P_non_gens = generic_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = generic_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            generic_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    
    #----------------------------------------
    # ode parameters and kwd para
    #----------------------------------------
    
    gens_ωref0_vref0_porder0_id_iq_vh_para = [
        [a_ωref0, a_vref, a_porder0, a_id, a_iq, a_vh ]
        for (a_ωref0,
             a_vref,
             a_porder0,
             a_id,
             a_iq,
             a_vh) in
            zip( ω_ref,
                 v_ref,
                 p_order,
                 gens_i_d,
                 gens_i_q,
                 gens_vh) ]
    
    #----------------------------------------        
    #----------------------------------------
    # ode 
    #----------------------------------------
    #----------------------------------------

    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh]
    
    ode_gens_plants_generic_model_func!(
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    algebraic_generic_pf_ΔPQ_mismatch_sol(
        du_algebraic_vars_view,
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
        kwd_para =
            algebraic_generic_model_sol_kwd_para )
            
    # ------------------------------------
    
    return nothing

end



function ode_generic_system_model_by_funcs_dynamics!(
    dx,
    x,
    model_dynamics_para,
    t;
    kwd_para =
        generic_model_dynamics_kwd_para  )

    #----------------------------------------

    (;generic_model_dynamics_para,
      plants_cb_paras_switches ) =
         model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     ode_plants_kwd_para,

     algebraic_generic_model_sol_kwd_para ) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              :id_iq_pg_vh_Idx,

             :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_sol_kwd_para))


    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))

    #----------------------------------------

    (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))
   
    idx_range =
        first(dyn_pf_vh_Idxs):last(dyn_pf_iq_Idxs)
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    # -------------------------------------

    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------
    
    du_algebraic_vars_view =
        @view dx[algebraic_var_Idx_in_state]
    
    u_algebraic_vars_view =
        @view x[algebraic_var_Idx_in_state]

    #----------------------------------------
    
    du_state_vars_view =
        @view dx[state_var_Idx_in_state]
    
    u_state_vars_view =
        @view x[state_var_Idx_in_state]
    #----------------------------------------
    
    # gens_δ = u_state_vars_view[
    #     δ_idx_in_state]
        
    # gens_ed_dash = u_state_vars_view[
    #     ed_dash_idx_in_state]
    
    # gens_eq_dash = u_state_vars_view[
    #     eq_dash_idx_in_state]

    
    gens_δ = x[
        δ_idx_in_state]
        
    gens_ed_dash = x[
        ed_dash_idx_in_state]
    
    gens_eq_dash = x[
        eq_dash_idx_in_state]
    
    #----------------------------------------
    
    # vh = @view x[vh_Idx_in_state]
    
    # θh = @view x[θh_Idx_in_state]
    
    # gens_i_d = @view x[id_Idx_in_state]
    
    # gens_i_q = @view x[ iq_Idx_in_state]

    
    vh = x[vh_Idx_in_state]
    
    θh = x[θh_Idx_in_state]
    
    gens_i_d = x[id_Idx_in_state]
    
    gens_i_q = x[ iq_Idx_in_state]
    
    #----------------------------------------
        
    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    # non_gens_nodes_vh =
    #     vh[ non_gens_nodes_idx ]

    # non_gens_nodes_θh =
    #     θh[ non_gens_nodes_idx ]
         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        generic_model_dynamics_para[
            dyn_ω_ref_Idx]
    
    v_ref = generic_model_dynamics_para[
        dyn_v_ref_Idx]
    
    p_order = generic_model_dynamics_para[
        dyn_p_order_Idx]
    
    P_non_gens = generic_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = generic_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            generic_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]
    
    #----------------------------------------
    #----------------------------------------
    # using per plant dynamics 
    #----------------------------------------
    #----------------------------------------

    (per_plant_gov_dyn_fun_sym_name,
     per_plant_generic_model_para,
     per_plant_generic_model_kwd_para) =
        get_generic_namedtuple_per_plant_para_wt_kwd_para(
             ωref0_vref0_porder0_id_iq_vh;
             kwd_para = ode_plants_kwd_para )

    for (
         a_du_state_vars_view,
         a_u_state_vars_view,
         
         a_gov_dyn_fun_sym_name,
         a_ωref0_vref0_porder0_id_iq_vh,
         a_per_kwd_para) in zip(
             
             du_states_vars_views,
             u_states_vars_views,
             
             per_plant_gov_dyn_fun_sym_name,
             per_plant_generic_model_para,
             per_plant_generic_model_kwd_para )
        
        if a_gov_dyn_fun_sym_name == :nothing
            
            ode_a_SC_plant_generic_model_func!(
                a_du_state_vars_view,
                a_u_state_vars_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    a_per_kwd_para )
        else

            ode_a_SM_plant_generic_model_func!(
                a_du_state_vars_view,
                a_u_state_vars_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    a_per_kwd_para )
            
        end
        

    end

    
    # #----------------------------------------
    # #----------------------------------------
    # # using aggregate plants dynamics 
    # #----------------------------------------
    # #----------------------------------------

    # ode_gens_plants_generic_model_func!(
    #     du_state_vars_view,
    #     u_state_vars_view,
    #     ωref0_vref0_porder0_id_iq_vh,
    #     t;
    #     kwd_para =
    #         ode_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    # pf_algebraic_generic_model_sol(
    #     dx[algebraic_var_Idx_in_state],
    #     x[algebraic_var_Idx_in_state],
    #     δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    #     kwd_para =
    #         algebraic_generic_model_sol_kwd_para )


    vh_θh_id_iq =
        [vh;
         θh;
         gens_i_d;
         gens_i_q]

    # ode_pf_algebraic_generic_model_sol_ext
    
    pf_sol = algebraic_generic_pf_ΔPQ_mismatch_sol(
        vh_θh_id_iq,
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
        kwd_para =
            algebraic_generic_model_sol_kwd_para )
    
        
    du_algebraic_vars_view .=
        pf_sol[ idx_range ] - u_algebraic_vars_view
        
    
    # ------------------------------------
    
    return nothing

end



function ode_generic_system_dynamics_by_ode_pf_funcs!(
    dx,
    x,
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
    t;
    kwd_para =
        generic_model_dynamics_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     ode_plants_kwd_para,

     algebraic_generic_model_sol_kwd_para ) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              :id_iq_pg_vh_Idx,

             :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_sol_kwd_para))


    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))

    #----------------------------------------

    (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))
   
    idx_range =
        first(dyn_pf_vh_Idxs):last(dyn_pf_iq_Idxs)
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    # -------------------------------------

    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------
    
    du_algebraic_vars_view =
        @view dx[algebraic_var_Idx_in_state]
    
    u_algebraic_vars_view =
        @view x[algebraic_var_Idx_in_state]

    #----------------------------------------
    
    du_state_vars_view =
        @view dx[state_var_Idx_in_state]
    
    u_state_vars_view =
        @view x[state_var_Idx_in_state]
    #----------------------------------------
    
    # gens_δ = u_state_vars_view[
    #     δ_idx_in_state]
        
    # gens_ed_dash = u_state_vars_view[
    #     ed_dash_idx_in_state]
    
    # gens_eq_dash = u_state_vars_view[
    #     eq_dash_idx_in_state]

    
    gens_δ = x[
        δ_idx_in_state]
        
    gens_ed_dash = x[
        ed_dash_idx_in_state]
    
    gens_eq_dash = x[
        eq_dash_idx_in_state]
    
    #----------------------------------------
    
    # vh = @view x[vh_Idx_in_state]
    
    # θh = @view x[θh_Idx_in_state]
    
    # gens_i_d = @view x[id_Idx_in_state]
    
    # gens_i_q = @view x[ iq_Idx_in_state]

    
    vh = x[vh_Idx_in_state]
    
    θh = x[θh_Idx_in_state]
    
    gens_i_d = x[id_Idx_in_state]
    
    gens_i_q = x[ iq_Idx_in_state]
    
    #----------------------------------------
        
    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    # non_gens_nodes_vh =
    #     vh[ non_gens_nodes_idx ]

    # non_gens_nodes_θh =
    #     θh[ non_gens_nodes_idx ]
         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_v_ref_Idx]
    
    p_order = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]
    
    P_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Png_Idx]
    
    Q_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]
    
    #----------------------------------------
    #----------------------------------------
    # using plant dynamics 
    #----------------------------------------
    #----------------------------------------
    
    ode_gens_plants_generic_model_func!(
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            ode_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    # pf_algebraic_generic_model_sol(
    #     dx[algebraic_var_Idx_in_state],
    #     x[algebraic_var_Idx_in_state],
    #     δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    #     kwd_para =
    #         algebraic_generic_model_sol_kwd_para )


    vh_θh_id_iq =
        [vh;
         θh;
         gens_i_d;
         gens_i_q]

    # ode_pf_algebraic_generic_model_sol_ext
    
    pf_sol = algebraic_generic_pf_ΔPQ_mismatch_sol(
        vh_θh_id_iq,
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
        kwd_para =
            algebraic_generic_model_sol_kwd_para )
    
        
    du_algebraic_vars_view .=
        pf_sol[ idx_range ] - u_algebraic_vars_view
        
    
    # ------------------------------------
    
    return nothing

end



function dae_generic_system_model_dynamics!(
    res,
    dx,
    x,
    model_dynamics_para,
    t;
    kwd_para =
        generic_model_dynamics_kwd_para  )

    #----------------------------------------

    (;generic_model_dynamics_para,
     plants_cb_paras_switches ) =
         model_dynamics_para
    
    #----------------------------------------

    (;ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx, #

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_fun_kwd_n2s_idxs, #
     dyn_pf_fun_kwd_net_idxs, #

     pf_generic_gens_para, #

     Ynet_wt_nodes_idx_wt_adjacent_nodes,
     dae_plants_kwd_para ) =
         NamedTupleTools.select(
             kwd_para,
             (:ωs,
              :loc_load_exist,
              :state_vars_idx,

              :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_generic_gens_para,
              
              :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              :dae_plants_kwd_para ))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))

    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    # -------------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))  

    #----------------------------------------
    
    (ra,
     X_d,
     X_q,     
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
             pf_generic_gens_para,
             (:ra,
              :X_d,
              :X_q,     
              :X_d_dash,
              :X_q_dash ))
    
    #----------------------------------------

   (Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        NamedTupleTools.select(
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            (:Ynet,
             :nodes_idx_with_adjacent_nodes_idx)) 
    
    #----------------------------------------
    #----------------------------------------

     ode_res_views = @view res[state_var_Idx_in_state]

     ode_dx_views  = @view  dx[state_var_Idx_in_state]

     ode_x_views   =  @view x[state_var_Idx_in_state]

    #----------------------------------------
    #----------------------------------------

    state_vars =  @view x[state_var_Idx_in_state]
    
    vh  =  @view x[vh_Idx_in_state]
    
    θh  =  @view x[θh_Idx_in_state]
    
    gens_i_d =  @view x[id_Idx_in_state]
    
    gens_i_q =  @view x[iq_Idx_in_state]
    
    #----------------------------------------
    
    gens_δ = state_vars[δ_idx_in_state]
        
    gens_ed_dash = state_vars[ed_dash_idx_in_state]
    
    gens_eq_dash = state_vars[eq_dash_idx_in_state]
        
    #----------------------------------------
    
    gens_vh = vh[ gens_nodes_idx ]

    gens_θh = θh[ gens_nodes_idx ]

    non_gens_nodes_vh = vh[ non_gens_nodes_idx ]

    non_gens_nodes_θh = θh[ non_gens_nodes_idx ]
            
    #----------------------------------------
    
    # vh = [
    #     idx ∈ gens_nodes_idx ?
    #         gens_vh[
    #             n2s_gens_idx[ idx ] ]  :
    #                 non_gens_nodes_vh[
    #                     n2s_non_gens_idx[ idx ] ]
    #            for idx in all_nodes_idx ]

    # θh = [
    #     idx ∈ gens_nodes_idx ?
    #         gens_θh[
    #             n2s_gens_idx[ idx] ] :
    #                 non_gens_nodes_θh[
    #                     n2s_non_gens_idx[ idx ] ]
    #     for idx in all_nodes_idx ]
    
    #----------------------------------------
    
    ω_ref =
        generic_model_dynamics_para[
            dyn_ω_ref_Idx]
    
    v_ref = generic_model_dynamics_para[
        dyn_v_ref_Idx]
    
    p_order = generic_model_dynamics_para[
        dyn_p_order_Idx]
    
    P_non_gens = generic_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = generic_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            generic_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    
    #----------------------------------------
    # ode parameters and kwd para
    #----------------------------------------
    
    gens_ωref0_vref0_porder0_id_iq_vh_para = [
        [a_ωref0, a_vref, a_porder0, a_id, a_iq, a_vh ]
        for (a_ωref0, a_vref,
             a_porder0,
             a_id, a_iq,
             a_vh) in
            zip( ω_ref, v_ref,
                 p_order,
                 gens_i_d, gens_i_q,
                 gens_vh) ]
    
    #----------------------------------------        
    #----------------------------------------
    # ode 
    #----------------------------------------
    #----------------------------------------

    # dae_gens_plants_generic_model_func!(
    #     res[state_var_Idx_in_state],
    #     dx[state_var_Idx_in_state],
    #     x[state_var_Idx_in_state],
    #     ωref0_vref0_porder0_id_iq_vh,
    #     t;
    #     kwd_para =
    #         plants_kwd_para )

    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]
    
    dae_gens_plants_generic_model_func!(
        ode_res_views,
        ode_dx_views,
        ode_x_views,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            dae_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    # -----------------------------------
    # active power mismatch
    # -----------------------------------
    
    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
           vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """
    
    res[vh_Idx_in_state] .= [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * (gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -          
          sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                              nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])]))  -
                                   P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]]) :
        (vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])]) ))
                        for nth_idx in all_nodes_idx ] 

    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------
    
    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    res[θh_Idx_in_state]  .= [ nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          ( vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -           
           sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]]) :
           (vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])]))) 
                              for nth_idx in all_nodes_idx ]

    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """
    
    res[id_Idx_in_state] .= [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin( gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """
    
    res[iq_Idx_in_state] .= [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[ n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  
    
    # ------------------------------------
    
    return nothing

end


function dae_generic_system_model_by_funcs_dynamics!(
    res,
    dx,
    x,
    model_dynamics_para,
    t;
    kwd_para =
        generic_model_dynamics_kwd_para  )

    #----------------------------------------

    (;generic_model_dynamics_para,
      plants_cb_paras_switches ) =
         model_dynamics_para
    
    #----------------------------------------

    (;
     gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx, 

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     # dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 
     # pf_generic_gens_para, 
     # Ynet_wt_nodes_idx_wt_adjacent_nodes,
     
     dae_plants_kwd_para,

     algebraic_generic_model_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              # :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              # :pf_generic_gens_para,              
              # :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              
              :dae_plants_kwd_para,

              :algebraic_generic_model_kwd_para))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))

    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))
    
    #----------------------------------------
    
   (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))

    #----------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx))
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------

    res_states_vars_views =
        [ view(res, idx) for idx in state_vars_idx]
    
    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------
    
    res_state_vars_view =
        @view res[state_var_Idx_in_state]
    
    du_state_vars_view =
        @view dx[state_var_Idx_in_state]
    
    u_state_vars_view =
        @view x[state_var_Idx_in_state]

    #----------------------------------------
    
    gens_δ = u_state_vars_view[
        δ_idx_in_state]
        
    gens_ed_dash = u_state_vars_view[
        ed_dash_idx_in_state]
    
    gens_eq_dash = u_state_vars_view[
        eq_dash_idx_in_state]    
    
    #----------------------------------------
    
    res_algebraic_vars_view =
        @view res[algebraic_var_Idx_in_state]

    # note du for algebraic_vars = 0
    
    du_algebraic_vars_view =
        @view dx[algebraic_var_Idx_in_state]
    
    u_algebraic_vars_view =
        @view x[algebraic_var_Idx_in_state]

    #----------------------------------------

    vh  =  @view x[ vh_Idx_in_state]
    
    θh  =  @view x[ θh_Idx_in_state]
    
    gens_i_d =  @view x[id_Idx_in_state]
    
    gens_i_q =  @view x[iq_Idx_in_state]


    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        vh[
            non_gens_nodes_idx ]

    non_gens_nodes_θh =
        θh[
            non_gens_nodes_idx ]
    
    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        generic_model_dynamics_para[
            dyn_ω_ref_Idx]
    
    v_ref = generic_model_dynamics_para[
        dyn_v_ref_Idx]
    
    p_order = generic_model_dynamics_para[
        dyn_p_order_Idx]
    
    P_non_gens = generic_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = generic_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            generic_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end


    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]
    
    #----------------------------------------
    #----------------------------------------
    # using per plant dynamics 
    #----------------------------------------
    #----------------------------------------

    (per_plant_gov_dyn_fun_sym_name,
     per_plant_generic_model_para,
     per_plant_generic_model_kwd_para) =
         get_generic_namedtuple_per_plant_para_wt_kwd_para(
             ωref0_vref0_porder0_id_iq_vh;
             kwd_para = dae_plants_kwd_para )

    for (a_res_state_vars_view,
         a_du_state_vars_view,
         a_u_state_vars_view,
         
         a_gov_dyn_fun_sym_name,
         a_ωref0_vref0_porder0_id_iq_vh,
         a_per_kwd_para) in zip(
             res_states_vars_views,
             du_states_vars_views,
             u_states_vars_views,
             
             per_plant_gov_dyn_fun_sym_name,
             per_plant_generic_model_para,
             per_plant_generic_model_kwd_para )
        
        if a_gov_dyn_fun_sym_name == :nothing
            
            dae_a_SC_plant_generic_model_func!(
                a_res_state_vars_view,
                a_du_state_vars_view,
                a_u_state_vars_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    a_per_kwd_para )
        else

            dae_a_SM_plant_generic_model_func!(
                a_res_state_vars_view,
                a_du_state_vars_view,
                a_u_state_vars_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    a_per_kwd_para )
            
        end
        

    end

    
    # #----------------------------------------
    # #----------------------------------------
    # # using aggregate plants dynamics 
    # #----------------------------------------
    # #----------------------------------------

    # dae_gens_plants_generic_model_func!(
    #     res_state_vars_view,
    #     du_state_vars_view,
    #     u_state_vars_view,
    #     ωref0_vref0_porder0_id_iq_vh,
    #     t;
    #     kwd_para =
    #         dae_plants_kwd_para )
    
    #----------------------------------------
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    # dae_algebraic_generic_model_func!
    
    algebraic_generic_pf_ΔPQ_mismatch!(
        res_algebraic_vars_view,
        u_algebraic_vars_view,
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
        kwd_para =
            algebraic_generic_model_kwd_para  )
    
    return nothing

end


function dae_generic_system_dynamics_by_dae_pf_funcs!(
    res,
    dx,
    x,
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
    t;
    kwd_para =
        generic_model_dynamics_kwd_para  )

    #----------------------------------------

    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;
     gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx, 

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     # dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 
     # pf_generic_gens_para, 
     # Ynet_wt_nodes_idx_wt_adjacent_nodes,
     
     dae_plants_kwd_para,

     algebraic_generic_model_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              # :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              # :pf_generic_gens_para,              
              # :Ynet_wt_nodes_idx_wt_adjacent_nodes,
              
              :dae_plants_kwd_para,

              :algebraic_generic_model_kwd_para ))

    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))

    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))
    
    #----------------------------------------
    
   (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))

    #----------------------------------------

    (gens_nodes_idx,
     non_gens_nodes_idx,
     gens_nodes_with_loc_loads_idx,
     all_nodes_idx) =
         NamedTupleTools.select(
             dyn_pf_fun_kwd_net_idxs,
             (:gens_nodes_idx,
              :non_gens_nodes_idx,
              :gens_nodes_with_loc_loads_idx,
              :all_nodes_idx))
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    #----------------------------------------

    res_states_vars_views =
        [ view(res, idx) for idx in state_vars_idx]
    
    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------
    
    res_state_vars_view =
        @view res[state_var_Idx_in_state]
    
    du_state_vars_view =
        @view dx[state_var_Idx_in_state]
    
    u_state_vars_view =
        @view x[state_var_Idx_in_state]

    #----------------------------------------
    
    gens_δ = u_state_vars_view[
        δ_idx_in_state]
        
    gens_ed_dash = u_state_vars_view[
        ed_dash_idx_in_state]
    
    gens_eq_dash = u_state_vars_view[
        eq_dash_idx_in_state]    
    
    #----------------------------------------
    
    res_algebraic_vars_view =
        @view res[algebraic_var_Idx_in_state]

    # note du for algebraic_vars = 0
    
    du_algebraic_vars_view =
        @view dx[algebraic_var_Idx_in_state]
    
    u_algebraic_vars_view =
        @view x[algebraic_var_Idx_in_state]

    #----------------------------------------

    vh  =  @view x[ vh_Idx_in_state]
    
    θh  =  @view x[ θh_Idx_in_state]
    
    gens_i_d =  @view x[id_Idx_in_state]
    
    gens_i_q =  @view x[iq_Idx_in_state]

    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        vh[
            non_gens_nodes_idx ]

    non_gens_nodes_θh =
        θh[
            non_gens_nodes_idx ]
    
    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_v_ref_Idx]
    
    p_order =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_p_order_Idx]
    
    P_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Png_Idx]
    
    Q_non_gens =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end


    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]
    
    #----------------------------------------
    #----------------------------------------
    # using per plant dynamics 
    #----------------------------------------
    #----------------------------------------

    # (per_plant_gov_dyn_fun_sym_name,
    #  per_plant_generic_model_para,
    #  per_plant_generic_model_kwd_para) =
    #      get_generic_namedtuple_per_plant_para_wt_kwd_para(
    #          ωref0_vref0_porder0_id_iq_vh;
    #          kwd_para = dae_plants_kwd_para )

    # for (a_res_state_vars_view,
    #      a_du_state_vars_view,
    #      a_u_state_vars_view,
         
    #      a_gov_dyn_fun_sym_name,
    #      a_ωref0_vref0_porder0_id_iq_vh,
    #      a_per_kwd_para) in zip(
    #          res_states_vars_views,
    #          du_states_vars_views,
    #          u_states_vars_views,
             
    #          per_plant_gov_dyn_fun_sym_name,
    #          per_plant_generic_model_para,
    #          per_plant_generic_model_kwd_para )
        
    #     if a_gov_dyn_fun_sym_name == :nothing
            
    #         dae_a_SC_plant_generic_model_func!(
    #             a_res_state_vars_view,
    #             a_du_state_vars_view,
    #             a_u_state_vars_view,
    #             a_ωref0_vref0_porder0_id_iq_vh,
    #             t;
    #             kwd_para =
    #                 a_per_kwd_para )
    #     else

    #         dae_a_SM_plant_generic_model_func!(
    #             a_res_state_vars_view,
    #             a_du_state_vars_view,
    #             a_u_state_vars_view,
    #             a_ωref0_vref0_porder0_id_iq_vh,
    #             t;
    #             kwd_para =
    #                 a_per_kwd_para )
            
    #     end
        

    # end

    
    # #----------------------------------------
    # #----------------------------------------
    # # using aggregate plants dynamics 
    # #----------------------------------------
    # #----------------------------------------

    dae_gens_plants_generic_model_func!(
        res_state_vars_view,
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            dae_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]

    # dae_algebraic_generic_model_func!
    
    algebraic_generic_pf_ΔPQ_mismatch!(
        res_algebraic_vars_view,
        u_algebraic_vars_view,
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
        kwd_para =
            algebraic_generic_model_kwd_para  )
    
    return nothing

end


# ---------------------------------------------------
# mass matrix
# ---------------------------------------------------



function mm_ode_generic_system_model_dynamics!(
    dx,
    x,
    model_dynamics_para,
    t;
    kwd_para =
        generic_model_dynamics_kwd_para  )

    #----------------------------------------


    (;generic_model_dynamics_para,
     plants_cb_paras_switches ) =
         model_dynamics_para
    
    #----------------------------------------

    (;ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx, #

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs, #
     dyn_pf_fun_kwd_net_idxs, #

     # pf_generic_gens_para, #

     ode_plants_kwd_para,

     algebraic_generic_model_sol_kwd_para ) =
         NamedTupleTools.select(
             kwd_para,
             (:ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              # :pf_generic_gens_para,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_sol_kwd_para))


    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))

    #----------------------------------------

    (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    # -------------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))  

    #----------------------------------------

    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------
    
    du_state_vars_view =
        @view dx[state_var_Idx_in_state]
    
    du_algebraic_vars_view =
        @view dx[algebraic_var_Idx_in_state]
    
    u_state_vars_view =
        @view x[state_var_Idx_in_state]
    
    u_algebraic_vars_view =
        @view x[algebraic_var_Idx_in_state]

    #----------------------------------------

    # ode_dx_views   = @view  dx[state_var_Idx_in_state]

    # ode_x_views    =  @view x[state_var_Idx_in_state]
    
    # algeb_dx_views =  @view dx[algebraic_var_Idx_in_state]

    # algeb_x_views  =  @view  x[algebraic_var_Idx_in_state]
    
    #----------------------------------------

    # state_vars = x[state_var_Idx_in_state]
    
    # vh  = x[vh_Idx_in_state]
    
    # θh  = x[θh_Idx_in_state]
    
    # gens_i_d = x[id_Idx_in_state]
    
    # gens_i_q = x[iq_Idx_in_state]

    
    vh  = u_algebraic_vars_view[
        dyn_pf_vh_Idxs]
    
    θh  = u_algebraic_vars_view[
        dyn_pf_θh_Idxs]
    
    gens_i_d = u_algebraic_vars_view[
        dyn_pf_id_Idxs]
    
    gens_i_q = u_algebraic_vars_view[
        dyn_pf_iq_Idxs]
         
    #----------------------------------------
    
    # gens_δ = x[δ_idx_in_state]
        
    # gens_ed_dash = x[ed_dash_idx_in_state]
    
    # gens_eq_dash = x[eq_dash_idx_in_state]

    
    gens_δ = u_state_vars_view[
        δ_idx_in_state]
        
    gens_ed_dash = u_state_vars_view[
        ed_dash_idx_in_state]
    
    gens_eq_dash = u_state_vars_view[
        eq_dash_idx_in_state]
    
    #----------------------------------------
    
    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    non_gens_nodes_vh =
        vh[ non_gens_nodes_idx ]

    non_gens_nodes_θh =
        θh[ non_gens_nodes_idx ]
            
    #----------------------------------------
    #----------------------------------------
    
    ω_ref =
        generic_model_dynamics_para[
            dyn_ω_ref_Idx]
    
    v_ref = generic_model_dynamics_para[
        dyn_v_ref_Idx]
    
    p_order = generic_model_dynamics_para[
        dyn_p_order_Idx]
    
    P_non_gens = generic_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = generic_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            generic_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    
    #----------------------------------------
    # ode parameters and kwd para
    #----------------------------------------
    
    gens_ωref0_vref0_porder0_id_iq_vh_para = [
        [a_ωref0, a_vref, a_porder0, a_id, a_iq, a_vh ]
        for (a_ωref0,
             a_vref,
             a_porder0,
             a_id,
             a_iq,
             a_vh) in
            zip( ω_ref,
                 v_ref,
                 p_order,
                 gens_i_d,
                 gens_i_q,
                 gens_vh) ]
    
    #----------------------------------------        
    #----------------------------------------
    # ode 
    #----------------------------------------
    #----------------------------------------

    # # ode_gens_plants_generic_model_func!(
    # #     dx[state_var_Idx_in_state],
    # #     x[state_var_Idx_in_state],
    # #     ωref0_vref0_porder0_id_iq_vh,
    # #     t;
    # #     kwd_para =
    # #         plants_kwd_para )

    
    # ode_gens_plants_generic_model_func!(
    #     ode_dx_views,
    #     ode_x_views,
    #     ωref0_vref0_porder0_id_iq_vh,
    #     t;
    #     kwd_para =
    #         plants_kwd_para )

    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh]
    
    ode_gens_plants_generic_model_func!(
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            ode_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    # pf_algebraic_generic_model_sol(
    #     dx[algebraic_var_Idx_in_state],
    #     x[algebraic_var_Idx_in_state],
    #     δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    #     kwd_para =
    #         algebraic_generic_model_sol_kwd_para )

    
    # ode_pf_algebraic_generic_model_sol(
    #     du_algebraic_vars_view,
    #     u_algebraic_vars_view,
    #     δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    #     kwd_para =
    #         algebraic_generic_model_sol_kwd_para )
    
    algebraic_generic_pf_ΔPQ_mismatch_sol(
        du_algebraic_vars_view,
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
        kwd_para =
            algebraic_generic_model_sol_kwd_para )
    
    # ------------------------------------
    
    return nothing

end



function mm_ode_generic_system_model_by_funcs_dynamics!(
    dx,
    x,
    model_dynamics_para,
    t;
    kwd_para =
        generic_model_dynamics_kwd_para  )

    #----------------------------------------


    (;generic_model_dynamics_para,
      plants_cb_paras_switches ) =
         model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     ode_plants_kwd_para,

     algebraic_generic_model_sol_kwd_para,
     algebraic_generic_model_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_sol_kwd_para,
              :algebraic_generic_model_kwd_para))


    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))

    #----------------------------------------

    (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))
   
    idx_range =
        first(dyn_pf_vh_Idxs):last(dyn_pf_iq_Idxs)
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    # -------------------------------------

    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------
    
    du_algebraic_vars_view =
        @view dx[algebraic_var_Idx_in_state]
    
    u_algebraic_vars_view =
        @view x[algebraic_var_Idx_in_state]

    #----------------------------------------
    
    du_state_vars_view =
        @view dx[state_var_Idx_in_state]
    
    u_state_vars_view =
        @view x[state_var_Idx_in_state]
    #----------------------------------------
    
    # gens_δ = u_state_vars_view[
    #     δ_idx_in_state]
        
    # gens_ed_dash = u_state_vars_view[
    #     ed_dash_idx_in_state]
    
    # gens_eq_dash = u_state_vars_view[
    #     eq_dash_idx_in_state]

    
    gens_δ = x[
        δ_idx_in_state]
        
    gens_ed_dash = x[
        ed_dash_idx_in_state]
    
    gens_eq_dash = x[
        eq_dash_idx_in_state]
    
    #----------------------------------------
    
    # vh = @view x[vh_Idx_in_state]
    
    # θh = @view x[θh_Idx_in_state]
    
    # gens_i_d = @view x[id_Idx_in_state]
    
    # gens_i_q = @view x[ iq_Idx_in_state]

    
    vh = x[vh_Idx_in_state]
    
    θh = x[θh_Idx_in_state]
    
    gens_i_d = x[id_Idx_in_state]
    
    gens_i_q = x[ iq_Idx_in_state]
    
    #----------------------------------------
        
    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    # non_gens_nodes_vh =
    #     vh[ non_gens_nodes_idx ]

    # non_gens_nodes_θh =
    #     θh[ non_gens_nodes_idx ]
         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        generic_model_dynamics_para[
            dyn_ω_ref_Idx]
    
    v_ref = generic_model_dynamics_para[
        dyn_v_ref_Idx]
    
    p_order = generic_model_dynamics_para[
        dyn_p_order_Idx]
    
    P_non_gens = generic_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = generic_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            generic_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]
    
    #----------------------------------------
    #----------------------------------------
    # using per plant dynamics 
    #----------------------------------------
    #----------------------------------------

    (per_plant_gov_dyn_fun_sym_name,
     per_plant_generic_model_para,
     per_plant_generic_model_kwd_para) =
         get_generic_namedtuple_per_plant_para_wt_kwd_para(
             ωref0_vref0_porder0_id_iq_vh;
             kwd_para = ode_plants_kwd_para )

    for (
         a_du_state_vars_view,
         a_u_state_vars_view,
         
         a_gov_dyn_fun_sym_name,
         a_ωref0_vref0_porder0_id_iq_vh,
         a_per_kwd_para) in zip(
             
             du_states_vars_views,
             u_states_vars_views,
             
             per_plant_gov_dyn_fun_sym_name,
             per_plant_generic_model_para,
             per_plant_generic_model_kwd_para )
        
        if a_gov_dyn_fun_sym_name == :nothing
            
            ode_a_SC_plant_generic_model_func!(
                a_du_state_vars_view,
                a_u_state_vars_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    a_per_kwd_para )
        else

            ode_a_SM_plant_generic_model_func!(
                a_du_state_vars_view,
                a_u_state_vars_view,
                a_ωref0_vref0_porder0_id_iq_vh,
                t;
                kwd_para =
                    a_per_kwd_para )
            
        end
        

    end

    
    # #----------------------------------------
    # #----------------------------------------
    # # using aggregate plants dynamics 
    # #----------------------------------------
    # #----------------------------------------

    # ode_gens_plants_generic_model_func!(
    #     du_state_vars_view,
    #     u_state_vars_view,
    #     ωref0_vref0_porder0_id_iq_vh,
    #     t;
    #     kwd_para =
    #         ode_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    # pf_algebraic_generic_model_sol(
    #     dx[algebraic_var_Idx_in_state],
    #     x[algebraic_var_Idx_in_state],
    #     δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    #     kwd_para =
    #         algebraic_generic_model_sol_kwd_para )


    vh_θh_id_iq =
        [vh;
         θh;
         gens_i_d;
         gens_i_q]

    # # ode_pf_algebraic_generic_model_sol_ext
    
    # pf_sol = algebraic_generic_pf_ΔPQ_mismatch_sol(
    #     vh_θh_id_iq,
    #     δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    #     kwd_para =
    #         algebraic_generic_model_sol_kwd_para )
    
        
    # du_algebraic_vars_view .=
    #     pf_sol[ idx_range ] -  u_algebraic_vars_view

    
    algebraic_generic_pf_ΔPQ_mismatch!(
        du_algebraic_vars_view,
        # x[algebraic_var_Idx_in_state],
        vh_θh_id_iq,
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
        kwd_para =
            algebraic_generic_model_kwd_para  )
    
    # ------------------------------------
    
    return nothing

end



function mm_ode_generic_system_dynamics_by_ode_pf_funcs!(
    dx,
    x,
    ω_ref_v_ref_p_order_Png_Qng_Pll_Qll,
    t;
    kwd_para =
        generic_model_dynamics_kwd_para  )

    #----------------------------------------


    # (;generic_model_dynamics_para,
    #   plants_cb_paras_switches ) =
    #      model_dynamics_para
    
    #----------------------------------------

    (;gens_nodes_idx,
     ωs,
     loc_load_exist,
     state_vars_idx,

     # id_iq_pg_vh_Idx,

     dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

     gens_state_vars_idx_in_state,
     state_vars_and_i_dq_Idx_in_state,
     state_algebraic_vars_Idx_in_state,

     dyn_pf_flat_vh_flat_θh_id_iq_Idx,

     dyn_pf_fun_kwd_n2s_idxs, 
     dyn_pf_fun_kwd_net_idxs, 

     ode_plants_kwd_para,

     algebraic_generic_model_sol_kwd_para,
     algebraic_generic_model_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (:gens_nodes_idx,
              :ωs,
              :loc_load_exist,
              :state_vars_idx,

              # :id_iq_pg_vh_Idx,

              :dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,

              :gens_state_vars_idx_in_state,
              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :dyn_pf_flat_vh_flat_θh_id_iq_Idx,

              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,
              
              :ode_plants_kwd_para,

              :algebraic_generic_model_sol_kwd_para,
              :algebraic_generic_model_kwd_para))


    #----------------------------------------
    
    (;dyn_ω_ref_Idx,
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_ω_ref_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (:dyn_ω_ref_Idx,
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    #----------------------------------------
    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))

    #----------------------------------------

    (;
     dyn_pf_vh_Idxs,
     dyn_pf_θh_Idxs,       
     dyn_pf_id_Idxs,
     dyn_pf_iq_Idxs) =
         NamedTupleTools.select(
             dyn_pf_flat_vh_flat_θh_id_iq_Idx,
             (:dyn_pf_vh_Idxs,
              :dyn_pf_θh_Idxs,       
              :dyn_pf_id_Idxs,
              :dyn_pf_iq_Idxs ))
   
    idx_range =
        first(dyn_pf_vh_Idxs):last(dyn_pf_iq_Idxs)
    
    #----------------------------------------
    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    # -------------------------------------

    du_states_vars_views =
        [ view(dx, idx) for idx in state_vars_idx]

    u_states_vars_views =
        [ view(x, idx) for idx in state_vars_idx]

    #----------------------------------------
    
    du_algebraic_vars_view =
        @view dx[algebraic_var_Idx_in_state]
    
    u_algebraic_vars_view =
        @view x[algebraic_var_Idx_in_state]

    #----------------------------------------
    
    du_state_vars_view =
        @view dx[state_var_Idx_in_state]
    
    u_state_vars_view =
        @view x[state_var_Idx_in_state]
    #----------------------------------------
    
    # gens_δ = u_state_vars_view[
    #     δ_idx_in_state]
        
    # gens_ed_dash = u_state_vars_view[
    #     ed_dash_idx_in_state]
    
    # gens_eq_dash = u_state_vars_view[
    #     eq_dash_idx_in_state]

    
    gens_δ = x[
        δ_idx_in_state]
        
    gens_ed_dash = x[
        ed_dash_idx_in_state]
    
    gens_eq_dash = x[
        eq_dash_idx_in_state]
    
    #----------------------------------------
    
    # vh = @view x[vh_Idx_in_state]
    
    # θh = @view x[θh_Idx_in_state]
    
    # gens_i_d = @view x[id_Idx_in_state]
    
    # gens_i_q = @view x[ iq_Idx_in_state]

    
    vh = x[vh_Idx_in_state]
    
    θh = x[θh_Idx_in_state]
    
    gens_i_d = x[id_Idx_in_state]
    
    gens_i_q = x[ iq_Idx_in_state]
    
    #----------------------------------------
        
    gens_vh =
        vh[ gens_nodes_idx ]

    gens_θh =
        θh[ gens_nodes_idx ]

    # non_gens_nodes_vh =
    #     vh[ non_gens_nodes_idx ]

    # non_gens_nodes_θh =
    #     θh[ non_gens_nodes_idx ]
         
    #----------------------------------------    
    #----------------------------------------
    
    ω_ref =
        ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
            dyn_ω_ref_Idx]
    
    v_ref = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_v_ref_Idx]
    
    p_order = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_p_order_Idx]
    
    P_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Png_Idx]
    
    Q_non_gens = ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            ω_ref_v_ref_p_order_Png_Qng_Pll_Qll[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------
    
    ωref0_vref0_porder0_id_iq_vh =
        [ω_ref;
         v_ref;
         p_order;
         gens_i_d;
         gens_i_q;
         gens_vh ]
    
    #----------------------------------------
    #----------------------------------------
    # using plant dynamics 
    #----------------------------------------
    #----------------------------------------
    
    ode_gens_plants_generic_model_func!(
        du_state_vars_view,
        u_state_vars_view,
        ωref0_vref0_porder0_id_iq_vh,
        t;
        kwd_para =
            ode_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    # pf_algebraic_generic_model_sol(
    #     dx[algebraic_var_Idx_in_state],
    #     x[algebraic_var_Idx_in_state],
    #     δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    #     kwd_para =
    #         algebraic_generic_model_sol_kwd_para )

    vh_θh_id_iq =
        [vh;
         θh;
         gens_i_d;
         gens_i_q]

    # # ode_pf_algebraic_generic_model_sol_ext
    
    # pf_sol = algebraic_generic_pf_ΔPQ_mismatch_sol(
    #     vh_θh_id_iq,
    #     δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
    #     kwd_para =
    #         algebraic_generic_model_sol_kwd_para )
    
        
    # du_algebraic_vars_view .=
    #     pf_sol[ idx_range ] -  u_algebraic_vars_view
        
    algebraic_generic_pf_ΔPQ_mismatch!(
        du_algebraic_vars_view,
        # x[algebraic_var_Idx_in_state],
        vh_θh_id_iq,
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
        kwd_para =
            algebraic_generic_model_kwd_para  )
    
    # ------------------------------------
    
    return nothing

end

#---------------------------------------------------
####################################################
#---------------------------------------------------


# ------------------------------------------------------
# ------------------------------------------------------
# Partitioned-explicit (PE) method.
# ------------------------------------------------------
# ------------------------------------------------------


#-----------------------------------------------------
# generic model 
#-----------------------------------------------------


function ode_a_gen_generic_model_func!(
    dx, x,
    p,
    t;
    kwd_para =
        kwd_para )

    # (vh, θh, i_d, i_q, V_ref) = p
    
    (vh, θh, V_ref, Tm) = p

    (gen_para,
     avr_para,
     ωs) =
         kwd_para
    
    (H,
     D,
     ra,
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gen_para,
             (:H,
              :D,
              :ra,
              :X_d,
              :X_q,

              :X_d_dash,
              :X_q_dash,

              :T_d_dash,
              :T_q_dash))
    
    (Ka, Ta,
     Ke, Te,
     Kf, Tf,
     Ae, Be,
     Tr ) =
         NamedTupleTools.select(
             avr_para,
             (:Ka, :Ta,
              :Ke, :Te,
              :Kf, :Tf,
              :Ae, :Be,
              :Tr))
        
    δ, ω, ed_dash, eq_dash, E_fd, R_f, V_R = x

    # Zdq = Z_dq(ra, X_d_dash, X_q_dash)
    # inv_Zdq = invZ_dq(ra, X_d_dash, X_q_dash)

    i_dq = invZ_dq(ra, X_d_dash, X_q_dash) * [
            ed_dash - vh * sin(δ - θh), eq_dash -
                vh * cos(δ - θh)]
    
    i_d = i_dq[1]

    i_q = i_dq[2]

    # if abs(dx[2]) < Δω_toleraance
        
    #     Δω = ωs - ω
        
    #     dx[1] = ω - ωs +  Δω
        
    # else
        
    #     dx[1] = ω - ωs 
        
    # end        

    dx[1] = ω - ωs
    
    dx[2] = ( ωs/(2*H)) * (Tm - ed_dash * i_d - eq_dash *
        i_q - (X_q_dash - X_d_dash) * i_d * i_q -
        D * (ω - ωs))
    
    dx[3] = (1/T_q_dash) *
        ( (X_q - X_q_dash) * i_q  - ed_dash)

    dx[4] = (1/T_d_dash) *
        (E_fd - (X_d - X_d_dash) * i_d - eq_dash)

    dx[5] = (1/Te) * (V_R - (Ke + Sevf(Ae, Be, E_fd)) *
        E_fd  )
    
    dx[6] = (1/Tf) * (Kf * E_fd /Tf -  R_f )

    
    dxs[7] = (1/Ta) * ( Ka * R_f - Ka * Kf * E_fd /Tf +
        Ka * (V_ref - vh) - V_R)
    
    return nothing    

end


function ode_generic_model_func!(
    dx, x,
    gens_vh_θh_V_ref_Tm_para,
    t;
    kwd_para =
        kwd_para )

    (;gens_para,
     avrs_para,
     ωs,
     state_vars_idx,
     ode_vh_θh_V_ref_Tm_Idx) =
         kwd_para

    (;ode_vh_Idx,
     ode_θh_Idx,
     ode_V_ref_Idx,
     ode_Tm_Idx ) =
         ode_vh_θh_V_ref_Tm_Idx

    #----------------------------------------
    
    gens_vh =
        gens_vh_θh_V_ref_Tm_para[
            ode_vh_Idx]

    gens_θh =
        gens_vh_θh_V_ref_Tm_para[
            ode_θh_Idx]
    
    V_ref =
        gens_vh_θh_V_ref_Tm_para[
            ode_V_ref_Idx]
    
    Tm =
        gens_vh_θh_V_ref_Tm_para[
            ode_Tm_Idx]

    #----------------------------------------
    
    (H,
     D,

     ra,
     
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
        gens_para
    
    (Ka, Ta,
     Ke, Te,
     Kf, Tf,
     Ae, Be,
     Tr ) =
        avrs_para
    

    for (gen_vh,
         gen_θh,
         gen_V_ref, gen_Tm,
         gen_H, gen_D, gen_ra,
         gen_X_d, gen_X_q,
         gen_X_d_dash, gen_X_q_dash,
         gen_T_d_dash, gen_T_q_dash,
         avr_Ka, avr_Ta,
         avr_Ke, avr_Te,
         avr_Kf, avr_Tf,
         avr_Ae, avr_Be,
         avr_Tr,
         state_var_idx) in
        zip( gens_vh,
             gens_θh, 
             V_ref, Tm,
             H, D, ra,
             X_d, X_q,
             X_d_dash, X_q_dash,
             T_d_dash, T_q_dash,
             Ka, Ta,
             Ke, Te,
             Kf, Tf,
             Ae, Be,
             Tr,
             state_vars_idx )

        gen_vh_θh_V_ref_Tm_para = [gen_vh;
               gen_θh;
               gen_V_ref;
               gen_Tm] 

        # (vh, θh, V_ref, Tm)
        
        gen_para = (gen_H, gen_D, gen_ra,
                    gen_X_d, gen_X_q,
                    gen_X_d_dash, gen_X_q_dash,
                    gen_T_d_dash, gen_T_q_dash, )
        
        avr_para = (avr_Ka, avr_Ta,
                    avr_Ke, avr_Te,
                    avr_Kf, avr_Tf,
                    avr_Ae, avr_Be,
                    avr_Tr)
        
        a_kwd_para = (;gen_para,
                      avr_para,
                      ωs)

        ode_a_gen_generic_model_func!(
            dx[state_var_idx],
            x[state_var_idx],
            gen_vh_θh_V_ref_Tm_para,
            t;
            kwd_para =
                a_kwd_para  )

    end

    return nothing    

end




function ode_a_gen_generic_model_by_ext_idq_func!(
    dx, x,
    gen_vh_id_iq_V_ref_Tm_para,
    t;
    kwd_para =
        gen_ode_kwd_para )
    
    vh, i_d, i_q, V_ref, Tm =
        gen_vh_id_iq_V_ref_Tm_para

    
    (gen_para,
     avr_para,
     ωs) =
         kwd_para

    
    (H,
     D,
     ra,
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gen_para,
             (:H,
              :D,
              :ra,
              :X_d,
              :X_q,

              :X_d_dash,
              :X_q_dash,

              :T_d_dash,
              :T_q_dash))
    
    (Ka, Ta,
     Ke, Te,
     Kf, Tf,
     Ae, Be,
     Tr ) =
         NamedTupleTools.select(
             avr_para,
             (:Ka, :Ta,
              :Ke, :Te,
              :Kf, :Tf,
              :Ae, :Be,
              :Tr))

    
    δ, ω, ed_dash, eq_dash, E_fd, R_f, V_R = x

    dx[1] = ω - ωs
    
    dx[2] = ( ωs/(2*H)) * (Tm - ed_dash * i_d - eq_dash *
        i_q - (X_q_dash - X_d_dash) * i_d * i_q -
        D * (ω - ωs))

    dx[3] = (1/T_q_dash) *
        ( (X_q - X_q_dash) * i_q  - ed_dash)

    dx[4] = (1/T_d_dash) *
        (E_fd - (X_d - X_d_dash) * i_d  - eq_dash)

    dx[5] = (1/Te) * (V_R - (Ke + Sevf(Ae, Be, E_fd)) *
        E_fd  )
    
    dx[6] = (1/Tf) * ((Kf/Tf) *  E_fd  -  R_f )

    
    dx[7] = (1/Ta) * ( Ka * R_f  - (Ka * Kf/Tf) * E_fd +
        Ka * (V_ref - vh) - V_R)
    
    return nothing    

end


#----------------------------------------
                                         

function ode_gens_generic_model_by_ext_idq_func!(
    dx, x,
    gens_vh_id_iq_V_ref_Tm_para,
    t;
    kwd_para =
        gens_ode_kwd_para  )

    (;generic_gens_para,
     generic_avrs_para,
     ωs,
     state_vars_idx,
     dyn_5_gens_type_paras_Idx ) =
         NamedTupleTools.select(
             kwd_para,
             (:generic_gens_para,
              :generic_avrs_para,
              :ωs,
              :state_vars_idx,
              :dyn_5_gens_type_paras_Idx))

    # (ode_vh_Idx,
    #  ode_id_Idx,
    #  ode_iq_Idx,
    #  ode_V_ref_Idx,
    #  ode_Tm_Idx ) =
    #      NamedTupleTools.select(
    #          dyn_5_gens_type_paras_Idx ,
    #          (:gens_para_1_Idxs,
    #           :gens_para_2_Idxs,       
    #           :gens_para_3_Idxs,
    #           :gens_para_4_Idxs,       
    #           :gens_para_5_Idxs))

    ode_vh_Idx =
        dyn_5_gens_type_paras_Idx.gens_para_1_Idx
    
    ode_id_Idx =
        dyn_5_gens_type_paras_Idx.gens_para_2_Idx
    
    ode_iq_Idx =
        dyn_5_gens_type_paras_Idx.gens_para_3_Idx
    
    ode_V_ref_Idx =
        dyn_5_gens_type_paras_Idx.gens_para_4_Idx
    
    ode_Tm_Idx =
        dyn_5_gens_type_paras_Idx.gens_para_5_Idx
    
    #----------------------------------------
    
    gens_vh =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_vh_Idx]
    
    gens_i_d =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_id_Idx]
    
    gens_i_q =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_iq_Idx]
    
    V_ref =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_V_ref_Idx]
    
    Tm =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_Tm_Idx]

    #----------------------------------------

    for (state_var_idx,
         gen_vh,
         gen_i_d,
         gen_i_q,
         gen_V_ref,
         gen_Tm,
         gen_para,
         avr_para) in
        zip(state_vars_idx,
            gens_vh,
            gens_i_d,
            gens_i_q,
            V_ref,
            Tm,
            generic_gens_para,
            generic_avrs_para)

        gen_vh_id_iq_V_ref_Tm_para =
            [gen_vh;
             gen_i_d;
             gen_i_q;
             gen_V_ref;
             gen_Tm] 
        
        a_kwd_para = (;gen_para,
                      avr_para,
                      ωs)

        ode_a_gen_generic_model_by_ext_idq_func!(
            dx[state_var_idx],
            x[state_var_idx],
            gen_vh_id_iq_V_ref_Tm_para,
            t;
            kwd_para =
                a_kwd_para  )

    end

    return nothing    

end


#-----------------------------------------------------
# generic model system dynamics
#-----------------------------------------------------


function ode_generic_model_dynamics!(
    dx,
    x,
    dae_SC_generic_model_dynamics_para,
    t;
    kwd_para =
        dae_SC_generic_model_dynamics_kwd_para  )

    #----------------------------------------

    (;ωs,
     gens_nodes_idx,         
     loc_load_exist,
     state_vars_idx,

     state_vars_and_i_dq_Idx_in_state,     
     state_algebraic_vars_Idx_in_state,

     gens_state_vars_idx_in_state,
     
     dyn_4_gens_type_paras_Idx,
     dyn_5_gens_type_paras_Idx,

     dyn_2_gens_paras_wt_Png_Qng_Pll_Qll_Idx,
     dyn_3_gens_paras_wt_Png_Qng_Pll_Qll_Idx,

     generic_gens_para,
     generic_govs_para,
     generic_avrs_para,

     algebraic_generic_model_sol_kwd_para,
     gens_ode_kwd_para ) =
         NamedTupleTools.select(
             kwd_para,
             (:ωs,
              :gens_nodes_idx,         
              :loc_load_exist,
              :state_vars_idx,

              :state_vars_and_i_dq_Idx_in_state,
              :state_algebraic_vars_Idx_in_state,

              :gens_state_vars_idx_in_state,

              :dyn_4_gens_type_paras_Idx,
              :dyn_5_gens_type_paras_Idx,

              :dyn_2_gens_paras_wt_Png_Qng_Pll_Qll_Idx,
              :dyn_3_gens_paras_wt_Png_Qng_Pll_Qll_Idx,

              :generic_gens_para,
              :generic_govs_para,
              :generic_avrs_para,

              :algebraic_generic_model_sol_kwd_para,
              :gens_ode_kwd_para ) )

    #----------------------------------------

    (;
     # state_var_Idx_in_state,
     vh_Idx_in_state,
     θh_Idx_in_state,
     id_Idx_in_state,
     iq_Idx_in_state) =
         NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
             (
                 # :state_var_Idx_in_state,
                 :vh_Idx_in_state,
                 :θh_Idx_in_state,
                 :id_Idx_in_state,
                 :iq_Idx_in_state))

    
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state) =
        state_algebraic_vars_Idx_in_state

    
    (dyn_V_ref_Idx,
     dyn_Tm_Idx,
     
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_2_gens_paras_wt_Png_Qng_Pll_Qll_Idx,
             (:gens_para_1_Idx,
              :gens_para_2_Idx,
              :dyn_P_non_gens_Idxs,
              :dyn_Q_non_gens_Idxs,
              :dyn_P_gens_loc_load_Idxs,
              :dyn_Q_gens_loc_load_Idxs))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             SC_generic_model_vars_wt_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    #----------------------------------------

    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             SC_generic_model_states_comp_idxs_in_Idx,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state
              ))
    
    # -------------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))  

    #----------------------------------------
    
   #  (ra,
   #   X_d,
   #   X_q,
   #   X_d_dash,
   #   X_q_dash ) =
   #       NamedTupleTools.select(
   #          pf_gens_para,
   #           (:ra,
   #            :X_d,
   #            :X_q,
   #            :X_d_dash,
   #            :X_q_dash ) )
    
    
   #  #----------------------------------------

   # (Ynet,
   #  nodes_idx_with_adjacent_nodes_idx) =
   #      NamedTupleTools.select(
   #          Ynet_wt_nodes_idx_wt_adjacent_nodes,
   #          (:Ynet,
   #           :nodes_idx_with_adjacent_nodes_idx)) 
    
    #----------------------------------------
    #----------------------------------------

    state_vars = x[state_var_Idx_in_state]
    
    vh  = x[vh_Idx_in_state]
    
    θh  = x[θh_Idx_in_state]
    
    gens_i_d = x[id_Idx_in_state]
    
    gens_i_q = x[iq_Idx_in_state]
    
    #----------------------------------------
    
    gens_δ = x[δ_idx_in_state]
        
    gens_ed_dash = x[ed_dash_idx_in_state]
    
    gens_eq_dash = x[eq_dash_idx_in_state]
        
    #----------------------------------------
    
    gens_vh =
        vh[gens_nodes_idx]

    gens_θh =
        θh[gens_nodes_idx]

    non_gens_nodes_vh =
        vh[non_gens_nodes_idx]

    non_gens_nodes_θh =
        θh[non_gens_nodes_idx]
    
    #----------------------------------------
        
    V_ref = dae_SC_generic_model_dynamics_para[
        dyn_V_ref_Idx]
    
    Tm = dae_SC_generic_model_dynamics_para[
        dyn_Tm_Idx]
    
    P_non_gens = dae_SC_generic_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = dae_SC_generic_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            dae_SC_generic_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            dae_SC_generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------    
    # single 
    #----------------------------------------

    for (state_var_idx,
         a_gen_vh,
         a_gen_i_d,
         a_gen_i_q,
         a_v_ref,
         a_τm,
         a_gen_para,
         a_avr_para) in
        zip(state_vars_idx,
            gens_vh,
            gens_i_d,
            gens_i_q,
            V_ref,
            Tm,
            generic_gens_para,
            generic_avrs_para)

        gen_ode_kwd_para =
            (a_gen_para,
             a_avr_para,
             ωs)
        
        # gen_vh_id_iq_V_ref_Tm_para =
        #     [a_gen_vh, a_gen_i_d, a_gen_i_q,
        #      a_v_ref, a_τm ]
        
        ode_a_gen_generic_model_by_ext_idq_func!(
            dx[state_var_idx],
            x[state_var_idx],
            [a_gen_vh, a_gen_i_d,
             a_gen_i_q, a_v_ref, a_τm ],
            t;
            kwd_para =
                gen_ode_kwd_para )
    end
    
    # #----------------------------------------
    # # ode parameters and kwd para
    # #----------------------------------------
    
    # gens_ode_kwd_para = [
    #     (gen_para =
    #         (gen_H, gen_D,
    #          gen_X_d, gen_X_q,
    #          gen_X_d_dash, gen_X_q_dash,
    #          gen_T_d_dash, gen_T_q_dash ),
    #      avr_para =
    #          (avr_Ka, avr_Ta,
    #           avr_Ke, avr_Te,
    #           avr_Kf, avr_Tf,
    #           avr_Ae, avr_Be,
    #           avr_Tr),
    #      ωs = ωs )
    #     for (gen_H, gen_D,
    #          gen_X_d, gen_X_q,
    #          gen_X_d_dash, gen_X_q_dash,
    #          gen_T_d_dash, gen_T_q_dash,
    #          avr_Ka, avr_Ta,
    #          avr_Ke, avr_Te,
    #          avr_Kf, avr_Tf,
    #          avr_Ae, avr_Be,
    #          avr_Tr) in
    #         zip(H, D,
    #          X_d, X_q,
    #          X_d_dash, X_q_dash,
    #          T_d_dash, T_q_dash,
    #          Ka, Ta,
    #          Ke, Te,
    #          Kf, Tf,
    #          Ae, Be,
    #          Tr )]    
    
    # gens_vh_id_iq_V_ref_Tm_para = [
    #     [a_vh, a_id, a_iq, a_vref, a_tm]
    #     for (a_vh, a_id, a_iq, a_vref, a_tm) in
    #         zip(gens_vh,
    #         gens_i_d,
    #         gens_i_q,
    #         V_ref,
    #             Tm) ]
    
    # #----------------------------------------        
    # #----------------------------------------
    # # ode 
    # #----------------------------------------
    # #----------------------------------------
    
    # ode_dx_views = [view(dx, idx)
    #                 for idx in state_vars_idx ]
    
    # ode_x_views  = [view( x, idx)
    #                 for idx in state_vars_idx ]
    
    # for (ode_res, ode_dx, ode_x,
    #      gen_vh_id_iq_V_ref_Tm_para,
    #      gen_ode_kwd_para) in
    #     zip(ode_res_views,
    #         ode_dx_views,
    #         ode_x_views,
    #         gens_vh_id_iq_V_ref_Tm_para,
    #         gens_ode_kwd_para)
        
    #     ode_a_gen_generic_model_by_ext_idq_func!(
    #         ode_dx,
    #         ode_x,
    #         gen_vh_id_iq_V_ref_Tm_para,
    #         t;
    #         kwd_para =
    #             gen_ode_kwd_para )
        
    # end

    
    
    return nothing

end


# ------------------------------------------------------
# ------------------------------------------------------
# Simultaneous implicit DAE (SI) method
# ------------------------------------------------------
# ------------------------------------------------------

#-----------------------------------------------------
# generic model 
#-----------------------------------------------------


function dae_a_gen_generic_model_func!(
    res,
    dx, x,
    p,
    t;
    kwd_para =
        kwd_para )

    # (vh, θh, i_d, i_q, V_ref) = p
    
    (vh, θh, V_ref, Tm) = p

    (gen_para,
     avr_para,
     ωs) =
         kwd_para
    
    (H,
     D,
     ra,
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gen_para,
             (:H,
              :D,
              :ra,
              :X_d,
              :X_q,

              :X_d_dash,
              :X_q_dash,

              :T_d_dash,
              :T_q_dash))
    
    (Ka, Ta,
     Ke, Te,
     Kf, Tf,
     Ae, Be,
     Tr ) =
         NamedTupleTools.select(
             avr_para,
             (:Ka, :Ta,
              :Ke, :Te,
              :Kf, :Tf,
              :Ae, :Be,
              :Tr))
        
    δ, ω, ed_dash, eq_dash, E_fd, R_f, V_R = x

    # Zdq = Z_dq(ra, X_d_dash, X_q_dash)
    # inv_Zdq = invZ_dq(ra, X_d_dash, X_q_dash)

    i_dq = invZ_dq(ra, X_d_dash, X_q_dash) * [
            ed_dash - vh * sin(δ - θh), eq_dash -
                vh * cos(δ - θh)]
    
    i_d = i_dq[1]

    i_q = i_dq[2]

    # if abs( dx[2] ) < Δω_toleraance
        
    #     Δω = ωs - ω
        
    #     res[1] = ω - ωs - dx[1] +  Δω
        
    # else
        
    #     res[1] = ω - ωs - dx[1]
        
    # end

    res[1] = ω - ωs - dx[1]
    
    res[2] = ( ωs/(2*H)) * (Tm - ed_dash * i_d - eq_dash *
        i_q - (X_q_dash - X_d_dash) * i_d * i_q -
        D * (ω - ωs))  - dx[2]
    
    res[3] = (1/T_q_dash) *
        ( (X_q - X_q_dash) * i_q  - ed_dash) - dx[3]

    res[4] = (1/T_d_dash) *
        (E_fd - (X_d - X_d_dash) * i_d - eq_dash) - dx[4]

    res[5] = (1/Te) * (V_R - (Ke + Sevf(Ae, Be, E_fd)) *
        E_fd  ) - dx[5]
    
    res[6] = (1/Tf) * (Kf * E_fd /Tf -  R_f ) - dx[6]

    
    res[7] = (1/Ta) * ( Ka * R_f - Ka * Kf * E_fd /Tf +
        Ka * (V_ref - vh) - V_R) - dx[7]
    
    return nothing    

end


function dae_generic_model_func!(
    res,
    dx, x,
    gens_vh_θh_V_ref_Tm_para,
    t;
    kwd_para =
        kwd_para )

    (;gens_para,
     avrs_para,
     ωs,
     state_vars_idx,
     ode_vh_θh_V_ref_Tm_Idx) =
         kwd_para

    (;ode_vh_Idx,
     ode_θh_Idx,
     ode_V_ref_Idx,
     ode_Tm_Idx ) =
         ode_vh_θh_V_ref_Tm_Idx

    #----------------------------------------
    
    gens_vh =
        gens_vh_θh_V_ref_Tm_para[
            ode_vh_Idx]

    gens_θh =
        gens_vh_θh_V_ref_Tm_para[
            ode_θh_Idx]
    
    V_ref =
        gens_vh_θh_V_ref_Tm_para[
            ode_V_ref_Idx]
    
    Tm =
        gens_vh_θh_V_ref_Tm_para[
            ode_Tm_Idx]

    #----------------------------------------
    
    (H,
     D,

     ra,
     
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
        gens_para
    
    (Ka, Ta,
     Ke, Te,
     Kf, Tf,
     Ae, Be,
     Tr ) =
        avrs_para
    

    for (gen_vh,
         gen_θh,
         gen_V_ref, gen_Tm,
         gen_H, gen_D, gen_ra,
         gen_X_d, gen_X_q,
         gen_X_d_dash, gen_X_q_dash,
         gen_T_d_dash, gen_T_q_dash,
         avr_Ka, avr_Ta,
         avr_Ke, avr_Te,
         avr_Kf, avr_Tf,
         avr_Ae, avr_Be,
         avr_Tr,
         state_var_idx) in
        zip( gens_vh,
             gens_θh, 
             V_ref, Tm,
             H, D, ra,
             X_d, X_q,
             X_d_dash, X_q_dash,
             T_d_dash, T_q_dash,
             Ka, Ta,
             Ke, Te,
             Kf, Tf,
             Ae, Be,
             Tr,
             state_vars_idx )

        gen_vh_θh_V_ref_Tm_para = [gen_vh;
               gen_θh;
               gen_V_ref;
               gen_Tm] 

        # (vh, θh, V_ref, Tm)
        
        gen_para = (gen_H, gen_D, gen_ra,
                    gen_X_d, gen_X_q,
                    gen_X_d_dash, gen_X_q_dash,
                    gen_T_d_dash, gen_T_q_dash, )
        
        avr_para = (avr_Ka, avr_Ta,
                    avr_Ke, avr_Te,
                    avr_Kf, avr_Tf,
                    avr_Ae, avr_Be,
                    avr_Tr)
        
        a_kwd_para = (;gen_para,
                      avr_para,
                      ωs)

        dae_a_gen_generic_model_func!(
            res[state_var_idx],
            dx[state_var_idx],
            x[state_var_idx],
            gen_vh_θh_V_ref_Tm_para,
            t;
            kwd_para =
                a_kwd_para  )

    end

    return nothing    

end




function dae_a_gen_generic_model_by_ext_idq_func!(
    res,
    dx, x,
    gen_vh_id_iq_V_ref_Tm_para,
    t;
    kwd_para =
        gen_ode_kwd_para )
    
    vh, i_d, i_q, V_ref, Tm =
        gen_vh_id_iq_V_ref_Tm_para

    
    (gen_para,
     avr_para,
     ωs) =
         kwd_para

    
    (H,
     D,
     ra,
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gen_para,
             (:H,
              :D,
              :ra,
              :X_d,
              :X_q,

              :X_d_dash,
              :X_q_dash,

              :T_d_dash,
              :T_q_dash))
    
    (Ka, Ta,
     Ke, Te,
     Kf, Tf,
     Ae, Be,
     Tr ) =
         NamedTupleTools.select(
             avr_para,
             (:Ka, :Ta,
              :Ke, :Te,
              :Kf, :Tf,
              :Ae, :Be,
              :Tr))
    
    δ, ω, ed_dash, eq_dash, E_fd, R_f, V_R = x
    
    # if abs( dx[2] ) < Δω_toleraance
        
    #     Δω = ωs - ω
        
    #     res[1] = ω - ωs - dx[1] +  Δω

    # else
        
    #     res[1] = ω - ωs - dx[1]
        
    # end

    res[1] = ω - ωs - dx[1]
    
    res[2] = ( ωs/(2*H)) * (Tm - ed_dash * i_d - eq_dash *
        i_q - (X_q_dash - X_d_dash) * i_d * i_q -
        D * (ω - ωs))  - dx[2]

    res[3] = (1/T_q_dash) *
        ( (X_q - X_q_dash) * i_q  - ed_dash) - dx[3]

    res[4] = (1/T_d_dash) *
        (E_fd - (X_d - X_d_dash) * i_d  - eq_dash) - dx[4]

    res[5] = (1/Te) * (V_R - (Ke + Sevf(Ae, Be, E_fd)) *
        E_fd  ) - dx[5]
    
    res[6] = (1/Tf) * ((Kf/Tf) *  E_fd  -  R_f ) - dx[6]

    
    res[7] = (1/Ta) * ( Ka * R_f  - (Ka * Kf/Tf) * E_fd +
        Ka * (V_ref - vh) - V_R) - dx[7]
    
    return nothing    

end


#----------------------------------------
#----------------------------------------
                                         

function dae_gens_generic_model_by_ext_idq_func!(
    res,
    dx, x,
    gens_vh_id_iq_V_ref_Tm_para,
    t;
    kwd_para =
        kwd_para  )

    (;generic_gens_para,
     generic_avrs_para,
     ωs,
     state_vars_idx,
     generic_5_gens_type_paras_Idx ) =
         kwd_para

    (ode_vh_Idx,
     ode_id_Idx,
     ode_iq_Idx,
     ode_V_ref_Idx,
     ode_Tm_Idx ) =
         generic_5_gens_type_paras_Idx

    #----------------------------------------
    
    gens_vh =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_vh_Idx]
    
    gens_i_d =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_id_Idx]
    
    gens_i_q =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_iq_Idx]
    
    V_ref =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_V_ref_Idx]
    
    Tm =
        gens_vh_id_iq_V_ref_Tm_para[
            ode_Tm_Idx]

    #----------------------------------------

    for (state_var_idx,
         gen_para,
         avr_para) in
        zip(state_vars_idx,
            generic_gens_para,
            generic_avrs_para)

        gen_vh_id_iq_V_ref_Tm_para = [gen_vh;
               gen_i_d;
               gen_i_q;
               gen_V_ref;
               gen_Tm] 
        
        a_kwd_para = (;gen_para,
                      avr_para,
                      ωs)

        dae_a_gen_generic_model_by_ext_idq_func!(
            res[state_var_idx],
            dx[state_var_idx],
            x[state_var_idx],
            gen_vh_id_iq_V_ref_Tm_para,
            t;
            kwd_para =
                a_kwd_para  )

    end

    return nothing    

end



#-----------------------------------------------------
# generic model system dynamicss
#-----------------------------------------------------

function dae_generic_model_dynamics!(
    res,
    dx,
    x,
    dae_SC_generic_model_dynamics_para,
    t;
    kwd_para =
        dae_SC_generic_model_dynamics_kwd_para  )

    #----------------------------------------

    (;
     ωs,
     loc_load_exist,
     state_vars_idx ,
     SC_generic_model_states_idx_in_state_Idx ,

     dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
     SC_generic_model_vars_wt_i_dq_Idx_in_state,
     SC_generic_model_states_comp_idxs_in_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,

     pf_gens_para,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     generic_5_gens_type_paras_Idx,
     generic_gens_para,
     generic_avrs_para,
     
     gens_ode_kwd_para) =
         NamedTupleTools.select(
             kwd_para,
             (:ωs,
              :loc_load_exist,
              :state_vars_idx ,
              :SC_generic_model_states_idx_in_state_Idx ,

              :dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
              :SC_generic_model_vars_wt_i_dq_Idx_in_state,
              :SC_generic_model_states_comp_idxs_in_Idx,
              
              :dyn_pf_fun_kwd_n2s_idxs,
              :dyn_pf_fun_kwd_net_idxs,

              :pf_gens_para,
              :Ynet_wt_nodes_idx_wt_adjacent_nodes,

              :generic_5_gens_type_paras_Idx,
              :generic_gens_para,
              :generic_avrs_para,

              :gens_ode_kwd_para) )

    #----------------------------------------
    
    # needed
    (;dyn_V_ref_Idx,
     dyn_Tm_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_V_ref_Tm_Png_Qng_Pll_Qll_Idx,
             (:dyn_V_ref_Idx,
              :dyn_Tm_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             SC_generic_model_vars_wt_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    #----------------------------------------

    
    (δ_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             SC_generic_model_states_comp_idxs_in_Idx,
             (:δ_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state
              ))
    
    # -------------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))  

    #----------------------------------------
    
    (ra,
     X_d,
     X_q,
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
            pf_gens_para,
             (:ra,
              :X_d,
              :X_q,
              :X_d_dash,
              :X_q_dash ) )
    
    
    #----------------------------------------

   (Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        NamedTupleTools.select(
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            (:Ynet,
             :nodes_idx_with_adjacent_nodes_idx)) 
    
    #----------------------------------------
    #----------------------------------------

    state_vars = x[state_var_Idx_in_state]
    
    vh  = x[vh_Idx_in_state]
    
    θh  = x[θh_Idx_in_state]
    
    gens_i_d = x[id_Idx_in_state]
    
    gens_i_q = x[iq_Idx_in_state]
    
    #----------------------------------------
    
    gens_δ = x[δ_idx_in_state]
        
    gens_ed_dash = x[ed_dash_idx_in_state]
    
    gens_eq_dash = x[eq_dash_idx_in_state]
        
    #----------------------------------------
    
    gens_vh =
        vh[gens_nodes_idx]

    gens_θh =
        θh[gens_nodes_idx]

    non_gens_nodes_vh =
        vh[non_gens_nodes_idx]

    non_gens_nodes_θh =
        θh[non_gens_nodes_idx]
    
    #----------------------------------------
    
    # vh = [
    #     idx ∈ gens_nodes_idx ?
    #         gens_vh[
    #             n2s_gens_idx[ idx ] ]  :
    #                 non_gens_nodes_vh[
    #                     n2s_non_gens_idx[ idx ] ]
    #            for idx in all_nodes_idx ]

    # θh = [
    #     idx ∈ gens_nodes_idx ?
    #         gens_θh[
    #             n2s_gens_idx[ idx] ] :
    #                 non_gens_nodes_θh[
    #                     n2s_non_gens_idx[ idx ] ]
    #     for idx in all_nodes_idx ]
    
    #----------------------------------------
        
    V_ref = dae_SC_generic_model_dynamics_para[
        dyn_V_ref_Idx]
    
    Tm = dae_SC_generic_model_dynamics_para[
        dyn_Tm_Idx]
    
    P_non_gens = dae_SC_generic_model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = dae_SC_generic_model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            dae_SC_generic_model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            dae_SC_generic_model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------    
    # single 
    #----------------------------------------

    for (state_var_idx,
         a_gen_vh,
         a_gen_i_d,
         a_gen_i_q,
         a_v_ref,
         a_τm,
         a_gen_para,
         a_avr_para) in
        zip(state_vars_idx,
            gens_vh,
            gens_i_d,
            gens_i_q,
            V_ref,
            Tm,
            generic_gens_para,
            generic_avrs_para)

        gen_ode_kwd_para =
            (a_gen_para,
             a_avr_para,
             ωs)
        
        # gen_vh_id_iq_V_ref_Tm_para =
        #     [a_gen_vh, a_gen_i_d, a_gen_i_q,
        #      a_v_ref, a_τm ]
        
        dae_a_gen_generic_model_by_ext_idq_func!(
            res[state_var_idx],
            dx[state_var_idx],
            x[state_var_idx],
            [a_gen_vh, a_gen_i_d,
             a_gen_i_q, a_v_ref, a_τm ],
            t;
            kwd_para =
                gen_ode_kwd_para )
    end
    
    # #----------------------------------------
    # # ode parameters and kwd para
    # #----------------------------------------
    
    # gens_ode_kwd_para = [
    #     (gen_para =
    #         (gen_H, gen_D,
    #          gen_X_d, gen_X_q,
    #          gen_X_d_dash, gen_X_q_dash,
    #          gen_T_d_dash, gen_T_q_dash ),
    #      avr_para =
    #          (avr_Ka, avr_Ta,
    #           avr_Ke, avr_Te,
    #           avr_Kf, avr_Tf,
    #           avr_Ae, avr_Be,
    #           avr_Tr),
    #      ωs = ωs )
    #     for (gen_H, gen_D,
    #          gen_X_d, gen_X_q,
    #          gen_X_d_dash, gen_X_q_dash,
    #          gen_T_d_dash, gen_T_q_dash,
    #          avr_Ka, avr_Ta,
    #          avr_Ke, avr_Te,
    #          avr_Kf, avr_Tf,
    #          avr_Ae, avr_Be,
    #          avr_Tr) in
    #         zip(H, D,
    #          X_d, X_q,
    #          X_d_dash, X_q_dash,
    #          T_d_dash, T_q_dash,
    #          Ka, Ta,
    #          Ke, Te,
    #          Kf, Tf,
    #          Ae, Be,
    #          Tr )]    
    
    # gens_vh_id_iq_V_ref_Tm_para = [
    #     [a_vh, a_id, a_iq, a_vref, a_tm]
    #     for (a_vh, a_id, a_iq, a_vref, a_tm) in
    #         zip(gens_vh,
    #         gens_i_d,
    #         gens_i_q,
    #         V_ref,
    #             Tm) ]
    
    # #----------------------------------------        
    # #----------------------------------------
    # # ode 
    # #----------------------------------------
    # #----------------------------------------
    
    # ode_res_views = [view(res, idx)
    #                  for idx in state_vars_idx ]
    
    # ode_dx_views = [view(dx, idx)
    #                 for idx in state_vars_idx ]
    
    # ode_x_views  = [view( x, idx)
    #                 for idx in state_vars_idx ]
    
    # for (ode_res, ode_dx, ode_x,
    #      gen_vh_id_iq_V_ref_Tm_para,
    #      gen_ode_kwd_para) in
    #     zip(ode_res_views,
    #         ode_dx_views,
    #         ode_x_views,
    #         gens_vh_id_iq_V_ref_Tm_para,
    #         gens_ode_kwd_para)
        
    #     dae_a_gen_generic_model_by_ext_idq_func!(
    #         ode_res,
    #         ode_dx,
    #         ode_x,
    #         gen_vh_id_iq_V_ref_Tm_para,
    #         t;
    #         kwd_para =
    #             gen_ode_kwd_para )
        
    # end

    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------


    # -----------------------------------
    # active power mismatch
    # -----------------------------------
    
    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """

    # # P_mismatch
    
    res[vh_Idx_in_state] .= [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          vh[ n2s_all_nodes_idx[ nth_idx]] * gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -
          P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]] -
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])])) :
        (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                             nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx] ]) ]))
                        for nth_idx in all_nodes_idx ] 


    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------

    
    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    # # Q_mismatch

    res[θh_Idx_in_state] .= [ nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]] -
           vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
           (vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            vh[ n2s_all_nodes_idx[nth_idx]] * sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                              nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])] ) ) 
                              for nth_idx in all_nodes_idx ]

    
    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """
    
    res[id_Idx_in_state] .= [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin( gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """
    
    res[iq_Idx_in_state] .= [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[ n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  
    
    # ------------------------------------
    
    return nothing

end



#---------------------------------------------------
#---------------------------------------------------
# Dynamic components initialisation and modelling
#---------------------------------------------------
#---------------------------------------------------


function ode_sauer_gen_fun_SM_2axis_cb_v6_init(
    gen_vh,
    gen_θh,
    P_g,
    Q_g,
    ωs;
    kwd_para = a_sauer_plant_init_kwd_para  )
    
    (gen_para,
     avr_para) =
         kwd_para
    
    
    (ra, X_d, X_q, X_d_dash, X_q_dash) = # gen_para
        NamedTupleTools.select(
            gen_para,
            (:ra,
             :X_d,
             :X_q,
             :X_d_dash,
             :X_q_dash))

    
    (Ka, Ta,
     Ke, Te,
     Kf, Tf,
     Ae, Be,
     Tr ) = # avr_para
         NamedTupleTools.select(
             avr_para,
             (:Ka, :Ta,
              :Ke, :Te,
              :Kf, :Tf,
              :Ae, :Be,
              :Tr))
    
    Ig = (P_g - im * Q_g) / (
        gen_vh * exp(-im * gen_θh))

    E_gen = gen_vh * exp(im * gen_θh) +
        (ra + im * X_q) * Ig

    δ = angle( E_gen )

    E = abs(E_gen)

    i_dq =  abs(Ig) * exp(im * (angle(Ig) - δ + π/2) )

    v_dq =  gen_vh * exp(im * (gen_θh - δ + π/2) )

    i_d = real(i_dq )

    i_q = imag(i_dq )

    v_d = real(v_dq )

    v_q = imag(v_dq )

    #   ed_dash = (X_q - X_q_dash) * i_q
    
    ed_dash = v_d  + ra * i_d - X_q_dash * i_q
    
    eq_dash = v_q  + ra * i_q + X_d_dash * i_d

    E_fd = vf_tilade = eq_dash + (X_d - X_d_dash) * i_d

    vr1 = V_R = (Ke + Sevf(Ae, Be,  E_fd)) *  E_fd

    vr2 = R_f = Kf/Tf * E_fd

    v_ref = V_ref =  gen_vh + (V_R / Ka)

    Tm = τm_tilade = ed_dash * i_d + eq_dash * i_q +
        (X_q_dash - X_d_dash) * i_d * i_q

    δ_deg = δ * 180.0 /pi
    
    return  (states_vars =
        (δ, ωs, ed_dash, eq_dash,
         vr1, vr2, vf_tilade),
             ref_vars = (; v_ref, τm_tilade),
             others = (;E_fd, V_R, R_f, V_ref, Tm, E,
                       i_d, i_q),
             i_dq = (;i_d, i_q),
             v_dq = (;v_d, v_q),
             diagnostics =
                 (; δ,
                  δ_deg,
                  ed_dash,
                  eq_dash,
                  vr1,
                  vr2,
                  vf_tilade,
                  v_ref,
                  τm_tilade,
                  i_d,
                  i_q,
                  v_d,
                  v_q) ) 

end


function ode_sauer_gens_fun_SM_2axis_cb_v6_init(
    gens_vh,
    gens_θh,
    pf_P_gens,
    pf_Q_gens,
    ωs;
    kwd_para =
        sauer_plants_init_kwd_para )

    (gens_para,
     avrs_para ) =
         kwd_para

    states_init_wt_ref =
        [ ode_sauer_gen_fun_SM_2axis_cb_v6_init(
            gen_vh,
            gen_θh,
            pf_P_gen,
            pf_Q_gen,
            ωs;
            kwd_para =
                (gen_para,
                 avr_para))
          for (gen_vh,
               gen_θh,
               pf_P_gen,
               pf_Q_gen,
               gen_para,
               avr_para ) in
              zip( gens_vh,
                   gens_θh,
                   pf_P_gens,
                   pf_Q_gens,
                   gens_para,
                   avrs_para) ]

    #----------------------------------------
    #----------------------------------------
    
    gens_plants_states_init =
        [[[ a_state for a_state in
               a_plant_state.states_vars]
          for a_plant_state in
              states_init_wt_ref ]...;]
    
    gens_plants_refs = [a_plant_state.ref_vars
         for a_plant_state in
                   states_init_wt_ref ]

    gens_plants_refs =
        NamedTupleTools.select(
            get_selected_vec_nt_to_vec_vec(
                gens_plants_refs,
                nothing;
                selections = (:v_ref, :τm_tilade ),
                vec_datatype = Float64 ),
            (:v_ref, :τm_tilade ) )
    #
    
    gens_plants_i_dq =
        [a_plant_state.i_dq
         for a_plant_state in
                   states_init_wt_ref ]
    
    gens_plants_i_dq =
        NamedTupleTools.select(
            get_selected_vec_nt_to_vec_vec(
                gens_plants_i_dq,
                nothing;
                selections = ( :i_d, :i_q ),
                vec_datatype = Float64 ),
            ( :i_d, :i_q ) )
    
    
    gens_plants_states_init_other_cal =
        [ a_plant_state.others
          for a_plant_state in
              states_init_wt_ref ]

    
    gens_plants_states_init_diagnostics =
        [ a_plant_state.diagnostics
          for a_plant_state in
              states_init_wt_ref ]
    
    
    
    #----------------------------------------

    return (;gens_plants_states_init,
            gens_plants_refs,
            gens_plants_i_dq,
            gens_plants_states_init_other_cal,
            states_init_wt_ref,
            gens_plants_states_init_diagnostics)
    
end


#---------------------------------------------
# ODE
#---------------------------------------------


function ode_sauer_gen_fun_SM_2axis_cb_v6_func!(
    dx, x,
    gen_vh_id_iq_v_ref0_porder0_para,
    t;
    kwd_para =
        gen_avr_kwd_para )
    
    (vh,
     i_d,
     i_q,
     v_ref0,
     porder0) =
        gen_vh_id_iq_v_ref0_porder0_para
    
    (gen_para,
     avr_para,
     ωs) =
         kwd_para
    
    (;H,
     D,
     
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gen_para,
             (:H,
              :D,
              
              :X_d,
              :X_q,
              
              :X_d_dash,
              :X_q_dash,
              
              :T_d_dash,
              :T_q_dash))

    
    (Ka, Ta,
     Ke, Te,
     Kf, Tf,
     Ae, Be,
     Tr ) =
         NamedTupleTools.select(
             avr_para,
             (:Ka, :Ta,
              :Ke, :Te,
              :Kf, :Tf,
              :Ae, :Be,
              :Tr))


    
    τm = porder0

    δ, ω, ed_dash, eq_dash, vr1, vr2, vf_tilade = x

    # Gen
    
    dx[1] = ω - ωs 
    
    # dx[2] = ( ωs /(2*H) ) * ( τm -
    #     ed_dash * i_d - eq_dash * i_q -
    #     (X_q_dash - X_d_dash) * i_d * i_q - D * (ω - ωs)) 

    
    dx[2] = ( ωs /(2*H) ) * ( τm -
        ed_dash * i_d - eq_dash * i_q -
        (X_q_dash - X_d_dash) * i_d * i_q) 
    
    dx[3] = (1/T_q_dash) *
        ( (X_q - X_q_dash) * i_q  - ed_dash ) 

    dx[4] = (1/T_d_dash) *
        (vf_tilade - (X_d - X_d_dash) * i_d  -
        eq_dash) 

    # AVR

    dx[5] = (1/Ta) * ( Ka * vr2 + Ka * (v_ref0 - vh) -
        (Ka * Kf/Tf) * vf_tilade  - vr1 )    

    dx[6] = (1/Tf) * ((Kf/Tf) * vf_tilade - vr2 )

    dx[7] = (1/Te) * ( vr1 - (Ke + Sevf(Ae, Be, vf_tilade)) * vf_tilade)     
    
    return nothing    

end



function ode_sauer_gens_fun_SM_2axis_cb_v6_func!(
    dx, x,
    vh_id_iq_vref0_porder0,
    t;
    kwd_para =
        ode_sauer_gens_plants_kwd_para )

    (state_vars_idx,

     # vh_id_iq_vref0_porder0_Idx,
     # model_dynamics_kwd_para,
     dyn_5_gens_paras_Idx,

     generic_gens_para,
     generic_avrs_para,

     avrs_cb_sw,

     ωs) =
        kwd_para

    # (vh_Idx,
    #  id_Idx,    
    #  iq_Idx,
    #  vref0_Idx,    
    #  porder0_Idx
    #  ) = vh_id_iq_vref0_porder0_Idx


    # vh = vh_id_iq_vref0_porder0[
    #     vh_Idx]
    
    # id =
    #     vh_id_iq_vref0_porder0[
    #         id_Idx]
    
    # iq =
    #     vh_id_iq_vref0_porder0[
    #         iq_Idx]
    
    # vref0 = vh_id_iq_vref0_porder0[
    #     vref0_Idx]
    
    # porder0 =
    #     vh_id_iq_vref0_porder0[
    #         porder0_Idx]
    
    (;
     gen_para_1_Idx,
     gen_para_2_Idx,
     gen_para_3_Idx,
     gen_para_4_Idx,
     gen_para_5_Idx) =
        dyn_5_gens_paras_Idx
    

    vh = vh_id_iq_vref0_porder0[
        gen_para_1_Idx]
    
    id =
        vh_id_iq_vref0_porder0[
            gen_para_2_Idx]
    
    iq =
        vh_id_iq_vref0_porder0[
            gen_para_3_Idx]
    
    vref0 = vh_id_iq_vref0_porder0[
        gen_para_4_Idx]
    
    porder0 =
        vh_id_iq_vref0_porder0[
            gen_para_5_Idx]
    
    #----------------------------------------
    
    for (state_var_idx,

         a_vh,
         a_id,
         a_iq,
         a_vref0,
         a_porder0,
             
         gen_para,
         avr_para,

         avr_cb_sw) in
        zip( state_vars_idx,

             vh,
             id,
             iq,
             vref0,
             porder0,
             
             generic_gens_para,
             generic_avrs_para,

             avrs_cb_sw )

            gen_vh_id_iq_v_ref0_porder0_para =
                [a_vh,
                 a_id,
                 a_iq,
                 a_vref0,
                 a_porder0
                 ]
        
        gen_avr_kwd_para =
            (gen_para,
             avr_para,
             ωs)

        ode_sauer_gen_fun_SM_2axis_cb_v6_func!(
            dx[state_var_idx],
            x[state_var_idx],
            gen_vh_id_iq_v_ref0_porder0_para,
            t;
            kwd_para =
                gen_avr_kwd_para )

    end
    
    #----------------------------------------

    return nothing
end



function ode_sauer_system_model_dynamics!(
    dx,
    x,
    model_dynamics_para,
    t;
    kwd_para =
        model_dynamics_kwd_para  )

    #----------------------------------------

    (;
     state_vars_idx,
     state_vars_and_i_dq_Idx_in_state,
     dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     gens_state_vars_idx_in_state,

     dyn_5_gens_paras_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     non_pre_ordered_pf_vars_Idxs,

    # vh_id_iq_vref0_porder0_Idx,

     pf_gens_para,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     ode_sauer_gens_plants_kwd_para,
     state_algebraic_vars_Idx_in_state,
     algebraic_generic_model_sol_kwd_para) =
        NamedTupleTools.select(
            kwd_para,
            (:state_vars_idx,
             :state_vars_and_i_dq_Idx_in_state,          
             :dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx ,
             :gens_state_vars_idx_in_state ,

             :dyn_5_gens_paras_Idx,

             :dyn_pf_fun_kwd_n2s_idxs,
             :dyn_pf_fun_kwd_net_idxs,
             :non_pre_ordered_pf_vars_Idxs,

             # :vh_id_iq_vref0_porder0_Idx,

             :pf_gens_para,
             :Ynet_wt_nodes_idx_wt_adjacent_nodes,

             :ode_sauer_gens_plants_kwd_para,
             :state_algebraic_vars_Idx_in_state,
             :algebraic_generic_model_sol_kwd_para))

    #----------------------------------------
    
    (;
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    
  
    (;state_var_Idx_in_state,
     algebraic_var_Idx_in_state ) =
         NamedTupleTools.select(
             state_algebraic_vars_Idx_in_state,
             (:state_var_Idx_in_state,
              :algebraic_var_Idx_in_state))
    
    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    #----------------------------------------

    
    (;δ_idx_in_state,
     ω_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ω_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    # -------------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))  

    #----------------------------------------
    
    (;ra,
     X_d,
     X_q,     
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
             pf_gens_para,
             (:ra,
              :X_d,
              :X_q,     
              :X_d_dash,
              :X_q_dash ))
    
    #----------------------------------------

   (Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        NamedTupleTools.select(
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            (:Ynet,
             :nodes_idx_with_adjacent_nodes_idx)) 
    
    #----------------------------------------

    state_vars = x[state_var_Idx_in_state]
    
    vh  = x[vh_Idx_in_state]
    
    θh  = x[θh_Idx_in_state]
    
    gens_i_d = x[id_Idx_in_state]
    
    gens_i_q = x[iq_Idx_in_state]
    
    #----------------------------------------
    
    gens_δ = x[δ_idx_in_state]
        
    gens_ed_dash = x[ed_dash_idx_in_state]
    
    gens_eq_dash = x[eq_dash_idx_in_state]
        
    #----------------------------------------
    
    gens_vh =
        vh[
            gens_nodes_idx ]

    gens_θh =
        θh[
            gens_nodes_idx ]

    non_gens_nodes_vh =
        flat_vh[
            non_gens_nodes_idx ]

    non_gens_nodes_θh =
        flat_θh[
            non_gens_nodes_idx ]
            
    #----------------------------------------
    
    v_ref = model_dynamics_para[
        dyn_v_ref_Idx]
    
    p_order = model_dynamics_para[
        dyn_p_order_Idx]
    
    P_non_gens = model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------

     # ode_res_views = @view res[state_var_Idx_in_state]

     # ode_dx_views  = @view  dx[state_var_Idx_in_state]

     # ode_x_views   =  @view x[state_var_Idx_in_state]
    
    #----------------------------------------
    # ode parameters and kwd para
    #----------------------------------------
    
    # gens_vh_id_iq_vref0_porder0_para_per_gen = [
    #     [a_vh, a_id, a_iq, a_vref, a_porder0 ]
    #     for (a_vh,
    #          a_id,
    #          a_iq,
    #          a_vref,
    #          a_porder0
    #          ) in
    #         zip(gens_vh, 
    #             gens_i_d,
    #             gens_i_q,
    #             v_ref,
    #             p_order
    #              ) ]

    #----------------------------------------        

    vh_id_iq_vref0_porder0 =
        [gens_vh; 
         gens_i_d;
         gens_i_q;
         v_ref;
         p_order]
    
    #----------------------------------------        
    #----------------------------------------
    # ode 
    #----------------------------------------
    #----------------------------------------
    
    ode_sauer_gens_fun_SM_2axis_cb_v6_func!(
        dx[state_var_Idx_in_state],
        x[state_var_Idx_in_state],
        vh_id_iq_vref0_porder0,
        t;
        kwd_para =
            ode_sauer_gens_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para =
        [gens_δ;
         gens_ed_dash;
         gens_eq_dash;
         P_non_gens;
         Q_non_gens;
         P_g_loc_load;
         Q_g_loc_load]
    
    pf_algebraic_generic_model_sol(
        dx[algebraic_var_Idx_in_state],
        x[algebraic_var_Idx_in_state],
        δ_ed_dash_eq_dash_Png_Qng_Pll_Qll_para;
        kwd_para =
            algebraic_generic_model_sol_kwd_para )
    
    # ------------------------------------
    
    return nothing

end



#---------------------------------------------
# DAE
#---------------------------------------------


function dae_sauer_gen_fun_SM_2axis_cb_v6_func!(
    res,
    dx, x,
    gen_vh_id_iq_v_ref0_porder0_para,
    t;
    kwd_para =
        gen_avr_kwd_para )
    
    (vh,
     i_d,
     i_q,
     v_ref0,
     porder0) =
        gen_vh_id_iq_v_ref0_porder0_para
    
    (gen_para,
     avr_para,
     ωs) =
         kwd_para
    
    (;H,
     D,
     
     X_d,
     X_q,
     
     X_d_dash,
     X_q_dash,
     
     T_d_dash,
     T_q_dash ) =
         NamedTupleTools.select(
             gen_para,
             (:H,
              :D,
              
              :X_d,
              :X_q,
              
              :X_d_dash,
              :X_q_dash,
              
              :T_d_dash,
              :T_q_dash))

    
    (Ka, Ta,
     Ke, Te,
     Kf, Tf,
     Ae, Be,
     Tr ) =
         NamedTupleTools.select(
             avr_para,
             (:Ka, :Ta,
              :Ke, :Te,
              :Kf, :Tf,
              :Ae, :Be,
              :Tr))


    
    τm = porder0

    δ, ω, ed_dash, eq_dash, vr1, vr2, vf_tilade = x

    # Gen
    
    # if abs( dx[ 2 ] ) < Δω_toleraance
        
    #     Δω = ωs - ω
        
    #     res[1] = ω - ωs - dx[1] +  Δω

    # else
        
    #     res[1] = ω - ωs - dx[1]
        
    # end

    res[1] = ω - ωs - dx[1]
    
    # res[2] = ( ωs /(2*H) ) * ( τm -
    #     ed_dash * i_d - eq_dash * i_q -
    #     (X_q_dash - X_d_dash) * i_d * i_q - D * (ω - ωs)) -
    #     dx[2]

    
    res[2] = ( ωs /(2*H) ) * ( τm -
        ed_dash * i_d - eq_dash * i_q -
        (X_q_dash - X_d_dash) * i_d * i_q) -
        dx[2]
    
    res[3] = (1/T_q_dash) *
        ( (X_q - X_q_dash) * i_q  - ed_dash) -
        dx[3]

    res[4] = (1/T_d_dash) *
        (vf_tilade - (X_d - X_d_dash) * i_d  -
        eq_dash) -
        dx[4]

    # AVR

    res[5] = (1/Ta) * ( Ka * vr2 + Ka * (v_ref0 - vh) -
        (Ka * Kf/Tf) * vf_tilade  - vr1 ) - dx[5]    

    res[6] = (1/Tf) * ((Kf/Tf) * vf_tilade - vr2 ) -
        dx[6]

    res[7] = (1/Te) * ( vr1 - (Ke + Sevf(Ae, Be, vf_tilade)) * vf_tilade) -
        dx[7]    
    
    return nothing    

end



function dae_sauer_gens_fun_SM_2axis_cb_v6_func!(
    res,
    dx, x,
    vh_id_iq_vref0_porder0,
    t;
    kwd_para =
        ode_sauer_gens_plants_kwd_para )

    (state_vars_idx,

     # vh_id_iq_vref0_porder0_Idx,
     # model_dynamics_kwd_para,
     dyn_5_gens_paras_Idx,

     generic_gens_para,
     generic_avrs_para,

     avrs_cb_sw,

     ωs) =
        kwd_para

    # (vh_Idx,
    #  id_Idx,    
    #  iq_Idx,
    #  vref0_Idx,    
    #  porder0_Idx
    #  ) = vh_id_iq_vref0_porder0_Idx


    # vh = vh_id_iq_vref0_porder0[
    #     vh_Idx]
    
    # id =
    #     vh_id_iq_vref0_porder0[
    #         id_Idx]
    
    # iq =
    #     vh_id_iq_vref0_porder0[
    #         iq_Idx]
    
    # vref0 = vh_id_iq_vref0_porder0[
    #     vref0_Idx]
    
    # porder0 =
    #     vh_id_iq_vref0_porder0[
    #         porder0_Idx]
    
    (;
     gen_para_1_Idx,
     gen_para_2_Idx,
     gen_para_3_Idx,
     gen_para_4_Idx,
     gen_para_5_Idx) =
        dyn_5_gens_paras_Idx
    

    vh = vh_id_iq_vref0_porder0[
        gen_para_1_Idx]
    
    id =
        vh_id_iq_vref0_porder0[
            gen_para_2_Idx]
    
    iq =
        vh_id_iq_vref0_porder0[
            gen_para_3_Idx]
    
    vref0 = vh_id_iq_vref0_porder0[
        gen_para_4_Idx]
    
    porder0 =
        vh_id_iq_vref0_porder0[
            gen_para_5_Idx]
    
    #----------------------------------------
    
    for (state_var_idx,

         a_vh,
         a_id,
         a_iq,
         a_vref0,
         a_porder0,
             
         gen_para,
         avr_para,

         avr_cb_sw) in
        zip( state_vars_idx,

             vh,
             id,
             iq,
             vref0,
             porder0,
             
             generic_gens_para,
             generic_avrs_para,

             avrs_cb_sw )

            gen_vh_id_iq_v_ref0_porder0_para =
                [a_vh,
                 a_id,
                 a_iq,
                 a_vref0,
                 a_porder0
                 ]
        
        gen_avr_kwd_para =
            (gen_para,
             avr_para,
             ωs)

        dae_sauer_gen_fun_SM_2axis_cb_v6_func!(
            res[state_var_idx],
            dx[state_var_idx],
            x[state_var_idx],
            gen_vh_id_iq_v_ref0_porder0_para,
            t;
            kwd_para =
                gen_avr_kwd_para )

    end
    
    #----------------------------------------

    return nothing
end



function dae_sauer_system_model_dynamics!(
    res,
    dx,
    x,
    model_dynamics_para,
    t;
    kwd_para =
        model_dynamics_kwd_para  )

    #----------------------------------------

    (;
     state_vars_idx,
     state_vars_and_i_dq_Idx_in_state,
     dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
     gens_state_vars_idx_in_state,

     dyn_5_gens_paras_Idx,

     dyn_pf_fun_kwd_n2s_idxs,
     dyn_pf_fun_kwd_net_idxs,
     non_pre_ordered_pf_vars_Idxs,

    # vh_id_iq_vref0_porder0_Idx,

     pf_gens_para,
     Ynet_wt_nodes_idx_wt_adjacent_nodes,

     ode_sauer_gens_plants_kwd_para) =
        NamedTupleTools.select(
            kwd_para,
            (:state_vars_idx,
             :state_vars_and_i_dq_Idx_in_state,          
             :dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx ,
             :gens_state_vars_idx_in_state ,

             :dyn_5_gens_paras_Idx,

             :dyn_pf_fun_kwd_n2s_idxs,
             :dyn_pf_fun_kwd_net_idxs,
             :non_pre_ordered_pf_vars_Idxs,

             # :vh_id_iq_vref0_porder0_Idx,

             :pf_gens_para,
             :Ynet_wt_nodes_idx_wt_adjacent_nodes,

             :ode_sauer_gens_plants_kwd_para ))

    #----------------------------------------
    
    (;
     dyn_v_ref_Idx,
     dyn_p_order_Idx,
     dyn_Png_Idx,
     dyn_Qng_Idx,
     dyn_Pll_Idx,
     dyn_Qll_Idx ) =
         NamedTupleTools.select(
             dyn_v_ref_p_order_Png_Qng_Pll_Qll_Idx,
             (
              :dyn_v_ref_Idx,
              :dyn_p_order_Idx,
              :dyn_Png_Idx,
              :dyn_Qng_Idx,
              :dyn_Pll_Idx,
              :dyn_Qll_Idx))

    #----------------------------------------
    
        (;state_var_Idx_in_state,
         vh_Idx_in_state,
         θh_Idx_in_state,
         id_Idx_in_state,
         iq_Idx_in_state) = NamedTupleTools.select(
             state_vars_and_i_dq_Idx_in_state,
                 (:state_var_Idx_in_state,
                  :vh_Idx_in_state,
                  :θh_Idx_in_state,
                  :id_Idx_in_state,
                  :iq_Idx_in_state))
    
    #----------------------------------------

    
    (;δ_idx_in_state,
     ω_idx_in_state,
     ed_dash_idx_in_state,
     eq_dash_idx_in_state ) =
         NamedTupleTools.select(
             gens_state_vars_idx_in_state,
             (:δ_idx_in_state,
              :ω_idx_in_state,
              :ed_dash_idx_in_state,
              :eq_dash_idx_in_state))
    
    # -------------------------------------

   (;
    n2s_slack_gens_idx,
    n2s_non_slack_gens_idx,
    n2s_gens_idx,
    n2s_non_gens_idx,
    n2s_gens_with_loc_load_idxs,
    n2s_all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_n2s_idxs,
            (:n2s_slack_gens_idx,
             :n2s_non_slack_gens_idx,
             :n2s_gens_idx,
             :n2s_non_gens_idx,
             :n2s_gens_with_loc_load_idxs,
             :n2s_all_nodes_idx))
    
    #----------------------------------------

   (;
    slack_gens_nodes_idx,
    non_slack_gens_nodes_idx,
    gens_nodes_idx,
    non_gens_nodes_idx,
    gens_nodes_with_loc_loads_idx,
    all_nodes_idx ) =
        NamedTupleTools.select(
            dyn_pf_fun_kwd_net_idxs,
            (:slack_gens_nodes_idx,
             :non_slack_gens_nodes_idx,
             :gens_nodes_idx,
             :non_gens_nodes_idx,
             :gens_nodes_with_loc_loads_idx,
             :all_nodes_idx))  

    #----------------------------------------
    
    (;ra,
     X_d,
     X_q,     
     X_d_dash,
     X_q_dash ) =
         NamedTupleTools.select(
             pf_gens_para,
             (:ra,
              :X_d,
              :X_q,     
              :X_d_dash,
              :X_q_dash ))
    
    #----------------------------------------

   (Ynet,
    nodes_idx_with_adjacent_nodes_idx) =
        NamedTupleTools.select(
            Ynet_wt_nodes_idx_wt_adjacent_nodes,
            (:Ynet,
             :nodes_idx_with_adjacent_nodes_idx)) 
    
    #----------------------------------------

    state_vars = x[state_var_Idx_in_state]
    
    flat_vh  = x[vh_Idx_in_state]
    
    flat_θh  = x[θh_Idx_in_state]
    
    gens_i_d = x[id_Idx_in_state]
    
    gens_i_q = x[iq_Idx_in_state]
    
    #----------------------------------------
    
    gens_δ = x[δ_idx_in_state]
        
    gens_ed_dash = x[ed_dash_idx_in_state]
    
    gens_eq_dash = x[eq_dash_idx_in_state]
        
    #----------------------------------------
    
    gens_vh =
        flat_vh[
            gens_nodes_idx ]

    gens_θh =
        flat_θh[
            gens_nodes_idx ]

    non_gens_nodes_vh =
        flat_vh[
            non_gens_nodes_idx ]

    non_gens_nodes_θh =
        flat_θh[
            non_gens_nodes_idx ]
            
    #----------------------------------------
    
    vh = [
        idx ∈ gens_nodes_idx ?
            gens_vh[
                n2s_gens_idx[ idx ] ]  :
                    non_gens_nodes_vh[
                        n2s_non_gens_idx[ idx ] ]
               for idx in all_nodes_idx ]

    θh = [
        idx ∈ gens_nodes_idx ?
            gens_θh[
                n2s_gens_idx[ idx] ] :
                    non_gens_nodes_θh[
                        n2s_non_gens_idx[ idx ] ]
        for idx in all_nodes_idx ]
    
    #----------------------------------------
    
    v_ref = model_dynamics_para[
        dyn_v_ref_Idx]
    
    p_order = model_dynamics_para[
        dyn_p_order_Idx]
    
    P_non_gens = model_dynamics_para[
        dyn_Png_Idx]
    
    Q_non_gens = model_dynamics_para[
        dyn_Qng_Idx]
    
    if loc_load_exist == true
        
        P_g_loc_load =
            model_dynamics_para[
                dyn_Pll_Idx]
        
        Q_g_loc_load =
            model_dynamics_para[
                dyn_Qll_Idx]        
    else

        P_g_loc_load  = []
        
        Q_g_loc_load  = []
        
    end

    #----------------------------------------

     # ode_res_views = @view res[state_var_Idx_in_state]

     # ode_dx_views  = @view  dx[state_var_Idx_in_state]

     # ode_x_views   =  @view x[state_var_Idx_in_state]
    
    #----------------------------------------
    # ode parameters and kwd para
    #----------------------------------------
    
    # gens_vh_id_iq_vref0_porder0_para_per_gen = [
    #     [a_vh, a_id, a_iq, a_vref, a_porder0 ]
    #     for (a_vh,
    #          a_id,
    #          a_iq,
    #          a_vref,
    #          a_porder0
    #          ) in
    #         zip(gens_vh, 
    #             gens_i_d,
    #             gens_i_q,
    #             v_ref,
    #             p_order
    #              ) ]

    #----------------------------------------        

    vh_id_iq_vref0_porder0 =
        [gens_vh; 
         gens_i_d;
         gens_i_q;
         v_ref;
         p_order]
    
    #----------------------------------------        
    #----------------------------------------
    # ode 
    #----------------------------------------
    #----------------------------------------
    
    dae_sauer_gens_fun_SM_2axis_cb_v6_func!(
        res[state_var_Idx_in_state],
        dx[state_var_Idx_in_state],
        x[state_var_Idx_in_state],
        vh_id_iq_vref0_porder0,
        t;
        kwd_para =
            ode_sauer_gens_plants_kwd_para )
    
    #----------------------------------------    
    #----------------------------------------
    # power flow equations
    #----------------------------------------
    #----------------------------------------

    # -----------------------------------
    # active power mismatch
    # -----------------------------------
    
    """
    gens:

    id * vh_i * sin(δ_i - θ_i) + iq * vh_i * cos(δ_i - θ_i) - Pl_i =
            vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    loads:

    Pl_i =
         -vh_i * ∑(vh_k * Y_ik * cos(θ_i - θ_k - β_ik ) )

    """
    
    res[vh_Idx_in_state] .= [ nth_idx ∈ non_gens_nodes_idx ?
         (P_non_gens[ n2s_non_gens_idx[ nth_idx]] +
          vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[ idx]] * abs(ynj) *
          cos( θh[ n2s_all_nodes_idx[nth_idx]] -  θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]]) ] ) ) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
         (vh[ n2s_all_nodes_idx[ nth_idx]] * (gens_i_d[ n2s_gens_idx[ nth_idx]] * 
          sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) +
          gens_i_q[ n2s_gens_idx[ nth_idx]] * 
          cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[ nth_idx]]) -          
          sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
                cos(θh[ n2s_all_nodes_idx[ nth_idx]] - θh[ n2s_all_nodes_idx[ idx]] - angle(ynj) )
                      for (ynj, idx) in
                          zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                              nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])]))  -
                                   P_g_loc_load[ n2s_gens_with_loc_load_idxs[nth_idx]]) :
        (vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
         sin( gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) +
         gens_i_q[ n2s_gens_idx[nth_idx]] * 
         cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] )  -
         sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
         cos(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[ nth_idx] ],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[ nth_idx]])]) ))
                        for nth_idx in all_nodes_idx ] 

    # -----------------------------------
    # reactive power mismatch
    # -----------------------------------
    
    """
    gens:

    id * vh_i * cos(δ_i - θ_i) - iq * vh_i * sin(δ_i - θ_i) - Ql_i =
            vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    loads:

    Ql_i =
         -vh_i * ∑(vh_k * Y_ik * sin(θ_i - θ_k - β_ik ) )

    """

    res[θh_Idx_in_state]  .= [ nth_idx ∈ non_gens_nodes_idx ?
        ( Q_non_gens[ n2s_non_gens_idx[ nth_idx ]] +
            vh[ n2s_all_nodes_idx[nth_idx]] * sum([ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                       zip( Ynet[ n2s_all_nodes_idx[nth_idx] ],
                            nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) :
                                nth_idx ∈ gens_nodes_with_loc_loads_idx ?
          ( vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
           cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
           gens_i_q[ n2s_gens_idx[nth_idx] ] * 
           sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]] ) -           
           sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
           sin( θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                   for (ynj, idx) in
                           zip(Ynet[ n2s_all_nodes_idx[nth_idx]],
                               nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])])) -
           Q_g_loc_load[n2s_gens_with_loc_load_idxs[nth_idx]]) :
           (vh[ n2s_all_nodes_idx[nth_idx]] * (gens_i_d[ n2s_gens_idx[nth_idx]] * 
            cos(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            gens_i_q[ n2s_gens_idx[nth_idx] ] * 
            sin(gens_δ[ n2s_gens_idx[nth_idx]] - θh[ n2s_all_nodes_idx[nth_idx]]) -
            sum( [ vh[ n2s_all_nodes_idx[idx]] * abs(ynj) *
            sin(θh[ n2s_all_nodes_idx[nth_idx]] - θh[ n2s_all_nodes_idx[idx]] - angle(ynj) )
                       for (ynj, idx) in
                           zip( Ynet[ n2s_all_nodes_idx[nth_idx]],
                                nodes_idx_with_adjacent_nodes_idx[ n2s_all_nodes_idx[nth_idx]])]))) 
                              for nth_idx in all_nodes_idx ]

    # -----------------------------------
    #  gens real part stator equation mismatch
    # -----------------------------------
    
    """ ed_dash - vh * sin(δ - θh) - ra * id +  Xq_dash * iq -   """
    
    res[id_Idx_in_state] .= [
        gens_ed_dash[ n2s_gens_idx[ idx ]] -
            gens_vh[ n2s_gens_idx[idx]] * sin( gens_δ[n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_d[ n2s_gens_idx[ idx]] + 
            X_q_dash[ n2s_gens_idx[idx]] * gens_i_q[n2s_gens_idx[idx]]
                    
        for idx in gens_nodes_idx ]  
    
    # -----------------------------------
    #  gens imag part stator equation mismatch
    # -----------------------------------

    """ eq_dash - vh * cos(δ - θh)  - ra * iq  - Xd_dash * id """
    
    res[iq_Idx_in_state] .= [
        gens_eq_dash[ n2s_gens_idx[ idx ]] - 
            gens_vh[ n2s_gens_idx[idx]] * cos(gens_δ[ n2s_gens_idx[idx]] -
            gens_θh[ n2s_gens_idx[idx]]) -
            ra[ n2s_gens_idx[ idx]] * gens_i_q[ n2s_gens_idx[ idx]] -
            X_d_dash[ n2s_gens_idx[idx]] * gens_i_d[ n2s_gens_idx[ idx]] 
        for idx in gens_nodes_idx ]  
    
    # ------------------------------------
    
    return nothing

end


#---------------------------------------------
#---------------------------------------------


# comment

