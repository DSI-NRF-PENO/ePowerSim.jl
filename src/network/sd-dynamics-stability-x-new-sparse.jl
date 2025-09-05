# (C) 2024 Power, Energy, Networks and Optimisation Research Group, Unisa, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

# AA Yusuff : yusufaa@unisa.ac.za


#########################################################
# ------------------------------------------------------
#  stability
# ------------------------------------------------------
#########################################################


"""

dx  = A × x + B × u

 y  = C × x + D × u


dx = A × x + B × u_network + C × u_setpoints

 y  = C × x + D × u

output y is incorporated into A as algebraic variables

Each controller's system matrix is formulated as


Ax_k_S = [
    Ax_cont_k_S_Sgen,
    Ax_cont_k_S_cont_1_S, ...,
    Ax_cont_k_S_cont_k_S, ...,
    Ax_cont_k_S_cont_m_S,
    
    Ax_cont_k_S_cont_1_A, ...,
    Ax_cont_k_S_cont_k_A, ...,
    Ax_cont_k_S_cont_m_A ] 

where Ax_cont_k_S_Sgen is the coupling of state variables of
controller k to the state variables of the generator it is embedded in.

Ax_cont_k_S_cont_m_S and Ax_cont_k_S_cont_m_A are the coupling of
the state variables of controler k to the state and algebraic
variables of controller m respectively.

Similarly, 

Ax_k_A = [
    Ax_cont_k_A_Sgen,
    Ax_cont_k_A_cont_1_S, ...,
    Ax_cont_k_A_cont_k_S, ...,
    Ax_cont_k_A_cont_m_S,
    
    Ax_cont_k_A_cont_1_A, ...,
    Ax_cont_k_A_cont_k_A, ...,
    Ax_cont_k_A_cont_m_A ] 


where Ax_cont_k_A_Sgen is the coupling of algbebraic variables of
controller k to the state variables of the generator it is embedded in.

Ax_cont_k_A_cont_m_S and  Ax_cont_k_A_cont_m_A are the coupling of
the algebraic variables of controler k to the state and algebraic
variables of controller m respectively.

The governing equation for a generator is also represented as

Ax_Sgen = [
    Ax_Sgen_cont_1_S,
    Ax_Sgen_cont_k_S, ...,
    Ax_Sgen_cont_m_S,
    
    Ax_Sgen_cont_1_A, ...,
    Ax_Sgen_cont_k_A, ...,
    Ax_Sgen_cont_m_A ] 

where Ax_Sgen_cont_m_S and Ax_Sgen_cont_m_A are the coupling of the
state variables of a generator to the state and algebraic variables of controller m.

The state variables of a generator are δ, ω, e_d_dash and e_q_dash.

The setpoints are  ωs, ω_ref0, v_ref0 and p_order0

The network injections are id,  iq,  ph,  vh


Implemented controllers for now:

gov_ieee_tgov1_cb
gov_t0_cb
gov_t1_cb
gov_t1_cb_sauer

avr_t0_cb
avr_t1_cb
avr_t1_cb_sauer

This can easily be extended to any controller.

"""


# ------------------------------------------------------
# ------------------------------------------------------
# utility functions
# ------------------------------------------------------
# ------------------------------------------------------


function zero_coupling_state_vars_to_state_vars(
    comp_A, comp_B )

    return zeros(
        length(comp_A.state_vars_syms),
        length(comp_B.state_vars_syms))
end


function zero_coupling_state_vars_to_im_alg_vars(
    comp_A, comp_B )

    return zeros(
        length(comp_A.state_vars_syms),
        length(comp_B.im_algebraic_vars_syms))
end


function zero_coupling_im_alg_vars_to_state_vars(
    comp_A, comp_B)

    return zeros(
        length(comp_A.im_algebraic_vars_syms),
        length(comp_B.state_vars_syms))
end


function zero_coupling_im_alg_vars_to_im_alg_vars(
    comp_A, comp_B)

    return zeros(
        length(comp_A.im_algebraic_vars_syms),
        length(comp_B.im_algebraic_vars_syms))
end


# ------------------------------------------------------
# ------------------------------------------------------
# control devices matrix
# ------------------------------------------------------
# ------------------------------------------------------


# ------------------------------------------------------
# GOV with_control_alg_vars = true
# ------------------------------------------------------


# ------------------------------------------------------
# GOV
# ------------------------------------------------------

 

function Ax_S_S_gov_ieee_tgov1_cb( gov )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values

    #      xg1               xg2    
    Ax = [-1/T1              0;
          (1 - T2/T3)/T3   (-1/T3)]
    

    return Ax

end


function Ax_S_Sgen_gov_ieee_tgov1_cb( gov; ω_ref0 = ω_ref0 )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values
    
    #     δ   u_gen_ω              e_d_dash e_q_dash
    
    Ac = [0  (-1/(T1*R *ω_ref0 ))     0       0;
          0   0                       0       0
          ]

    return Ac

end


function Bx_S_gov_ieee_tgov1_cb( gov )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values

    #      id  iq  ph  vh
    
     Bx = [0   0   0   0;
           0   0   0   0]

    return Bx

end


function Cx_S_gov_ieee_tgov1_cb( gov; ω_ref0 = ω_ref0 )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values

    
    #     ωs        ω_ref0           v_ref0  p_order0
    
    Cx = [0     1/(T1 * R * ω_ref0 )   0      1/T1;
          0            0               0       0
          ]

    return Cx

end


# ------------------------------------------------------
# Algebraic form for gov_ieee_tgov1_cb
#
# algebraic vars:
#    :τm_tilade, :ω_ref, :phat_in, :p_in
#
# state vars:
#    :xg1, :xg2, :xg3
# ------------------------------------------------------



"""
Coupling between state variables and algebraic variables
of gov `gov_ieee_tgov1`

     τm_tilade    #ω_ref

xg1   0        1/(R * ω_ref0)  

xg2   0          0

"""
function Ax_S_A_gov_ieee_tgov1_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values

    if ref_in_state == false

        #  τm_tilade

        Ac = [0   ;
              0   ]

        return Ac
        
    else    

        #    τm_tilade    #ω_ref

        Ac = [0         1/(R * ω_ref0);
              0          0]


        return Ac
        
    end

end


"""
Coupling between algebraic variables and state variables
of gov `gov_ieee_tgov1`

            xg1       xg2

 τm_tilade   T2/T3     1       

# ω_ref      0         0


"""
function Ax_A_S_gov_ieee_tgov1_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values

    if ref_in_state == false

        #      xg1        xg2

        Ac = [T2/T3        1]


        return Ac
        
    else

        #      xg1        xg2

        Ac = [T2/T3        1
              0            0]


        return Ac
                
    end
    
    

end



"""
Algebraic variables of gov `gov_ieee_tgov1`

          τm_tilade     ω_ref

 τm_tilade   -1         Dt/ωs     

  ω_ref      0           -1 


"""
function Ax_A_A_gov_ieee_tgov1_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values

    if ref_in_state == false

        #            τm_tilade

        Ax_algb = [    -1    ]     


        return Ax_algb
        
    else

        #            τm_tilade,   ω_ref

        Ax_algb = [    -1         Dt/ω_ref0 ;
                        0           -1 ]     

        return Ax_algb
        
    end
    

end


"""
Coupling between algebraic variables of gov `gov_ieee_tgov1`,
and state equations of generator

             δ   u_gen_ω       e_d_dash e_q_dash

 τm_tilade   0  (-Dt/ω_ref0)          0       0     

 ω_ref       0    0               0       0


"""
function Ax_A_Sgen_gov_ieee_tgov1_cb(
    gov;  ω_ref0 = ω_ref0, ref_in_state = false )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values

    if ref_in_state == false
    
        #     δ   u_gen_ω  e_d_dash  e_q_dash

        Ac = [0  (-Dt/ω_ref0)     0        0
              ]

        return Ac
        
    else
    
        #     δ   u_gen_ω    e_d_dash e_q_dash

        Ac = [0  (-Dt/ω_ref0)       0       0;
              0   0             0       0
              ]

        return Ac
        
    end

end


"""
Coupling between algebraic variables of gov `gov_ieee_tgov1`,
and network injections  id,  iq,  ph,  vh,

              id  iq  ph  vh

 τm_tilade   0    0   0   0     

 ω_ref       0    0   0   0


"""
function Bx_A_gov_ieee_tgov1_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values

    if ref_in_state == false

        #      id  iq  ph  vh

        Bx = [ 0    0   0  0 ]

        return Bx
        
    else
    
        #      id  iq  ph  vh

        Bx = [ 0    0   0  0;  
               0    0   0  0]

        return Bx
        
    end

end


"""
Coupling between algebraic variables of gov `gov_ieee_tgov1`,
and set points ωs, ω_ref0, v_ref0, p_order0

            ωs  ω_ref0  v_ref0  p_order0

 τm_tilade   0  Dt/(ωs)   0        0     

 ω_ref       0    1       0        0


"""
function Cx_A_gov_ieee_tgov1_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values

    if ref_in_state == false

        #      ωs        ω_ref0    v_ref0  p_order0

        Cx = [ 0          Dt/(ω_ref0)     0     0 ]


        return Cx
        
    else

        #      ωs        ω_ref0    v_ref0  p_order0

        Cx = [ 0          Dt/(ω_ref0)     0     0;
               0           1          0     0 ]


        return Cx
        
    end
    

end

# ------------------------------------------------------


function stab_Ac_gen_ieee_tgov1_cb(
    gov; ωs = ωs )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values
    
    #       u_gen_ω                e_d_dash e_q_dash
    Ac = [  -1/(T1 * R  * ωs )     0       0;
             0                        0       0
          ]

    return Ac

end

# ------------------------------------------------------

function output_coeff_gov_states_gov_ieee_tgov1_cb(
    gov; ωs = ωs )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values

    #       xg1      xg2  
    out = [(T2/T3)   1  ]

    return out

end


function output_coeff_gen_states_gov_ieee_tgov1_cb(
    gov, gen; ωs = ωs )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values
    
    out = zeros(1, length(gen.state_vars_syms))
    
    ω_idx = gen.dict_state_syms[:ω]

    out[1, ω_idx]  = -Dt/ω_ref0

    return out

end


function output_coeff_inputs_gov_ieee_tgov1_cb(
    gov; ωs = ωs )

    T1, T2, T3, Dt, p_max,p_min, R = gov.param_values
    
    input_vector  = [:ω_ref, :τm,  :v_ref ]

    out = [Dt/ω_ref0 0  0]

    return out

end

# ------------------------------------------------------
# ------------------------------------------------------

function Ax_S_S_gov_t0_cb( gov  )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    #    xg1                      xg2          xg3 
    Ax = [-1/Ts                  0           0;
          (1 - T3/Tc)/Tc       (-1/Tc)           0;
          (1/T5)*(1-T4/T5)*(T3/Tc) (1/T5)*(1-T4/T5) (-1/T5) ]

    return Ax

end



function Ax_S_Sgen_gov_t0_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false
    
        #     δ     u_gen_ω       e_d_dash  e_q_dash

        Ac = [0       0              0        0;
              0       0              0        0;
              0       0              0        0
              ]

        return Ac
        
    else

        #     δ   u_gen_ω              e_d_dash e_q_dash

        Ac = [0  (-1/(Ts * R * ω_ref0))   0       0;
              0   0                       0       0;
              0   0                       0       0
              ]

        return Ac

    end

end



function Bx_S_gov_t0_cb( gov )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    #      id  iq  ph  vh
    
     Bx = [0   0   0   0;
           0   0   0   0;
           0   0   0   0]

    return Bx

end


function Cx_S_gov_t0_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false

        #     ωs      ω_ref0       v_ref0  p_order0

        Cx = [0   1/(Ts*R *ω_ref0)   0      0;
              0          0           0      0;
              0          0           0      0
              ]

        return Cx
        
    else
    
        #     ωs           ω_ref0       v_ref  p_order0

        Cx = [0        1/(Ts*R*ω_ref0)   0    1/Ts;
              0             0            0      0;
              0             0            0      0
              ]

        return Cx
        
    end

end



# ------------------------------------------------------
# Algebraic form for gov_t0_cb
#
# algebraic vars:
#    :τm_tilade, :ω_ref, :phat_in, :p_in
#
# state vars:
#    :xg1, :xg2, :xg3
# ------------------------------------------------------

"""
Coupling between state variables and algebraic variables
of gov `gov_t0_cb`

     τm_tilade     ω_ref    phat_in   p_in

xg1   0             0        1/Ts      0       

xg2   0             0         0        0

xg3   0             0         0        0

"""
function Ax_S_A__gov_t0_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false

        #   τm_tilade phat_in   p_in

        Ac = [0       1/Ts        0 ;
              0        0          0 ;
              0        0          0    ]

        return Ac
        
    else
    
      # τm_tilade ω_ref phat_in   p_in
        
        Ac = [0     0   1/Ts  0  ;
              0     0    0    0   ;
              0     0    0    0 ]

        return Ac
        
    end

end


"""
Coupling between algebraic variables and state variables
of gov `gov_t0_cb`

            xg1                   xg2     xg3

 τm_tilade   (T4/T5)* (T3/Tc)    T4/T5     1    

 ω_ref       0                     0       0

phat_in      0                     0       0

 p_in        0                     0       0

"""
function Ax_A_S_gov_t0_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false

        #      xg1                   xg2     xg3

        Ac = [ (T4/T5)* (T3/Tc)    T4/T5     1;
               0                     0       0;
               0                     0       0]


        return Ac
        
    else
    
        #      xg1                   xg2     xg3

        Ac = [ (T4/T5)* (T3/Tc)    T4/T5     1;
               0                     0       0;
               0                     0       0;
               0                     0       0]


        return Ac
        
    end

end


"""
Algebraic variables of gov `gov_t0_cb`

            τm_tilade  #ω_ref      phat_in   p_in

 τm_tilade   -1           0           0       0

 #ω_ref       0          -1           0       0

 phat_in      0       1/(R*ωs)       (-1)     0

 p_in         0           0           1      (-1)

"""
function Ax_A_A_gov_t0_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false

        #            τm_tilade   phat_in   p_in

        Ax_algb = [   -1         0       0;
                      0         -1       0;
                      0          1       -1 ]     


        return Ax_algb
        
    else

        #            τm_tilade   ω_ref        phat_in   p_in
  
        Ax_algb = [   -1           0            0       0;
                      0          -1             0       0;
                      0       1/(R*ω_ref0 )    -1       0;
                      0           0             1       -1 ]     

        return Ax_algb
        
    end

end


"""
Coupling between algebraic variables of gov `gov_t0_cb`,
and state equations of generator

            δ      u_gen_ω     e_d_dash   e_q_dash

 τm_tilade  0        0            0          0    

 # ω_ref    0        0            0          0

 phat_in    0    -1/(R*ω)         0          0

 p_in       0        0            0          0


"""
function Ax_A_Sgen_gov_t0_cb(
    gov;  ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false
    
        #     δ   u_gen_ω           e_d_dash   e_q_dash

        Ac = [0        0               0          0;
              0   (-1/(R*ω_ref0))      0          0;
              0        0               0          0]

        return Ac
        
    else
    
        #     δ   u_gen_ω          e_d_dash   e_q_dash

        Ac = [0        0               0          0;
              0        0               0          0;
              0   (-1/(R*ω_ref0))      0          0;
              0        0               0          0]

        return Ac
        
    end

end



"""
Coupling between algebraic variables of gov `gov_t0_cb`,
and network injections  id,  iq,  ph,  vh,

              id  iq  ph  vh

 τm_tilade   0    0   0   0     

 #ω_ref      0    0   0   0

 phat_in     0    0   0   0

 p_in        0    0   0   0


"""
function Bx_A_gov_t0_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    if ref_in_state == false

        #      id  iq  ph  vh

        Bx = [ 0    0   0  0;
               0    0   0  0;
               0    0   0  0]


        return Bx
        
    else

        T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

        #      id  iq  ph  vh

        Bx = [ 0    0   0  0;
               0    0   0  0;
               0    0   0  0;
               0    0   0  0]    

        return Bx
        
    end

end


"""
Coupling between algebraic variables of gov `gov_t0_cb`,
and set points ωs, ω_ref0, v_ref0, p_order0

            ωs   #ω_ref0      v_ref0  p_order0

 τm_tilade   0     0            0        0     

 #ω_ref      0     1            0        0

 phat_in     0  1/(R*ωs)        0        1

 p_in        0     0            0        0


"""
function Cx_A_gov_t0_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false

        #      ωs     ω_ref0        v_ref0  p_order0

        Cx = [ 0        0              0       0;
               0   1/(R*ω_ref0)        0       1;
               0        0              0       0]


        return Cx
        
    else
    
        #      ωs     ω_ref0     v_ref0  p_order0

        Cx = [ 0        0           0       0;
               0        1           0       0;
               0   1/(R*ω_ref0)     0       1;
               0        0           0       0]


        return Cx
        
    end

end


# ------------------------------------------------------


function stab_Ac_gen_t0_cb( gov; ωs = ωs )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values
    
    #       u_gen_ω       e_d_dash e_q_dash
    Ac = [ -1/(Ts * R * ωs)   0       0;
            0                 0       0;
            0                 0       0
          ]

    return Ac

end

# ------------------------------------------------------

function output_coeff_gov_states_gov_t0_cb( gov  )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    #       xg1                  xg2    xg3 
    out = [ (T4/T5) * (T3/Tc)  (T4/T5)   1 ]

    return out

end


function output_coeff_gen_states_gov_t0_cb( gov, gen; ωs = ωs )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

  
    return zeros(1, length(gen.state_vars_syms)) 

end


function output_coeff_inputs_gov_t0_cb( gov; ωs = ωs )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    input_vector  = [:ω_ref, :τm,  :v_ref ]
  
    return zeros(1, length(input_vector))

end



# ------------------------------------------------------
# ------------------------------------------------------


function Ax_S_S_gov_t1_cb( gov  )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    #    xg1                        xg2            xg3 
    Ax = [-1/Ts                        0               0;
          (1 - T3/Tc)/Tc          ( -1/Tc)             0;
          (1/T5)*(1-T4/T5)*(T3/Tc) (1/T5)*(1-T4/T5)  (-1/T5)]

    return Ax

end


function Ax_S_Sgen_gov_t1_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false

        #     δ   u_gen_ω   e_d_dash e_q_dash

        Ac = [0    0            0       0;
              0    0            0       0;
              0    0            0       0
              ]

        return Ac
        
    else
    
        #     δ   u_gen_ω              e_d_dash e_q_dash
        Ac = [0  (-1/(Ts * R * ω_ref0 ))   0       0;
              0   0                        0       0;
              0   0                        0       0
              ]

        return Ac
        
    end

end



function Bx_S_gov_t1_cb( gov  )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    #      id  iq  ph  vh
     Bx = [0   0   0   0;
           0   0   0   0;
           0   0   0   0]

    return Bx

end



function Cx_S_gov_t1_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false
    
        #     ωs      ω_ref0       v_ref0  p_order0

        Cx = [0          0            0      0;
              0          0            0      0 ;
              0          0            0      0
              ]

        return Cx
        
    else
    
        #     ωs     ω_ref0          v_ref  p_order0

        Cx = [0     1/(Ts*R*ω_ref0)    0     1/Ts;
              0           0            0      0  ;
              0           0            0      0
              ]

        return Cx
        
    end

end


# ------------------------------------------------------
# Algebraic form for gov_t1_cb
#
# algebraic vars:
#    :τm_tilade, :ω_ref, :phat_in, :p_in
#
# state vars:
#    :xg1, :xg2, :xg3
# ------------------------------------------------------

"""
Coupling between state variables and algebraic variables
of gov `gov_t1_cb`

     τm_tilade     #ω_ref    phat_in 

xg1   0             0         1    

xg2   0             0         0     

xg3   0             0         0     

"""
function Ax_S_A_gov_t1_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false

        # τm_tilade  phat_in

        Ac = [0      1;
              0      0;
              0      0]

        return Ac
        
    else
    
       # τm_tilade ω_ref  phat_in

        Ac = [0     0    1;
              0     0    0;
              0     0    0]

        return Ac
        
    end

end


"""
Coupling between algebraic variables and state variables
of gov `gov_t1_cb`

            xg1                   xg2     xg3

 τm_tilade   (T4/T5)* (T3/Tc)    T4/T5     1    

 # ω_ref     0                     0       0

phat_in      0                     0       0


"""
function Ax_A_S_gov_t1_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false

        #      xg1                  xg2     xg3

        Ac = [ (T4/T5)* (T3/Tc)    T4/T5     1;
               0                     0       0 ]


        return Ac
        
    else
    
        #      xg1                   xg2     xg3

        Ac = [ (T4/T5)* (T3/Tc)    T4/T5     1;
               0                     0       0;
               0                     0       0 ]    

        return Ac

    end

end



"""
Algebraic variables of gov `gov_t1_cb`

            τm_tilade   #ω_ref      phat_in

 τm_tilade   -1           0           0   

 #ω_ref       0          -1           0   

 phat_in      0           0          -1  


"""
function Ax_A_A_gov_t1_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false

        #       τm_tilade     phat_in

        Ax_algb = [   -1       0;  
                      0       -1]

        return Ax_algb
        
    else
    
        #          τm_tilade   ω_ref      phat_in

        Ax_algb = [   -1        0              0;  
                      0        -1              0;  
                      0         0             -1]

        return Ax_algb

    end

end


"""
Coupling between algebraic variables of gov `gov_t1_cb`,
and state equations of generator

            δ      u_gen_ω     e_d_dash   e_q_dash

 τm_tilade  0        0            0          0    

 #ω_ref     0        0            0          0

 phat_in    0   -1/(R*ωs)         0          0


"""
function Ax_A_Sgen_gov_t1_cb(
    gov;  ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false

        #     δ   u_gen_ω         e_d_dash   e_q_dash

        Ac = [0        0             0          0;
              0   -1/(R*ω_ref0)      0          0 ]

        return Ac
        
    else
    
        #     δ   u_gen_ω       e_d_dash   e_q_dash

        Ac = [0        0           0          0;
              0        0           0          0;
              0    -1/(R*ω_ref0)   0          0 ]

        return Ac

    end

end


"""
Coupling between algebraic variables of gov `gov_t1_cb`,
and network injections  id,  iq,  ph,  vh,

              id  iq  ph  vh

 τm_tilade    0    0   0   0     

 #ω_ref       0    0   0   0

 phat_in      0    0   0   0


"""
function Bx_A_gov_t1_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false

        #      id  iq  ph  vh

        Bx = [ 0    0   0  0;  
               0    0   0  0
               ]

        return Bx
        
    else
    
        #      id  iq  ph  vh

        Bx = [ 0    0   0  0;  
               0    0   0  0;
               0    0   0  0
               ]

        return Bx

    end

end



"""
Coupling between algebraic variables of gov `gov_ieee_tgov1`,
and set points ωs, ω_ref0, v_ref0, p_order0

            ωs  ω_ref0  v_ref0  p_order0

 τm_tilade   0    0       0        0     

 #ω_ref      0    1       0        0

 phat_in     0  1/(R*ωs)  0        1


"""
function Cx_A_gov_t1_cb(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    if ref_in_state == false

        #      ωs      ω_ref0      v_ref0  p_order0

        Cx = [ 0          0           0     0;
               0       1/(R*ω_ref0)   0     1
               ]


        return Cx

    else
    
        #      ωs      ω_ref0        v_ref0  p_order0

        Cx = [ 0          0             0     0;
               0          1             0     0;
               0       1/(R*ω_ref0 )    0     1
               ]

        return Cx

    end

end


# ------------------------------------------------------


function stab_Ac_gen_t1_cb( gov; ωs = ωs )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values
    
    #       u_gen_ω        e_d_dash e_q_dash
    Ac = [ -1/(Ts * R * ωs)   0       0;
            0                 0       0;
            0                 0       0
          ]

    return Ac

end


function output_coeff_gov_states_gov_t1_cb( gov   )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values


    #       xg1                  xg2    xg3 
    out = [ (T4/T5) * (T3/Tc)  (T4/T5)   1 ]

    return out

end


function output_coeff_gen_states_gov_t1_cb( gov, gen; ωs = ωs )

    T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values
  
    return zeros(1, length(gen.state_vars_syms)) 

end


function output_coeff_inputs_gov_t1_cb( gov; ωs = ωs )

     T3, T4, T5, Tc, Ts, p_max, p_min, R = gov.param_values

    input_vector  = [:ω_ref, :τm,  :v_ref ]
  
    return zeros(1, length(input_vector))

end


# ------------------------------------------------------
# ------------------------------------------------------


function Ax_S_S_gov_t1_cb_sauer( gov  )

    Tc, Ts, p_max,p_min, R = gov.param_values

    #    xg1         xg2
    
    Ax = [-1/Ts      0;
          1/Tc     (-1/Tc)]

    return Ax

end


function Ax_S_Sgen_gov_t1_cb_sauer( gov;ω_ref0 = ω_ref0  )

    Tc, Ts, p_max,p_min, R = gov.param_values
    
    #     δ   u_gen_ω               e_d_dash e_q_dash
    
    Ac = [0  (-1/(Ts * R * ω_ref0 ))   0       0;
          0   0                        0       0]

    return Ac

end



function Bx_S_gov_t1_cb_sauer(gov )

    Tc, Ts, p_max,p_min, R = gov.param_values

    #      id  iq  ph  vh
    
     Bx = [0   0   0   0;
           0   0   0   0]

    return Bx

end


function Cx_S_gov_t1_cb_sauer( gov; ω_ref0 = ω_ref0 )

    Tc, Ts, p_max,p_min, R = gov.param_values

    #    ωs        ω_ref0          v_ref0  p_order0
    
    Cx = [0   1/(Ts * R * ω_ref0 )    0      1/Ts  ;
          0           0               0        0 ]

    return Cx

end


# ------------------------------------------------------
# Algebraic form for gov_t1_cb_sauer
#
# algebraic vars:
#    :τm_tilade, :ω_ref, :phat_in, :p_in
#
# state vars:
#    :xg1, :xg2, :xg3
# ------------------------------------------------------

"""
Coupling between state variables and algebraic variables
of gov `gov_t1_cb_sauer`

           ω_ref            τm_tilade            

xg1    1/(R * ω_ref0)         0   

xg2         0                 0    


"""
function Ax_S_A_gov_t1_cb_sauer(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    Tc, Ts, p_max,p_min, R = gov.param_values

    if  ref_in_state == false
    
        #   τm_tilade

        Ac = [0; 
              0

              ]

        return Ac

    else
    
           # ω_ref            τm_tilade

        Ac = [1/(R * ω_ref0)    0; 
              0                 0

              ]

        return Ac

    end

end


"""
Coupling between algebraic variables and state variables
of gov `gov_t1_cb_sauer`

            xg1     xg2 

 ω_ref       0        0  

τm_tilade    0        1

"""
function Ax_A_S_gov_t1_cb_sauer(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    Tc, Ts, p_max,p_min, R = gov.param_values

    if ref_in_state == false

        #      xg1     xg2    

        Ac = [  0      1]


        return Ac
        
    else
    
        #      xg1    xg2    

        Ac = [ 0        0;
               0        1]


        return Ac
        
    end

end


"""
Algebraic variables of gov `gov_t1_cb_sauer`

             #ω_ref   τm_tilade


 # ω_ref        -1       0

 τm_tilade      0       -1


"""
function Ax_A_A_gov_t1_cb_sauer(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    Tc, Ts, p_max,p_min, R = gov.param_values

    if ref_in_state == false

        #         τm_tilade  

        Ax_algb = [  -1 ]

        return Ax_algb
        
    else

        #           ω_ref   τm_tilade  

        Ax_algb = [   -1      0;
                       0     -1]

        return Ax_algb
        
    end

end


"""
Coupling between algebraic variables of gov `gov_t1_cb_sauer`,
and state equations of generator

            δ      u_gen_ω     e_d_dash   e_q_dash

 ω_ref      0        0            0          0

 τm_tilade  0        0            0          0

"""
function Ax_A_Sgen_gov_t1_cb_sauer(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    Tc, Ts, p_max,p_min, R = gov.param_values

    if  ref_in_state == false
    
        #     δ   u_gen_ω       e_d_dash   e_q_dash

        Ac = [ 0        0           0          0]

        return Ac

    else
        
        #     δ   u_gen_ω       e_d_dash   e_q_dash

        Ac = [0        0           0          0;
              0        0           0          0]

        return Ac

    end

end



"""
Coupling between algebraic variables of gov `gov_t1_cb_sauer`,
and network

            id      iq            pg         vh

 ω_ref      0        0            0          0

 τm_tilade  0        0            0          0

"""
function Bx_A_gov_t1_cb_sauer(
    gov; ω_ref0 = ω_ref0, ref_in_state = false )

    Tc, Ts, p_max,p_min, R = gov.param_values
    
    if ref_in_state == false

        #      id  iq  ph  vh

        Bx = [  0    0   0  0]

        return Bx
        
    else

        #      id  iq  ph  vh

        Bx = [ 0    0   0  0;
               0    0   0  0]

        return Bx

    end

end


"""
Coupling between algebraic variables of gov `gov_t1_cb_sauer`,
and set points

            ωs     ω_ref0    v_ref0  p_order0

 #ω_ref      0        1         0        0

τm_tilade    0        0         0        0

"""
function Cx_A_gov_t1_cb_sauer( gov;
                                  ω_ref0 = ω_ref0,
                                  ref_in_state = false )

    Tc, Ts, p_max,p_min, R = gov.param_values

    if ref_in_state == false

        #      ωs      ω_ref0      v_ref0  p_order0

        Cx = [ 0          0           0     0]


        return Cx
        
    else
    
        #      ωs      ω_ref0      v_ref0  p_order0

        Cx = [ 0          1           0     0;
               0          0           0     0]


        return Cx

    end

end



# ------------------------------------------------------


function stab_Ac_gov_gen_sauer( gov; ωs = ωs  )

    Tc, Ts, p_max,p_min, R = gov.param_values
    
    #       u_gen_ω         e_d_dash e_q_dash
    Ac = [ -1/(Ts * R * ωs )   0       0;
            0                  0       0]

    return Ac

end

# ------------------------------------------------------


function output_coeff_gov_states_gov_sauer( gov  )

    Tc, Ts, p_max,p_min, R = gov.param_values

    #    xg1       xg2    
    out = [0     1]

    return out

end


function output_coeff_gen_states_gov_sauer( gov, gen; ωs = ωs )

    Tc, Ts, p_max,p_min, R = gov.param_values
    
  
    return zeros(1, length(gen.state_vars_syms)) 

end


function output_coeff_inputs_gov_sauer( gov; ωs = ωs )

    Tc, Ts, p_max,p_min, R = gov.param_values

    input_vector  = [:ω_ref, :τm,  :v_ref ]
  
    return zeros(1, length(input_vector))

end


# ------------------------------------------------------
# AVR
# ------------------------------------------------------


function Ax_S_S_avr_t0_cb( avr, vf_tilade, vr1 )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    #     vm,    vr1       vr2    vf_tilade
    Ax = [-1/Tr     0      0        0;
          -Ka/Ta  (-1/Ta)  Ka/Ta    (-(Ka * Kf/Tf)/Ta) ;
          0        0       (-1/Tf)  (Kf/Tf)/Tf;
          0        1/Te     0       -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te
          ]    
    
    return Ax

end


function Ax_S_Sgen_avr_t0_cb( avr )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values
    
    #     δ   u_gen_ω   e_d_dash e_q_dash
    Ac = [0     0            0       0;
          0     0            0       0;
          0     0            0       0;
          0     0            0       0
          ]

    return Ac

end


function Bx_S_avr_t0_cb( avr )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    #     id  iq  ph   vh
    Bx = [0   0   0    1/Tr; 
          0   0   0    0;    
          0   0   0    0;    
          0   0   0    0   
          ]
    
    return Bx

end


function Cx_S_avr_t0_cb( avr )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    #    ωs    ω_ref0   v_ref0  p_order0
    
    Cx = [0        0      0       0 ;    
          0        0      Ka/Ta   0 ;
          0        0      0       0;    
          0        0      0       0
          ]
    
    return Cx

end


# ------------------------------------------------------
# Algebraic form for avr_t0_cb!
#
# algebraic vars:
#    v_ref, vr1_hat 
#     
#
# state vars:
#   vm, vr1, vr2, vf_tilade,
#
# ------------------------------------------------------

"""
Coupling between state variables and algebraic variables
of avr `avr_t0_cb!`

             v_ref           vr1_hat     vf

 vm            0              0           0

 vr1           0              0           0

 vr2           0              0           0

 vf_tilade     0              1/Te        0


"""
function Ax_S_A_avr_t0_cb( avr; ref_in_state = false )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #      vr1_hat   vf

        Ac = [    0      0; 
                  0      0;
                  0      0;
                  1/Te   0
              ]

        return Ac
        
    else

        #   v_ref   vr1_hat    vf

        Ac = [0        0       0;    
              0        0       0;
              0        0       0;
              0        1/Te    0
              ]

        return Ac

    end

end



"""
Coupling between algebraic variables and state variables
of avr `avr_t0_cb`

         vm     vr1    vr2   vf_tilade  
         
v_ref     0      0      0      0

vr1_hat   0      1      0      0

vf        0      0      0      1


"""
function Ax_A_S_avr_t0_cb( avr; ref_in_state = false  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #      vm   vr1     vr2   vf_tilade    

        Ac = [  0     1      0       0;
                0     0      0       1]

        return Ac
        
    else
    
        #      vm   vr1     vr2   vf_tilade    

        Ac = [  0     0      0       0;
                0     1      0       0;
                0     0      0       1]

        return Ac

    end

end



"""
Algebraic variables of avr `avr_t0_cb`

             v_ref    vr1_hat   vf

v_ref         -1        0       0
 
vr1_hat        0       -1       0

vf             0        0      -1


"""
function Ax_A_A_avr_t0_cb( avr;  ref_in_state = false  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #           vr1_hat  vf

        Ax_algb = [   -1     0;
                       0    -1 ]

        return Ax_algb
        
    else
    
        #       v_ref   vr1_hat  vf

        Ax_algb = [ -1      0     0;
                     0     -1     0;
                     0      0    -1 ]

        return Ax_algb

    end

end



"""
Coupling between algebraic variables of avr `avr_t0_cb`,
and state equations of generator

            δ      u_gen_ω     e_d_dash   e_q_dash

 v_ref      0        0            0          0

 vr1_hat    0        0            0          0

 vf         0        0            0          0 

"""
function Ax_A_Sgen_avr_t0_cb( avr;  ref_in_state = false )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false
    
        #     δ   u_gen_ω    e_d_dash   e_q_dash

        Ac = [0       0         0          0;
              0       0         0          0]

        return Ac
        
    else
    
        #     δ   u_gen_ω    e_d_dash   e_q_dash

        Ac = [0       0         0          0;
              0       0         0          0;
              0       0         0          0]

        return Ac

    end

end



"""
Coupling between algebraic variables of avr `avr_t0_cb`,
and network

            id      iq            pg         vh

 v_ref      0        0            0          0

 vr1_hat    0        0            0          0

 vf         0        0            0          0

"""
function Bx_A_avr_t0_cb( avr; ref_in_state = false )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #     id  iq  ph   vh
        Bx = [
              0   0   0    0;    
              0   0   0    0   
              ]

        return Bx
        
    else
    
        #     id  iq  ph   vh
        Bx = [
              0   0   0    0;    
              0   0   0    0;
              0   0   0    0           
              ]

        return Bx

    end

end



"""
Coupling between algebraic variables of avr `avr_t0_cb`,
and set points

             ωs    ω_ref0       v_ref0   p_order0

 v_ref      0        0            1          0

 vr1_hat    0        0            0          0

 vf         0        0            0          0

"""
function Cx_A_avr_t0_cb( avr; ref_in_state = false )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #     ωs    ω_ref0   v_ref  p_order0

        Cx = [ 0        0      0       0 ;
               0        0      0       0 
              ]

        return Cx
        
    else
           
        #     ωs    ω_ref0   v_ref  p_order0

        Cx = [0        0      1       0 ;    
              0        0      0       0 ;
              0        0      0       0
              ]

        return Cx
        
    end

end


# ------------------------------------------------------


function stab_Ac_avr_gen_t0_cb( avr )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values
    
    #     u_gen_ω   e_d_dash e_q_dash
    Ac =[   0            0       0;
            0            0       0;
            0            0       0;
            0            0       0
            ]

    return Ac

end

# ------------------------------------------------------


function output_coeff_avr_states_avr_t0_cb( avr )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    #      vm   vr1  vr2   vf_tilade    
    out = [0     0    0        1]

    return out

end


function output_coeff_gen_states_avr_t0_cb( avr, gen )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values
  
    return zeros(1, length(gen.state_vars_syms)) 

end


function output_coeff_inputs_avr_t0_cb( avr )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values
    input_vector  = [:ω_ref, :τm,  :v_ref ]
  
    return zeros(1, length(input_vector))

end


# ------------------------------------------------------

function Ax_S_S_avr_t1_cb(avr, vf_tilade, vr1 )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    #     vm,     vr1     vr2    vf_tilade
    Ax = [-1/Tr     0      0         0;
          -Ka/Ta  (-1/Ta)  Ka/Ta     (-Ka * Kf)/(Tf*Ta);
          0        0       (-1/Tf)    Kf/Tf^2 ;
          0        1/Te    0         -1*(Ke + Sevf(Ae,Be,vf_tilade))/Te
          ]    
    
    return Ax

end


function Ax_S_Sgen_avr_t1_cb( avr )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values
    
    #     δ   u_gen_ω   e_d_dash e_q_dash
    Ac = [0     0            0       0;
          0     0            0       0;
          0     0            0       0;
          0     0            0       0
          ]

    return Ac

end



function Bx_S_avr_t1_cb( avr  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    #     id, iq, ph,   vh 
    Bx = [0    0   0    1/Tr; 
          0    0   0    0;    
          0    0   0    0;    
          0    0   0    0    
          ]
    
    return Bx

end


function Cx_S_avr_t1_cb( avr )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    #   ωs  ω_ref0  v_ref0  p_order0
    
    Cx = [0     0     0         0;
          0     0     Ka/Ta     0;
          0     0     0         0;
          0     0     0         0
          ]
    
    return Cx

end


# ------------------------------------------------------
# Algebraic form for avr_t1_cb!
#
# algebraic vars:
#    v_ref, vr1_hat 
#     
#
# state vars:
#   vm, vr1, vr2, vf_tilade,
#
# ------------------------------------------------------

"""
Coupling between state variables and algebraic variables
of avr `avr_t1_cb!`

             v_ref   vf

 vm            0     0

 vr1           0     0

 vr2           0     0

 vf_tilade     0     0


"""
function Ax_S_A_avr_t1_cb( avr; ref_in_state = false )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #       vf

        Ac = [  0; 
                0;
                0;
                0   
              ]

        return Ac
        
    else

        #     v_ref   vf
        
        Ac = [0      0; 
              0      0;
              0      0;
              0      0   
              ]

        return Ac

    end

end



"""
Coupling between algebraic variables and state variables
of avr `avr_t1_cb`

         vm     vr1    vr2   vf_tilade  

          0      0      0      0
v_ref

 vf       0      0      0      1

"""
function Ax_A_S_avr_t1_cb( avr; ref_in_state = false )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #      vm   vr1   vr2   vf_tilade    

        Ac = [  0     0     0       1]


        return Ac
        
    else
    
        #      vm   vr1   vr2   vf_tilade    

        Ac = [  0     0      0       0;
                0     0      0       1]


        return Ac

    end

end



"""
Algebraic variables of avr `avr_t1_cb`

             v_ref   vf

v_ref         -1     0

vf             0    -1 
 

"""
function Ax_A_A_avr_t1_cb( avr; ref_in_state = false  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #           vf 

        Ax_algb = [ -1 ]

        return Ax_algb
        
    else
    
        #       v_ref     vf 

        Ax_algb = [ -1    0;
                     0   -1]

        return Ax_algb

    end

end


"""
Coupling between algebraic variables of avr `avr_t1_cb`,
and state equations of generator

            δ      u_gen_ω     e_d_dash   e_q_dash

 v_ref      0        0            0          0

 vf         0        0            0          0


"""
function Ax_A_Sgen_avr_t1_cb( avr; ref_in_state = false )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values
    
    if ref_in_state == false
    
        #     δ   u_gen_ω    e_d_dash   e_q_dash

        Ac = [0       0         0          0]

        return Ac
        
    else
        #     δ   u_gen_ω    e_d_dash   e_q_dash

        Ac = [0       0         0          0;
              0       0         0          0]

        return Ac

    end
    
    

end



"""
Coupling between algebraic variables of avr `avr_t1_cb`,
and network

            id      iq            pg         vh

 v_ref      0        0            0          0

 vf         0        0            0          0


"""
function Bx_A_avr_t1_cb( avr; ref_in_state = false )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #      id      iq     ph     vh

        Bx = [  0       0      0      0
              ]

        return Bx
        
    else

        #     id     iq     ph     vh

        Bx = [
            0       0      0      0 ;
            0       0      0      0
              ]

        return Bx

    end
    

end



"""
Coupling between algebraic variables of avr `avr_t1_cb`,
and set points

             ωs    ω_ref0       v_ref0   p_order0

 v_ref      0        0            1          0

 vf         0        0            0          0


"""
function Cx_A_avr_t1_cb( avr; ref_in_state = false )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #     ωs    ω_ref0  v_ref0  p_order0

        Cx = [ 0       0      0       0
              ]

        return Cx
        
    else
    
        #     ωs    ω_ref0  v_ref0  p_order0

        Cx = [
            0        0      1       0;
            0        0      0       0
              ]

        return Cx

    end

end


# ------------------------------------------------------


function stab_Ac_avr_gen_t1_cb( avr )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values
    
    #      u_gen_ω   e_d_dash e_q_dash
    Ac = [   0            0       0;
             0            0       0;
             0            0       0;
             0            0       0
          ]

    return Ac

end

# ------------------------------------------------------

function output_coeff_avr_states_avr_t1_cb( avr  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    #      vm   vr1  vr2   vf_tilade    
    out = [0     0    0        1]

    return out

end


function output_coeff_gen_states_avr_t1_cb( avr, gen  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

  
    return zeros(1, length(gen.state_vars_syms)) 

end


function output_coeff_inputs_avr_t1_cb( avr  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values
    
    input_vector  = [:ω_ref, :τm,  :v_ref ]
  
    return zeros(1, length(input_vector))

end


# ------------------------------------------------------
# ------------------------------------------------------

function Ax_S_S_avr_t1_cb_sauer(avr, vf_tilade, vr1 )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    #     vr1    vr2    vf_tilade
    # Ax = [-threshold_limits(vr1, V_R_max, V_R_min)/(vr1 * Ta)  Ka/Ta  -(Ka * Kf/Tf)/Ta;
    #       0      -1/Tf  (Kf/Tf)/Tf;
    #       threshold_limits(vr1, V_R_max, V_R_min)/(vr1 * Te)    0     (-Ke + Sevf(Ae, Be, vf_tilade))/Te
    #       ] 

    
    #     vr1    vr2    vf_tilade
    
    Ax = [# 1/Ta  Ka/Ta    (-(Ka * Kf/Tf)/Ta);
          -1/Ta  Ka/Ta    -Ka * Kf/(Tf*Ta);
          0      (-1/Tf)  (Kf/Tf)/Tf;
          1/Te    0       -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te
          ]      
    return Ax

end


function Ax_S_Sgen_avr_t1_cb_sauer( avr  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values
    
    #     δ   u_gen_ω   e_d_dash e_q_dash
    Ac = [0     0            0       0;
          0     0            0       0;
          0     0            0       0
          ]

    return Ac

end



function Bx_S_avr_t1_cb_sauer( avr )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    #     id, iq, ph,   vh    
    Bx = [0   0   0  (-Ka/Ta);
          0   0   0    0;
          0   0   0    0
          ]
    
    return Bx

end


function Cx_S_avr_t1_cb_sauer( avr )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    #   ωs   ω_ref0  v_ref0  p_order0
    
    Cx = [0     0    Ka/Ta    0;
          0     0      0     0;
          0     0      0     0
          ]
    
    return Cx

end



# ------------------------------------------------------
# Algebraic form for avr_t1_cb_sauer!
#
# algebraic vars:
#    v_ref, vr1_hat 
#     
#
# state vars:
#   vm, vr1, vr2, vf_tilade,
#
# ------------------------------------------------------

"""
Coupling between state variables and algebraic variables
of avr `avr_t1_cb_sauer`

             v_ref      vf


 vr1           0        0

 vr2           0        0

 vf_tilade     0        0


"""
function Ax_S_A_avr_t1_cb_sauer( avr; ref_in_state = false )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #      vf

        Ac = [  0;
                0;
                0 
              ]

        return Ac
        
    else
    
        #   v_ref    vf

        Ac = [0     0;
              0     0;
              0     0 
              ]

        return Ac

    end

end



"""
Coupling between algebraic variables and state variables
of avr `avr_t1_cb_sauer`

        vr1    vr2   vf_tilade  

          
v_ref    0      0      0

vf       0      0      1


"""
function Ax_A_S_avr_t1_cb_sauer( avr; ref_in_state = false  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #      vr1   vr2   vf_tilade    

        Ac = [   0      0       1]


        return Ac
        
    else
    
        #      vr1   vr2   vf_tilade    

        Ac = [  0      0       0;
                0      0       1]


        return Ac

    end

end



"""
Algebraic variables of avr `avr_t1_cb_sauer`

             v_ref   vf

v_ref         -1     0

vf             0    -1
 

"""
function Ax_A_A_avr_t1_cb_sauer( avr; ref_in_state = false  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #           vf

        Ax_algb = [ -1 ]

        return Ax_algb
        
    else
    
        #       v_ref    vf

        Ax_algb = [ -1   0;
                    0   -1]

        return Ax_algb

    end

end


"""
Coupling between algebraic variables of avr `avr_t1_cb_sauer`
and state equations of generator

            δ      u_gen_ω     e_d_dash   e_q_dash

 v_ref      0        0            0          0

 vf         0        0            0          0


"""
function Ax_A_Sgen_avr_t1_cb_sauer( avr; ref_in_state = false  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false
    
        #     δ   u_gen_ω    e_d_dash   e_q_dash

        Ac = [0       0         0          0]


        return Ac
        
    else
    
        #     δ   u_gen_ω    e_d_dash   e_q_dash

        Ac = [0       0         0          0;
              0       0         0          0]

        return Ac
        
    end

end



"""
Coupling between algebraic variables of avr `avr_t1_cb_sauer`,
and network

            id      iq            pg         vh

 v_ref      0        0            0          0

 vf         0        0            0          0


"""
function Bx_A_avr_t1_cb_sauer( avr; ref_in_state = false )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #     id  iq  ph   vh

        Bx = [ 0   0   0    0 ]


        return Bx
        
    else
    
        #     id  iq  ph   vh
        Bx = [
               0   0   0    0 ;
               0   0   0    0
              ]

        return Bx

    end

end


"""
Coupling between algebraic variables of avr `avr_t1_cb_sauer`,
and set points

             ωs    ω_ref0       v_ref   p_order0

 v_ref       0        0            1          0

 vf          0        0            0          0


"""
function Cx_A_avr_t1_cb_sauer( avr; ref_in_state = false )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    if ref_in_state == false

        #     ωs    ω_ref0   v_ref0  p_order0

        Cx = [ 0        0      0       0]


        return Cx
        
    else
    
        #     ωs    ω_ref0   v_ref0  p_order0

        Cx = [ 0        0      1       0 ;
               0        0      0       0
              ]

        return Cx

    end

end


# ------------------------------------------------------


function stab_Ac_avr_gen_t1_cb_sauer( avr  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values
    
    #      u_gen_ω   e_d_dash e_q_dash
    Ac = [   0            0       0;
             0            0       0;
             0            0       0
          ]

    return Ac

end



function output_coeff_avr_states_avr_t1_cb_sauer( avr  )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values
    
    #      vr1  vr2   vf_tilade    
    out = [0     0        1]

    return out

end


function output_coeff_gen_states_avr_t1_cb_sauer( avr, gen )

    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values
    
    return zeros(1, length(gen.state_vars_syms)) 

end


function output_coeff_inputs_avr_t1_cb_sauer( avr  )
    
    Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be =
        avr.param_values

    input_vector  = [:ω_ref, :τm,  :v_ref ]
  
    return zeros(1, length(input_vector))

end


# ------------------------------------------------------
# Generators nodes matrix
# ------------------------------------------------------


function Ax_gen(gen, D, H, Ωb, ωs, ra, xℓ,
                X_d, X_q, X_d_dash, X_q_dash,
                X_d_2dash, X_q_2dash,
                T_d_dash, T_q_dash,
                T_d_2dash, T_q_2dash )

    ωs = gen.ωs
    
    Ax = zeros( length(gen.state_vars_syms),
                length(gen.state_vars_syms) )

    δ_idx       = gen.dict_state_syms[:δ]
    ω_idx       = gen.dict_state_syms[:ω]
    ed_dash_idx = gen.dict_state_syms[:ed_dash]
    eq_dash_idx = gen.dict_state_syms[:eq_dash]

    Ax[δ_idx, ω_idx]             = 1
    Ax[ω_idx, ω_idx]             = (-D * ωs /(2 * H))
    Ax[ed_dash_idx, ed_dash_idx] = -1/T_q_dash
    Ax[eq_dash_idx, eq_dash_idx] = -1/T_d_dash    

    # #     δ        ω        ed_dash     eq_dash
    
    # Ax = [0        1           0            0;
    #       0 (-D*ωs/(2 * H))    0            0;
    #       0        0         -1/T_q_dash    0;
    #       0        0           0           -1/T_d_dash
    #       ]

    return Ax 
    
end


function Ax_gen(D, H, Ωb, ωs, ra, xℓ,
                X_d, X_q, X_d_dash, X_q_dash,
                X_d_2dash, X_q_2dash,
                T_d_dash, T_q_dash,
                T_d_2dash, T_q_2dash )

    #     δ        ω          ed_dash     eq_dash
    
    Ax = [0        1           0              0;
          0 (-D * ωs/(2 * H))  0              0;
          0        0           (-1/T_q_dash)  0;
          0        0           0           (-1/T_d_dash)
          ]

    return Ax 
    
end


function Bx_gen(gen, D, H, Ωb, ωs, ra, xℓ,
                X_d, X_q, X_d_dash, X_q_dash,
                X_d_2dash, X_q_2dash,
                T_d_dash, T_q_dash,
                T_d_2dash, T_q_2dash )

    input_vector = [:id, :iq , :ph,  :vh]
    
    input_size = length( input_vector )

    Bx = zeros( length(gen.state_vars_syms),
                input_size )
    
    dict_input_syms_Idx::OrderedDict{Symbol, Int64} =
        OrderedDict(sym => ind for (ind, sym) in
                        enumerate( input_vector ) )
    
    id_idx = dict_input_syms_Idx[:id]
    iq_idx = dict_input_syms_Idx[:iq]
    ph_idx = dict_input_syms_Idx[:ph]
    vh_idx = dict_input_syms_Idx[:vh]
    
    δ_idx       = gen.dict_state_syms[:δ]
    ω_idx       = gen.dict_state_syms[:ω]
    ed_dash_idx = gen.dict_state_syms[:ed_dash]
    eq_dash_idx = gen.dict_state_syms[:eq_dash]

    Bx[ω_idx, ph_idx]       = -1* ωs/(2*H)
    Bx[ed_dash_idx, iq_idx] = (X_q - X_q_dash)/T_q_dash
    Bx[eq_dash_idx, id_idx] = -(X_d - X_d_dash)/T_d_dash
    
    # #     id      iq                           ph    vh
    
    # Bx = [ 0      0                            0      0;
    #        0      0                        -ωs/(2*H)  0;
    #        0      (X_q - X_q_dash)/T_q_dash    0      0;
    #       -(X_d - X_d_dash)/T_d_dash  0        0      0
    #        ]

    return Bx 

end



function Bx_gen( D, H, Ωb, ωs, ra, xℓ,
                 X_d, X_q, X_d_dash, X_q_dash,
                 X_d_2dash, X_q_2dash,
                 T_d_dash, T_q_dash,
                 T_d_2dash, T_q_2dash )

    #     id      iq                           ph          vh
    
    Bx = [ 0      0                            0            0;
           0      0                        (-1 * ωs/(2*H))  0;
           0      (X_q - X_q_dash)/T_q_dash    0            0;
          (X_d - X_d_dash)/T_d_dash  0         0            0
           ]

    return Bx 

end


function Cx_gen( gen, D, H, Ωb, ωs, ra, xℓ,
                 X_d, X_q, X_d_dash, X_q_dash,
                 X_d_2dash, X_q_2dash,
                 T_d_dash, T_q_dash,
                 T_d_2dash, T_q_2dash )

    ωs = gen.ωs

    input_vector = [:ωs, :ω_ref0, :v_ref0, :p_order0 ]
    
    input_size = length( input_vector )

    Cx = zeros( length(gen.state_vars_syms),
                input_size )
    
    dict_input_syms_Idx::OrderedDict{Symbol, Int64} =
        OrderedDict(sym => ind for (ind, sym) in
                        enumerate( input_vector ) )
    
    ωs_idx       = dict_input_syms_Idx[:ωs]
    ω_ref_idx    = dict_input_syms_Idx[:ω_ref0]
    v_ref_idx    = dict_input_syms_Idx[:v_ref0]
    p_order0_idx = dict_input_syms_Idx[:p_order0]
    
    δ_idx       = gen.dict_state_syms[:δ]
    ω_idx       = gen.dict_state_syms[:ω]
    ed_dash_idx = gen.dict_state_syms[:ed_dash]
    eq_dash_idx = gen.dict_state_syms[:eq_dash]

    """
    The output of gov is used to control τm, hence 1*ωs/(2*H)
    that should be in [dω_idx, τm_idx ] is set to zero

    """

    Cx[δ_idx, ωs_idx]  = -1
    Cx[ω_idx, ωs_idx]  = D * ωs /(2*H)
    
    # #     ωs       ω_ref0      v_ref0         p_order0  
    # Cx = [-1          0          0             0;
    #       ωs*D/(2*H)  0          0             0;
    #       0           0          0             0;
    #       0           0          0             0]
    
    return Cx 

end



function Cx_gen( D, H, Ωb, ωs, ra, xℓ,
                 X_d, X_q, X_d_dash, X_q_dash,
                 X_d_2dash, X_q_2dash,
                 T_d_dash, T_q_dash,
                 T_d_2dash, T_q_2dash )

    """
    The output of a gov is used to control τm, hence 1*ωs/(2*H)
    that should be in [dω_idx, τm_idx ] is set to zero

    """
    #     ωs           ω_ref0        v_ref0      p_order0
    Cx = [-1              0            0            0;
          D * ωs/(2*H)    0            0            0;
          0               0            0            0;
          0               0            0            0 ]

    return Cx 

end


function Ac_gen_gov( gen, gov )
    
    ac_gov = zeros(length(gen.state_vars_syms),
               length(gov.state_vars_syms) )

    ωs    = gen.ωs
    
    ω_idx = gen.dict_state_syms[:ω]
    
    H     = gen.H

    D     = gen.D

    """

    In some cases, the output may not be a pure state
    var,

    e.g. governors

    """

    if gov.output_sig_syms[ 1 ] ∈ gov.state_vars_syms

        out_sym_idx = gov.dict_state_syms[
            gov.output_sig_syms[ 1 ] ]
        
        ac_gov[ω_idx, out_sym_idx] = ωs /(2 * H)

    else
        for (a_state, a_state_coeff ) ∈ zip(
            gov.state_vars_syms,
            gov.output_sig_state_coeff)

            a_state_idx = gov.dict_state_syms[ a_state  ]

            ac_gov[ω_idx, a_state_idx] =
                a_state_coeff *  ωs /(2 * H)
        end

    end
        
    
    
    return ac_gov 

end



"""
Gen state vars coupling to state vars of gov
"""
function Ac_S_gen_gov( gen, gov )
    
    ac_gov = zeros(
        length( gen.state_vars_syms),
        length( gov.state_vars_syms) )

    ωs    = gen.ωs
    
    ω_idx = gen.dict_state_syms[:ω]
    
    H     = gen.H
    
    if gov.output_sig_syms[ 1 ] ∈ gov.state_vars_syms
        
        out_sym_idx = gov.dict_state_syms[
            gov.output_sig_syms[ 1 ] ]
        
        ac_gov[ω_idx, out_sym_idx] = ωs /(2 * H)
        
    end
    
    
    return ac_gov 

end


"""

Gen state vars coupling to im_algberaic_vars of gov

"""
function Ac_A_gen_gov( gen, gov )
    
    ac_gov = zeros(length( gen.state_vars_syms ),
               length( gov.im_algebraic_vars_syms ) )

    ωs    = gen.ωs
    
    H     = gen.H
        
    ω_idx = gen.dict_state_syms[:ω]
    
    if gov.output_sig_syms[ 1 ] ∈ gov.im_algebraic_vars_syms 

         τm_sym_idx =
             gov.dict_im_algebraic_vars_syms[
                 gov.output_sig_syms[ 1 ] ]
    
        ac_gov[ω_idx, τm_sym_idx] = ωs /(2 * H)

     end
    
    
    # if gov.output_sig_syms[ 1 ] ∈ gov.algebraic_vars_syms 
    #     out_sym_idx =
    #         gov.dict_im_algebraic_vars_syms[
    #             gov.output_sig_syms[ 1 ] ]
    #     ac_gov[ω_idx, out_sym_idx] = ωs /(2 * H)        
    # end
    
    return ac_gov 

end




"""

Gen state vars coupling to state of avr

"""
function Ac_gen_avr( gen, avr )
        
    ac_avr = zeros(length(gen.state_vars_syms),
               length(avr.state_vars_syms) )
    
    eq_dash_idx = gen.dict_state_syms[:eq_dash]
    T_d_dash    = gen.T_d_dash
    

   """
    In some cases, the output may not be a pure state var

    """


    if avr.output_sig_syms[ 1 ] ∈ avr.state_vars_syms

        vf_sym_idx  = avr.dict_state_syms[
            avr.output_sig_syms[ 1 ]]
        
        ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash

    else
        for (a_state, a_state_coeff ) ∈ zip(
            avr.state_vars_syms,
            avr.output_sig_state_coeff)

            a_state_idx = avr.dict_state_syms[ a_state  ]

            ac_avr[eq_dash_idx, a_state_idx] =
                a_state_coeff * 1/T_d_dash
        end

    end
    

    return ac_avr

end


"""

Gen state vars coupling to state of avr

"""
function Ac_S_gen_avr( gen, avr )
        
    ac_avr = zeros(length(gen.state_vars_syms),
               length(avr.state_vars_syms) )
    
    eq_dash_idx = gen.dict_state_syms[:eq_dash]
    T_d_dash    = gen.T_d_dash

     
    if avr.output_sig_syms[ 1 ] ∈ avr.state_vars_syms
        
        vf_sym_idx  = avr.dict_state_syms[
            avr.output_sig_syms[ 1 ] ]

        ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash
                
    end     

    return ac_avr

end


"""

Gen state vars coupling to im_algebraic_vars of avr

"""
function Ac_A_gen_avr( gen, avr )
        
    ac_avr = zeros(length(gen.state_vars_syms),
               length(avr.im_algebraic_vars_syms) )
    
    eq_dash_idx = gen.dict_state_syms[:eq_dash]
    T_d_dash    = gen.T_d_dash

    vf_sym_idx  = avr.dict_im_algebraic_vars_syms[
        avr.output_sig_syms[ 1 ] ]
    
    ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash
    
    # if avr.output_sig_syms[ 1 ] ∈ avr.algebraic_vars_syms     
    #     vf_sym_idx  = avr.dict_im_algebraic_vars_syms[
    #         avr.output_sig_syms[ 1 ] ]
    #     ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash        
    # end
     
    return ac_avr

end

#-----------------------------------------
#-----------------------------------------



"""
Gen state vars coupling to state vars of gov
"""
function SM_Ax_gen_S_gov_S( gen, gov )
    
    ac_gov = zeros(
        length( gen.state_vars_syms),
        length( gov.state_vars_syms) )

    ωs    = gen.ωs
    
    ω_idx = gen.dict_state_syms[:ω]
    
    H     = gen.H
    
    if gov.output_sig_syms[ 1 ] ∈ gov.state_vars_syms
        
        out_sym_idx = gov.dict_state_syms[
            gov.output_sig_syms[ 1 ] ]
        
        ac_gov[ω_idx, out_sym_idx] = ωs /(2 * H)
        
    end
    
    
    return ac_gov 

end


"""

Gen state vars coupling to im_algberaic_vars of gov

"""
function SM_Ax_gen_S_gov_A( gen, gov )
    
    ac_gov = zeros(length( gen.state_vars_syms ),
               length( gov.im_algebraic_vars_syms ) )

    ωs    = gen.ωs
    
    H     = gen.H
        
    ω_idx = gen.dict_state_syms[:ω]
    
    if gov.output_sig_syms[ 1 ] ∈ gov.im_algebraic_vars_syms 

         τm_sym_idx =
             gov.dict_im_algebraic_vars_syms[
                 gov.output_sig_syms[ 1 ] ]
    
        ac_gov[ω_idx, τm_sym_idx] = ωs /(2 * H)

     end
        
    return ac_gov 

end





"""

Gen state vars coupling to state of avr

"""
function SM_Ax_gen_S_avr_S( gen, avr )
        
    ac_avr = zeros(length(gen.state_vars_syms),
               length(avr.state_vars_syms) )
    
    eq_dash_idx = gen.dict_state_syms[:eq_dash]
    T_d_dash    = gen.T_d_dash

     
    if avr.output_sig_syms[ 1 ] ∈ avr.state_vars_syms
        
        vf_sym_idx  = avr.dict_state_syms[
            avr.output_sig_syms[ 1 ] ]

        ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash
                
    end     

    return ac_avr

end


"""

Gen state vars coupling to im_algebraic_vars of avr

"""
function SM_Ax_gen_S_avr_A( gen, avr )
        
    ac_avr = zeros(length(gen.state_vars_syms),
               length(avr.im_algebraic_vars_syms) )
    
    eq_dash_idx = gen.dict_state_syms[:eq_dash]
    T_d_dash    = gen.T_d_dash

    if avr.output_sig_syms[ 1 ] ∈ avr.im_algebraic_vars_syms

        vf_sym_idx  = avr.dict_im_algebraic_vars_syms[
            avr.output_sig_syms[ 1 ] ]

        ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash

    end
     
    return ac_avr

end



#-----------------------------------------
#-----------------------------------------





function Ax_gen_τm_vf( gen, D, H, Ωb, ωs, ra, xℓ,
                       X_d, X_q, X_d_dash, X_q_dash,
                       X_d_2dash, X_q_2dash,
                       T_d_dash, T_q_dash,
                       T_d_2dash, T_q_2dash )

    input_vector = [:τm,  :vf ]
    
    input_size = length( input_vector )

    Cx_τm_vf = zeros( length(gen.state_vars_syms),
                input_size )
    
    dict_input_syms_Idx::OrderedDict{Symbol, Int64} =
        OrderedDict(sym => ind
                    for (ind, sym) in enumerate( input_vector ) )
    
    τm_idx      = dict_input_syms_Idx[:τm]
    vf_idx      = dict_input_syms_Idx[:vf]
    
    δ_idx       = gen.dict_state_syms[:δ]
    ω_idx       = gen.dict_state_syms[:ω]
    ed_dash_idx = gen.dict_state_syms[:ed_dash]
    eq_dash_idx = gen.dict_state_syms[:eq_dash]

    Cx_τm_vf[δ_idx, τm_idx] = 0
    Cx_τm_vf[ω_idx, τm_idx] =  ωs /(2 * H)
    Cx_τm_vf[eq_dash_idx, vf_idx] = 1/T_d_dash
    
    # #                 τm          vf
    # Cx_τm_vf = [       0           0;
    #               ( ωs /(2 * H))   0;
    #                    0           0;
    #                    0       1/T_d_dash
    #       ]

    return Cx_τm_vf
    
end


function Ax_gen_τm_vf( D, H, Ωb, ωs, ra, xℓ,
                       X_d, X_q, X_d_dash, X_q_dash,
                       X_d_2dash, X_q_2dash,
                       T_d_dash, T_q_dash,
                       T_d_2dash, T_q_2dash )

    #            τm             vf
    Cx_τm_vf = [0               0;
              (ωs/(2 * H))      0;
                0               0;
                0              1/T_d_dash
          ]

    return Cx_τm_vf 
    
end



# ------------------------------------------------------
# Synchronous Generators nodes matrix
# ------------------------------------------------------


function Ax_SC_gen(gen, D, H, Ωb, ωs, ra, xℓ,
                   X_d, X_q, X_d_dash, X_q_dash,
                   X_d_2dash, X_q_2dash,
                   T_d_dash, T_q_dash,
                   T_d_2dash, T_q_2dash )

    ωs = gen.ωs
    
    Ax = zeros( length(gen.state_vars_syms),
                length(gen.state_vars_syms) )

    δ_idx       = gen.dict_state_syms[:δ]
    ω_idx       = gen.dict_state_syms[:ω]
    ed_dash_idx = gen.dict_state_syms[:ed_dash]
    eq_dash_idx = gen.dict_state_syms[:eq_dash]

    Ax[δ_idx, ω_idx]             = 1
    Ax[ω_idx, ω_idx]             = (-D * ωs /(2 * H))
    Ax[ed_dash_idx, ed_dash_idx] = -1/T_q_dash
    Ax[eq_dash_idx, eq_dash_idx] = -1/T_d_dash    

    # #     δ        ω         ed_dash     eq_dash
    # Ax = [0        1           0            0;
    #       0 (-D*ωs/(2 * H))    0            0;
    #       0        0          -1/T_q_dash   0;
    #       0        0           0            -1/T_d_dash
    #       ]

    return Ax 
    
end


function Ax_SC_gen(D, H, Ωb, ωs, ra, xℓ,
                   X_d, X_q, X_d_dash, X_q_dash,
                   X_d_2dash, X_q_2dash,
                   T_d_dash, T_q_dash,
                   T_d_2dash, T_q_2dash )

    #     δ        ω          ed_dash     eq_dash
    Ax = [0        1           0            0;
          0 (-D * ωs/(2 * H))  0            0;
          0        0           (-1/T_q_dash)  0;
          0        0           0           (-1/T_d_dash)
          ]

    return Ax 
    
end


function Bx_SC_gen(gen, D, H, Ωb, ωs, ra, xℓ,
                   X_d, X_q, X_d_dash, X_q_dash,
                   X_d_2dash, X_q_2dash,
                   T_d_dash, T_q_dash,
                   T_d_2dash, T_q_2dash )

    input_vector = [:id, :iq , :ph,  :vh]
    
    input_size = length( input_vector )

    Bx = zeros( length(gen.state_vars_syms),
                input_size )
    
    dict_input_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( input_vector ) )
    
    id_idx = dict_input_syms_Idx[:id]
    iq_idx = dict_input_syms_Idx[:iq]
    ph_idx = dict_input_syms_Idx[:ph]
    vh_idx = dict_input_syms_Idx[:vh]
    
    δ_idx       = gen.dict_state_syms[:δ]
    ω_idx       = gen.dict_state_syms[:ω]
    ed_dash_idx = gen.dict_state_syms[:ed_dash]
    eq_dash_idx = gen.dict_state_syms[:eq_dash]

    Bx[ed_dash_idx, iq_idx] = (X_q - X_q_dash)/T_q_dash
    Bx[eq_dash_idx, id_idx] = -(X_d - X_d_dash)/T_d_dash
  
    # [ω_idx, ph_idx] is zero for sync condeser
    
    # #     id      iq                           ph    vh
    # Bx = [ 0      0                            0     0;
    #        0      0                            0     0;
    #        0      (X_q - X_q_dash)/T_q_dash    0     0;
    #       -(X_d - X_d_dash)/T_d_dash  0         0     0
    #        ]
    
    return Bx 

end



function Bx_SC_gen( D, H, Ωb, ωs, ra, xℓ,
                    X_d, X_q, X_d_dash, X_q_dash,
                    X_d_2dash, X_q_2dash,
                    T_d_dash, T_q_dash,
                    T_d_2dash, T_q_2dash )

    #     id      iq                           ph    vh
    
    Bx = [ 0      0                            0     0;
           0      0                            0     0;
           0      (X_q - X_q_dash)/T_q_dash    0     0;
          (X_d - X_d_dash)/T_d_dash  0         0     0
           ]
    
    return Bx 

end


function Cx_SC_gen( gen, D, H, Ωb, ωs, ra, xℓ,
                    X_d, X_q, X_d_dash, X_q_dash,
                    X_d_2dash, X_q_2dash,
                    T_d_dash, T_q_dash,
                    T_d_2dash, T_q_2dash )

    ωs = gen.ωs

    input_vector = [:ωs, :ω_ref,  :v_ref, :p_order0 ]
    
    input_size = length( input_vector )

    Cx = zeros( length(gen.state_vars_syms),
                input_size )
    
    dict_input_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( input_vector ) )
    
    ωs_idx       = dict_input_syms_Idx[:ωs]
    ω_ref_idx    = dict_input_syms_Idx[:ω_ref]
    v_ref_idx    = dict_input_syms_Idx[:v_ref]
    p_order0_idx = dict_input_syms_Idx[:p_order0]

    
    δ_idx       = gen.dict_state_syms[:δ]
    gen_ω_idx   = gen.dict_state_syms[:ω]
    ed_dash_idx = gen.dict_state_syms[:ed_dash]
    eq_dash_idx = gen.dict_state_syms[:eq_dash]

    Cx[δ_idx,       ωs_idx]    = -1
    Cx[gen_ω_idx,   ωs_idx]    =  D * ωs /(2*H)
    
    # #     ωs       ω_ref        v_ref       p_order0
    
    # Cx = [-1          0          0             0;
    #       ωs*D/(2*H)  0          0             0;
    #       0           0          0             0;
    #       0           0          0             0]
    
    return Cx 

end



function Cx_SC_gen( D, H, Ωb, ωs, ra, xℓ,
                    X_d, X_q, X_d_dash, X_q_dash,
                    X_d_2dash, X_q_2dash,
                    T_d_dash, T_q_dash,
                    T_d_2dash, T_q_2dash )

    # #    ωs          ω_ref        v_ref       p_order0
    
    Cx = [-1             0       0            0;
          D * ωs/(2*H)   0       0            0;
          0              0       0            0;
          0              0       0            0]
    
    return Cx 

end



function Ac_SC_gen_gov( gen, gov )


    """
    This function should never be called
    for synch generator there is no governor, hence
    ac_gov = 0
    """    
    
    return nothing

end

# ------------------------------------------------------

function Ac_SC_gen_avr( gen, avr )
        
    ac_avr = zeros(length(gen.state_vars_syms),
               length(avr.state_vars_syms) )
    
    eq_dash_idx = gen.dict_state_syms[:eq_dash]
    T_d_dash    = gen.T_d_dash
    

   """
    In some cases, the output may not be a pure state var
    """
    if avr.output_sig_syms[ 1 ] ∈ avr.state_vars_syms
        
        vf_sym_idx  = avr.dict_state_syms[
            avr.output_sig_syms[ 1 ]]
        
        ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash
        
    else
        for (a_state, a_state_coeff ) ∈ zip(
            avr.state_vars_syms, avr.output_sig_state_coeff)
            
            a_state_idx = avr.dict_state_syms[ a_state  ]
            
            ac_avr[eq_dash_idx, a_state_idx] =
                a_state_coeff * 1/T_d_dash
        end
        
    end
     
    return ac_avr

end



function Ac_SC_S_gen_avr( gen, avr )
        
    ac_avr = zeros(length(gen.state_vars_syms),
               length(avr.state_vars_syms) )
    
    eq_dash_idx = gen.dict_state_syms[:eq_dash]
    T_d_dash    = gen.T_d_dash

    
    if avr.output_sig_syms[ 1 ] ∈ avr.state_vars_syms
        
        vf_sym_idx  = avr.dict_im_state_vars_syms[
            avr.output_sig_syms[ 1 ] ]

        ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash
                
    end
     

    return ac_avr

end


function Ac_SC_A_gen_avr( gen, avr )
        
    ac_avr = zeros(length(gen.state_vars_syms),
               length(avr.im_algebraic_vars_syms) )
    
    eq_dash_idx = gen.dict_state_syms[:eq_dash]
    T_d_dash    = gen.T_d_dash

    vf_sym_idx  = avr.dict_im_algebraic_vars_syms[
        avr.output_sig_syms[ 1 ] ]
    
    ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash
    
    # if avr.output_sig_syms[ 1 ] ∈ avr.im_algebraic_vars_syms  
    #     vf_sym_idx  = avr.dict_im_algebraic_vars_syms[
    #         avr.output_sig_syms[ 1 ] ]
    #     ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash    
    # end
     

    return ac_avr

end


#-----------------------------------------
#-----------------------------------------



"""
Gen state vars coupling to state vars of gov
"""
function SC_Ax_gen_S_gov_S( gen, gov )

    nothing
    
end


"""

Gen state vars coupling to im_algberaic_vars of gov

"""
function SC_Ax_gen_S_gov_A( gen, gov )

    nothing

end



"""

Gen state vars coupling to state of avr

"""
function SC_Ax_gen_S_avr_S( gen, avr )
        
    ac_avr = zeros(length(gen.state_vars_syms),
               length(avr.state_vars_syms) )
    
    eq_dash_idx = gen.dict_state_syms[:eq_dash]
    T_d_dash    = gen.T_d_dash

     
    if avr.output_sig_syms[ 1 ] ∈ avr.state_vars_syms
        
        vf_sym_idx  = avr.dict_state_syms[
            avr.output_sig_syms[ 1 ] ]

        ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash
                
    end     

    return ac_avr

end


"""

Gen state vars coupling to im_algebraic_vars of avr

"""
function SC_Ax_gen_S_avr_A( gen, avr )
        
    ac_avr = zeros(length(gen.state_vars_syms),
               length(avr.im_algebraic_vars_syms) )
    
    eq_dash_idx = gen.dict_state_syms[:eq_dash]
    T_d_dash    = gen.T_d_dash

    if avr.output_sig_syms[ 1 ] ∈ avr.im_algebraic_vars_syms

        vf_sym_idx  = avr.dict_im_algebraic_vars_syms[
            avr.output_sig_syms[ 1 ] ]

        ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash

    end
     
    return ac_avr

end



#-----------------------------------------
#-----------------------------------------


function Ax_SC_gen_τm_vf( gen, D, H, Ωb, ωs, ra, xℓ,
                          X_d, X_q, X_d_dash, X_q_dash,
                          X_d_2dash, X_q_2dash,
                          T_d_dash, T_q_dash,
                          T_d_2dash, T_q_2dash )


    input_vector = [:τm,  :vf ]
    
    input_size = length( input_vector )

    Cx_τm_vf = zeros( length(gen.state_vars_syms),
                input_size )
    
    dict_input_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( input_vector ) )
    
    τm_idx      = dict_input_syms_Idx[:τm]
    vf_idx      = dict_input_syms_Idx[:vf]
    
    δ_idx       = gen.dict_state_syms[:δ]
    ω_idx       = gen.dict_state_syms[:ω]
    ed_dash_idx = gen.dict_state_syms[:ed_dash]
    eq_dash_idx = gen.dict_state_syms[:eq_dash]

    Cx_τm_vf[δ_idx, τm_idx] = 0  # 1
    Cx_τm_vf[ω_idx, τm_idx] = 0  # ωs /(2 * H)
    Cx_τm_vf[eq_dash_idx, vf_idx] = 1/T_d_dash
  

    # #                 τm    vf
    # Cx_τm_vf = [       0      0;
    #                    0      0;
    #                    0      0;
    #                    0     1/T_d_dash
    #       ]
    
    return Cx_τm_vf
    
end


function Ax_SC_gen_τm_vf( D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash )


    #            τm        vf
    Cx_τm_vf = [0           0;
                0           0;
                0           0;
                0           1/T_d_dash
          ]
    
    return Cx_τm_vf 
    
end



# ------------------------------------------------------
# Generators nodes stability matrix
# ------------------------------------------------------


function stab_Ax_gen(gen, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash )

    ωs = gen.ωs
    
    Ax = zeros( length(gen.stab_state_vars_syms),
                  length(gen.stab_state_vars_syms) )

    ω_idx       = gen.dict_stab_state_syms[:ω]
    ed_dash_idx = gen.dict_stab_state_syms[:ed_dash]
    eq_dash_idx = gen.dict_stab_state_syms[:eq_dash]

    Ax[ω_idx, ω_idx]             = (-D * ωs /(2 * H))
    Ax[ed_dash_idx, ed_dash_idx] = -1/T_q_dash
    Ax[eq_dash_idx, eq_dash_idx] = -1/T_d_dash
    

    # #           ω     ed_dash     eq_dash
    # Ax =[
    #      (-D/(2 * H))  0            0;
    #             0      -1/T_q_dash  0;
    #             0      0           -1/T_d_dash
    #          ]

    return Ax 
    
end


function stab_Ax_gen(D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash )

    #           ω          ed_dash     eq_dash
    stab_Ax =[
         (-D * ωs/(2 * H))  0            0;
                0           (-1/T_q_dash)  0;
                0           0           (-1/T_d_dash)
            ]

    return stab_Ax 
    
end


function stab_Bx_gen(gen, D, H, Ωb, ωs, ra, xℓ, X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash, T_d_dash, T_q_dash, T_d_2dash, T_q_2dash )

    input_vector = [:id, :iq , :ph,  :vh]
    
    input_size = length( input_vector )

    stab_Bx = zeros( length(gen.stab_state_vars_syms),
                input_size )
    
    dict_input_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( input_vector ) )
    
    id_idx = dict_input_syms_Idx[:id]
    iq_idx = dict_input_syms_Idx[:iq]
    ph_idx = dict_input_syms_Idx[:ph]
    vh_idx = dict_input_syms_Idx[:vh]
    
    ω_idx       = gen.dict_stab_state_syms[:ω]
    ed_dash_idx = gen.dict_stab_state_syms[:ed_dash]
    eq_dash_idx = gen.dict_stab_state_syms[:eq_dash]

    stab_Bx[ω_idx, ph_idx]       = -1* ωs/(2*H)
    stab_Bx[ed_dash_idx, iq_idx] = (X_q - X_q_dash)/T_q_dash
    stab_Bx[eq_dash_idx, id_idx] = (X_d - X_d_dash)/T_d_dash
    
    # #     id      iq                           ph    vh
    # Bx = [ 
    #        0      0                        -1/(2*H)  0;
    #        0      (X_q - X_q_dash)/T_q_dash    0     0;
    #       (X_d - X_d_dash)/T_d_dash  0         0     0
    #        ]

    return stab_Bx 

end



function stab_Bx_gen(
    D, H, Ωb, ωs, ra, xℓ,
    X_d, X_q, X_d_dash, X_q_dash,
    X_d_2dash, X_q_2dash, T_d_dash, T_q_dash,
    T_d_2dash, T_q_2dash )

    #     id      iq                           ph         vh
 stab_Bx = [ 
           0      0                        (-1 * ωs/(2*H))  0;
           0      (X_q - X_q_dash)/T_q_dash    0          0;
          (X_d - X_d_dash)/T_d_dash  0         0          0
           ]

    return stab_Bx 

end


function stab_Cx_gen(
    gen,
    D, H, Ωb, ωs, ra, xℓ,
    X_d, X_q, X_d_dash, X_q_dash,
    X_d_2dash, X_q_2dash, T_d_dash, T_q_dash,
    T_d_2dash, T_q_2dash )

    ωs = gen.ωs

    input_vector = [:ω_ref, :τm,  :v_ref ]
    
    input_size = length( input_vector )

    stab_Cx = zeros( length(gen.stab_state_vars_syms),
                input_size )
    
    dict_input_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( input_vector ) )
    
    ω_ref_idx = dict_input_syms_Idx[:ω_ref]
    τm_idx    = dict_input_syms_Idx[:τm]
    v_ref_idx = dict_input_syms_Idx[:v_ref]
    
    ω_idx       = gen.dict_stab_state_syms[:ω]
    ed_dash_idx = gen.dict_stab_state_syms[:ed_dash]
    eq_dash_idx = gen.dict_stab_state_syms[:eq_dash]

    stab_Cx[ω_idx, ω_ref_idx]       = -D * ωs /(2*H)
    stab_Cx[ω_idx, τm_idx]          = 1 * ωs /(2*H)
    stab_Cx[eq_dash_idx, v_ref_idx] = 1/T_d_dash
    
    # #     ω           τm        v_ref
    # Cx = [
    #      -ωs*D/(2*H)  ωs/(2*H)   0;
    #       0           0          0;
    #       0           0          1/T_d_dash ]
    
    return stab_Cx 

end



function stab_Cx_gen(
    D, H, Ωb, ωs, ra, xℓ,
    X_d, X_q, X_d_dash, X_q_dash, X_d_2dash, X_q_2dash,
    T_d_dash, T_q_dash, T_d_2dash, T_q_2dash )

    #           ωs                 τm          v_ref
    stab_Cx = [
                 -D * ωs/(2*H)      1*ωs/(2*H)   0;
                 0                  0            0;
                 0                  0            1/T_d_dash]

    return stab_Cx 

end



function stab_Ac_gen_gov( gen, gov )
    
    stab_ac_gov = zeros(length(gen.stab_state_vars_syms),
               length(gov.state_vars_syms) )

    ωs    = gen.ωs
    
    ω_idx = gen.dict_stab_state_syms[:ω]
    
    H     = gen.H
    
    """
    In some cases, the output may not be a pure state var,
    e.g. governors
    """
    if gov.output_sig_syms[ 1 ] ∈ gov.state_vars_syms
        
        out_sym_idx = gov.dict_state_syms[gov.output_sig_syms[ 1 ] ]
        stab_ac_gov[ω_idx, out_sym_idx] = ωs /(2 * H)
        
    else
        for (a_state, a_state_coeff ) ∈ zip(gov.state_vars_syms, gov.output_sig_state_coeff)
            
            a_state_idx = gov.dict_state_syms[ a_state  ]
            
            stab_ac_gov[ω_idx, a_state_idx] = a_state_coeff * ωs /(2 * H)
        end
        
    end
    
    return stab_ac_gov 
    

end


function stab_Ac_gen_avr( gen, avr )
        
    stab_ac_avr = zeros(length(gen.stab_state_vars_syms),
               length(avr.state_vars_syms) )
    
    eq_dash_idx = gen.dict_stab_state_syms[:eq_dash]

    T_d_dash    = gen.T_d_dash    

   """
    In some cases, the output may not be a pure state var
    """
    if avr.output_sig_syms[ 1 ] ∈ avr.state_vars_syms
        
        vf_sym_idx  = avr.dict_state_syms[avr.output_sig_syms[ 1 ]]
        stab_ac_avr[eq_dash_idx, vf_sym_idx] = 1/T_d_dash
        
    else
        for (a_state, a_state_coeff ) ∈ zip(avr.state_vars_syms, avr.output_sig_state_coeff)
            
            a_state_idx = avr.dict_state_syms[ a_state  ]
            
            stab_ac_avr[eq_dash_idx, a_state_idx] = a_state_coeff * 1/T_d_dash
        end
        
    end
     
    return stab_ac_avr    

end


function stab_Ax_gen_τm_vf(
    gen,
    D, H, Ωb, ωs, ra, xℓ, X_d, X_q,
    X_d_dash, X_q_dash, X_d_2dash, X_q_2dash,
    T_d_dash, T_q_dash, T_d_2dash, T_q_2dash )


    input_vector = [:τm,  :vf ]
    
    input_size = length( input_vector )

    stab_Cx_τm_vf = zeros( length(gen.stab_state_vars_syms),
                input_size )
    
    dict_input_syms_Idx::OrderedDict{Symbol, Int64} = OrderedDict(sym => ind for (ind, sym) in enumerate( input_vector ) )
    
    τm_idx      = dict_input_syms_Idx[:τm]
    vf_idx      = dict_input_syms_Idx[:vf]
    
    ω_idx       = gen.dict_stab_state_syms[:ω]
    ed_dash_idx = gen.dict_stab_state_syms[:ed_dash]
    eq_dash_idx = gen.dict_stab_state_syms[:eq_dash]

    stab_Cx_τm_vf[ω_idx, τm_idx] = ωs /(2 * H)
    stab_Cx_τm_vf[eq_dash_idx, vf_idx] = 1/T_d_dash
    
    # #                 τm    vf
    # Cx_τm_vf = [
    #             (1/(2 * H))   0;
    #                    0      0;
    #                    0     1/T_d_dash
    #       ]

    return stab_Cx_τm_vf
    
end


function stab_Ax_gen_τm_vf(
    D, H, Ωb, ωs, ra, xℓ,
    X_d, X_q, X_d_dash, X_q_dash,
    X_d_2dash, X_q_2dash, T_d_dash, T_q_dash,
    T_d_2dash, T_q_2dash )

    #                 τm        vf
    stab_Cx_τm_vf = [
                     (ωs/(2 * H))  0;
                       0           0;
                       0           1/T_d_dash
          ]

    return stab_Cx_τm_vf
    
end

#--------------------------------------------------------
# Generation plants system matrix
#--------------------------------------------------------


function only_gen_in_plant(plant )

    return (:Gen ∈ propertynames(plant)) &&  (:Gov ∉ propertynames(plant)) &&  (:Exc ∉ propertynames(plant) )

end


function only_gen_in_plant(plant; only_gen = false )

    return (:Gen ∈ propertynames(plant)) &&  (:Gov ∉ propertynames(plant)) &&  (:Exc ∉ propertynames(plant) ) || only_gen == true

end


function only_gen_avr_in_plant(plant)

    return (:Gen ∈ propertynames(plant)) &&  (:Exc ∈ propertynames(plant))  && (:Gov ∉ propertynames(plant) )  

end



function only_gen_gov_in_plant(plant)

    return (:Gen ∈ propertynames(plant)) &&  (:Gov ∈  propertynames(plant))  && (:Exc ∉  propertynames(plant))
    
end



function gen_gov_avr_in_plant(plant)

    (:Gen ∈ propertynames(plant)) &&  (:Gov ∈  propertynames(plant))  && (:Exc ∈  propertynames(plant))
end

# ------------------------------------------------------
# only gen
# ------------------------------------------------------


function get_only_gen_plant_system_matrices(
    Gen )

    return (Ax = Gen.Ax( Gen, Gen.param_matrix_values...),
            Bx = Gen.Bx( Gen, Gen.param_matrix_values...),
            Cx = Gen.Cx( Gen, Gen.param_matrix_values...),
            Ax_τm_vf = Gen.Ax_τm_vf(
                Gen, Gen.param_matrix_values...)) 

end


function get_only_gen_plant_system_matrices(
    plant;
    only_gen=true )

    (m_Ax,
     m_Bx,
     m_Cx,
     m_Ax_τm_vf) =
        get_only_gen_plant_system_matrices(
            plant.Gen )
    
    return (Ax = sparse( m_Ax ),
            Bx = sparse( m_Bx ),
            Cx = sparse( m_Cx ),
            Ax_τm_vf =  m_Ax_τm_vf ) 

end


function get_only_gen_plant_system_matrices!(
    (Ax_view, Bx_view, Cx_view),
    plant;
    only_gen = true,
    τm_vf_view = nothing  )

    m_Ax, m_Bx, m_Cx, m_Ax_τm_vf =
        get_only_gen_plant_system_matrices(
            plant.Gen;
            only_gen = only_gen  )

    Ax_view[:,:] .= m_Ax
    
    Bx_view[:,:] .=  m_Bx
    
    Cx_view[:,:] .= m_Cx
    
    τm_vf_view[:,:] .= m_Ax_τm_vf
    
    return nothing

end


# ------------------------------------------------------
# Only gen and gov
# ------------------------------------------------------


function get_only_gen_gov_plant_system_matrices(
    Gen,
    Gov;
    with_control_alg_vars =
        false )

    Ax_gen = Gen.Ax(
        Gen,
        Gen.param_matrix_values... )
    
    Bx_gen = Gen.Bx(
        Gen,
        Gen.param_matrix_values...)
    
    Cx_gen = Gen.Cx(
        Gen,
        Gen.param_matrix_values...)

    Bx_gov = Gov.Bx( Gov )
    
    Cx_gov = Gov.Cx( Gov )
    
    if with_control_alg_vars == false


        Ax_gen_gov = Gen.Ax_gen_gov( Gen, Gov )

        Ax_gov = Gov.Ax( Gov )

        Ac_gov = Gov.Ac( Gov )

        Ax = sparse(
            vcat(
                hcat( Ax_gen, Ax_gen_gov ),
                hcat( Ac_gov, Ax_gov )) )

        Bx = sparse( vcat(Bx_gen, Bx_gov) )

        Cx = sparse( vcat(Cx_gen, Cx_gov) )


        return (; Ax, Bx, Cx)
    else

        #---------------------------------------
        
        Ax_gov = Gov.Ax( Gov )
        Ac_gov = Gov.Ac( Gov )
        Bx_gov = Gov.Bx( Gov )
        Cx_gov = Gov.Cx( Gov )

        #---------------------------------------
        
        Ax_gen_S_gov_S = Gen.Ax_gen_S_gov_S(
            Gen, Gov )
        
        Ax_gen_S_gov_A = Gen.Ax_gen_S_gov_S(
            Gen, Gov )

        gov_Ac_S_with_A    = Gov.Ac_S_with_A(
            Gov; ω_ref0 = ωs )
        
        gov_Ac_A_with_S    = Gov.Ac_A_with_S(
            Gov; ω_ref0 = ωs)
        
        gov_Ax_A           = Gov.Ax_A(
            Gov; ω_ref0 = ωs)
        
        gov_Ac_A_with_Sgen = Gov.Ac_A_with_Sgen(
            Gov; ω_ref0 = ωs)

        gov_Bx_A = Gov.Bx_A(
            Gov; ω_ref0 = ωs )

        gov_Cx_A = Gov.Cx_A(
            Gov; ω_ref0 = ωs )
        
        Ax = sparse(
            vcat(
                hcat(Ax_gen,
                     Ax_gen_S_gov_S,
                     Ax_gen_S_gov_A ),
                hcat(Ac_gov,
                     Ax_gov,
                     gov_Ac_S_with_A ),
                hcat(gov_Ac_A_with_Sgen,
                     gov_Ac_A_with_S,
                     gov_Ax_A ) ) )
    
        Bx = sparse( vcat(Bx_gen, Bx_gov, gov_Bx_A) )

        Cx = sparse( vcat(Cx_gen, Cx_gov, gov_Cx_A) )    

        return (; Ax, Bx, Cx)        
        
    end
    
    
end



function get_only_gen_gov_plant_system_matrices(
    plant;
    with_control_alg_vars =
        false )

    Ax_gen = plant.Gen.Ax(
        plant.Gen,
        plant.Gen.param_matrix_values... )
    
    Bx_gen = plant.Gen.Bx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    Cx_gen = plant.Gen.Cx(
        plant.Gen,
        plant.Gen.param_matrix_values...)

    Bx_gov = plant.Gov.Bx( plant.Gov )
    
    Cx_gov = plant.Gov.Cx( plant.Gov )
    
    if with_control_alg_vars == false 

        Ax_gen_gov = plant.Gen.Ax_gen_gov(
            plant.Gen, plant.Gov )

        Ax_gov = plant.Gov.Ax( plant.Gov )

        Ac_gov = plant.Gov.Ac( plant.Gov )


        Ax = sparse(
            vcat(
                hcat(Ax_gen, Ax_gen_gov),
                hcat(Ac_gov,Ax_gov)) )
    
        Bx = sparse( vcat(Bx_gen, Bx_gov) )

        Cx = sparse( vcat(Cx_gen, Cx_gov) )    

        return (; Ax, Bx, Cx)
    else
        
        #---------------------------------------
        
        Ax_gov = plant.Gov.Ax( plant.Gov )
        Ac_gov = plant.Gov.Ac( plant.Gov )
        Bx_gov = plant.Gov.Bx( plant.Gov )
        Cx_gov = plant.Gov.Cx( plant.Gov )


        #---------------------------------------
        
        Ax_gen_S_gov_S = plant.Gen.Ax_gen_S_gov_S(
            plant.Gen, plant.Gov )
        
        Ax_gen_S_gov_A = plant.Gen.Ax_gen_S_gov_S(
            plant.Gen, plant.Gov )

        gov_Ac_S_with_A    = plant.Gov.Ac_S_with_A(
            plant.Gov; ω_ref0 = ωs )
        
        gov_Ac_A_with_S    = plant.Gov.Ac_A_with_S(
            plant.Gov; ω_ref0 = ωs)
        
        gov_Ax_A           = plant.Gov.Ax_A(
            plant.Gov; ω_ref0 = ωs)
        
        gov_Ac_A_with_Sgen = plant.Gov.Ac_A_with_Sgen(
            plant.Gov; ω_ref0 = ωs)

        gov_Bx_A = plant.Gov.Bx_A(
            plant.Gov; ω_ref0 = ωs )

        gov_Cx_A = plant.Gov.Cx_A(
            plant.Gov; ω_ref0 = ωs )
        
        Ax = sparse(
            vcat(
                hcat(Ax_gen,
                     Ax_gen_S_gov_S,
                     Ax_gen_S_gov_A ),
                hcat(Ac_gov,
                     Ax_gov,
                     gov_Ac_S_with_A ),
                hcat(gov_Ac_A_with_Sgen,
                     gov_Ac_A_with_S,
                     gov_Ax_A ) ) )
    
        Bx = sparse( vcat(Bx_gen, Bx_gov, gov_Bx_A) )

        Cx = sparse( vcat(Cx_gen, Cx_gov, gov_Cx_A) )    

        return (; Ax, Bx, Cx)        
        
    end
            
end


function get_only_gen_gov_plant_system_matrices!(
    ( Ax_view, Bx_view, Cx_view ),
    plant
    ; with_control_alg_vars =
        false )

    Ax_gen = plant.Gen.Ax(
        plant.Gen,
        plant.Gen.param_matrix_values... )

    Bx_gen = plant.Gen.Bx(
        plant.Gen,
        plant.Gen.param_matrix_values...)

    Cx_gen = plant.Gen.Cx(
        plant.Gen,
        plant.Gen.param_matrix_values...)

    Bx_gov = plant.Gov.Bx( plant.Gov )

    Cx_gov = plant.Gov.Cx( plant.Gov )

    if with_control_alg_vars == false

        Ax_gen_gov = plant.Gen.Ax_gen_gov(
            plant.Gen, plant.Gov )

        Ax_gov = plant.Gov.Ax( plant.Gov )

        Ac_gov = plant.Gov.Ac( plant.Gov )

        Ax_view[:,:] .= vcat(
            hcat(Ax_gen, Ax_gen_gov),
            hcat(Ac_gov,Ax_gov))

        Bx_view[:,:] .= vcat(Bx_gen, Bx_gov)

        Cx_view[:,:] .= vcat(Cx_gen, Cx_gov)    
        
    else

        #---------------------------------------
        
        Ax_gov = plant.Gov.Ax( plant.Gov )
        Ac_gov = plant.Gov.Ac( plant.Gov )
        Bx_gov = plant.Gov.Bx( plant.Gov )
        Cx_gov = plant.Gov.Cx( plant.Gov )

        #---------------------------------------
        
        Ax_gen_S_gov_S = plant.Gen.Ax_gen_S_gov_S(
            plant.Gen, plant.Gov )
        
        Ax_gen_S_gov_A = plant.Gen.Ax_gen_S_gov_S(
            plant.Gen, plant.Gov )

        gov_Ac_S_with_A    = plant.Gov.Ac_S_with_A(
            plant.Gov; ω_ref0 = ωs )
        
        gov_Ac_A_with_S    = plant.Gov.Ac_A_with_S(
            plant.Gov; ω_ref0 = ωs)
        
        gov_Ax_A           = plant.Gov.Ax_A(
            plant.Gov; ω_ref0 = ωs)
        
        gov_Ac_A_with_Sgen = plant.Gov.Ac_A_with_Sgen(
            plant.Gov; ω_ref0 = ωs)

        gov_Bx_A = plant.Gov.Bx_A(
            plant.Gov; ω_ref0 = ωs )

        gov_Cx_A = plant.Gov.Cx_A(
            plant.Gov; ω_ref0 = ωs )

        Ax_view[:,:] .= vcat(
                hcat(Ax_gen,
                     Ax_gen_S_gov_S,
                     Ax_gen_S_gov_A ),
                hcat(Ac_gov,
                     Ax_gov,
                     gov_Ac_S_with_A ),
                hcat(gov_Ac_A_with_Sgen,
                     gov_Ac_A_with_S,
                     gov_Ax_A ) )
    
        Bx_view[:,:] .= vcat(Bx_gen, Bx_gov, gov_Bx_A) 

        Cx_view[:,:] .= vcat(Cx_gen, Cx_gov, gov_Cx_A)     
        
    end
    
    return nothing
    
end

# ------------------------------------------------------
# only gen and avr
# ------------------------------------------------------

function get_only_gen_avr_pure_plant_system_matrices(
    plant )

    Ax_gen = plant.Gen.Ax(
        plant.Gen,
        plant.Gen.param_matrix_values... )
    
    Bx_gen = plant.Gen.Bx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    Cx_gen = plant.Gen.Cx(
        plant.Gen,
        plant.Gen.param_matrix_values...)

    Ax_gen_S_avr_S = plant.Gen.Ax_gen_S_avr_S(
        plant.Gen, plant.Exc )

    avr_Ax_S_S = plant.Exc.Ax_S_S(
        plant.Exc,
        plant.dummy_vf_tilade,
        plant.dummy_vr1 )

    avr_Ax_S_Sgen = plant.Exc.Ax_S_Sgen( plant.Exc )

    avr_Bx_S = plant.Exc.Bx_S( plant.Exc )

    avr_Cx_S = plant.Exc.Cx_S( plant.Exc )

    Ax =
        vcat(
            hcat(Ax_gen,     Ax_gen_S_avr_S),
            hcat(avr_Ax_S_Sgen, avr_Ax_S_S)) 

    Bx =  vcat(Bx_gen, avr_Bx_S) 

    Cx =  vcat(Cx_gen, avr_Cx_S) 

    return (; Ax, Bx, Cx)
end


function get_only_gen_avr_pure_plant_system_matrices(
    Gen, Exc, dummy_vr1, dummy_vf_tilade )

    Ax_gen = Gen.Ax(
        Gen,
        Gen.param_matrix_values... )
    
    Bx_gen = Gen.Bx(
        Gen,
        Gen.param_matrix_values...)
    
    Cx_gen = Gen.Cx(
        Gen,
        Gen.param_matrix_values...)
    
    Ax_Sgen_avr_S = Gen.Ax_gen_avr(
        Gen, Exc )

    avr_Ax_S_S = Exc.Ax(
        Exc,
        dummy_vf_tilade,
        dummy_vr1 )

    avr_Ax_S_Sgen = Exc.Ac( Exc )

    avr_Bx_S = Exc.Bx( Exc )

    avr_Cx_S = Exc.Cx( Exc )

    Ax = vcat(
        hcat(Ax_gen, Ax_Sgen_avr_S),
        hcat(avr_Ax_S_Sgen, avr_Ax_S_S))

    Bx = vcat(Bx_gen, avr_Bx_S)

    Cx = vcat(Cx_gen, avr_Cx_S)

    return (; Ax, Bx, Cx)
    
end



function get_only_gen_avr_im_plant_system_matrices(
    plant; ref_in_state = false )

    Ax_gen = plant.Gen.Ax(
        plant.Gen,
        plant.Gen.param_matrix_values... )
    
    Bx_gen = plant.Gen.Bx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    Cx_gen = plant.Gen.Cx(
        plant.Gen,
        plant.Gen.param_matrix_values...)

    Ax_gen_S_avr_S = plant.Gen.Ax_gen_S_avr_S(
        plant.Gen,
        plant.Exc )

    Ax_gen_S_avr_A = plant.Gen.Ax_gen_S_avr_A(
        plant.Gen,
        plant.Exc )

    #---------------------------------------

    avr_Ax_S_S = plant.Exc.Ax_S_S(
        plant.Exc,
        plant.dummy_vf_tilade,
        plant.dummy_vr1)

    avr_Ax_S_Sgen = plant.Exc.Ax_S_Sgen( plant.Exc )

    avr_Bx_S = plant.Exc.Bx_S( plant.Exc )

    avr_Cx_S = plant.Exc.Cx_S( plant.Exc )

    #---------------------------------------

    avr_Ax_S_A = plant.Exc.Ax_S_A(
        plant.Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_S = plant.Exc.Ax_A_S(
        plant.Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_Sgen = plant.Exc.Ax_A_Sgen(
        plant.Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_A  = plant.Exc.Ax_A_A(
        plant.Exc;
        ref_in_state = ref_in_state )

    avr_Bx_A = plant.Exc.Bx_A(
        plant.Exc;
        ref_in_state = ref_in_state )

    avr_Cx_A = plant.Exc.Cx_A(
        plant.Exc;
        ref_in_state = ref_in_state )

    gen_m = hcat(Ax_gen, Ax_gen_S_avr_S, Ax_gen_S_avr_A )

    avr_S_m = hcat(avr_Ax_S_Sgen, avr_Ax_S_S, avr_Ax_S_A )

    avr_A_m =hcat(avr_Ax_A_Sgen, avr_Ax_A_S, avr_Ax_A_A )

    Ax = vcat(gen_m,
              avr_S_m,
              avr_A_m)              

    Bx = vcat(Bx_gen,
              avr_Bx_S,
              avr_Bx_A)

    Cx = vcat(Cx_gen,
              avr_Cx_S,
              avr_Cx_A )

    return (; Ax, Bx, Cx)                
            
end



function get_only_gen_avr_im_plant_system_matrices(
    Gen,
    Exc,
    dummy_vr1,
    dummy_vf_tilade;
    ref_in_state = false )

    Ax_gen = Gen.Ax(
        Gen,
        Gen.param_matrix_values... )
    
    Bx_gen = Gen.Bx(
        Gen,
        Gen.param_matrix_values...)
    
    Cx_gen = Gen.Cx(
        Gen,
        Gen.param_matrix_values...)

    Ax_gen_S_avr_S = Gen.Ax_gen_S_avr_S(
        Gen,
        Exc )

    Ax_gen_S_avr_A = Gen.Ax_gen_S_avr_A(
        Gen,
        Exc )

    #---------------------------------------

    avr_Ax_S_S = Exc.Ax_S_S(
        Exc,
        dummy_vf_tilade,
        dummy_vr1)

    avr_Ax_S_Sgen = Exc.Ax_S_Sgen( Exc )

    avr_Bx_S = Exc.Bx_S( Exc )

    avr_Cx_S = Exc.Cx_S( Exc )

    #---------------------------------------

    avr_Ax_S_A = Exc.Ax_S_A(
        Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_S = Exc.Ax_A_S(
        Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_Sgen = Exc.Ax_A_Sgen(
        Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_A  = Exc.Ax_A_A(
        Exc;
        ref_in_state = ref_in_state )

    avr_Bx_A = Exc.Bx_A(
        Exc;
        ref_in_state = ref_in_state )

    avr_Cx_A = Exc.Cx_A(
        Exc;
        ref_in_state = ref_in_state )

    gen_m = hcat(Ax_gen, Ax_gen_S_avr_S, Ax_gen_S_avr_A )

    avr_S_m = hcat(avr_Ax_S_Sgen, avr_Ax_S_S,  avr_Ax_S_A )

    avr_A_m =hcat(avr_Ax_A_Sgen, avr_Ax_A_S, avr_Ax_A_A )

    Ax = vcat(gen_m,
              avr_S_m,
              avr_A_m)              

    Bx = vcat(Bx_gen,
              avr_Bx_S,
              avr_Bx_A)

    Cx = vcat(Cx_gen,
              avr_Cx_S,
              avr_Cx_A )

    return (; Ax, Bx, Cx)                
    
end


function get_only_gen_avr_pure_or_im_plant_system_matrices(
    plant;
    with_control_alg_vars = false,
    ref_in_state = false )

    # if with_control_alg_vars == with_control_alg_vars
    
    if with_control_alg_vars == false
        
        return get_only_gen_avr_pure_plant_system_matrices(
            plant )
    else
        
        return get_only_gen_avr_im_plant_system_matrices(
            plant;
            ref_in_state =
                ref_in_state)
    end
      
end


function get_only_gen_avr_pure_or_im_plant_system_matrices(
    Gen, Exc, dummy_vr1, dummy_vf_tilade;
    with_control_alg_vars = false,
    ref_in_state = false )

    if with_control_alg_vars == false
        
        return get_only_gen_avr_pure_plant_system_matrices(
            Gen, Exc, dummy_vr1, dummy_vf_tilade )
    else
        return get_only_gen_avr_im_plant_system_matrices(
            Gen,
            Exc,
            dummy_vr1,
            dummy_vf_tilade;
            ref_in_state =
                ref_in_state )
    end
    

end


function get_only_gen_avr_plant_system_matrices(
    plant;
    with_control_alg_vars = false,
    ref_in_state = false )

    v_Ax, v_Bx, v_Cx =
        get_only_gen_avr_pure_or_im_plant_system_matrices(
        plant;
            with_control_alg_vars =
                with_control_alg_vars,
            ref_in_state =
                ref_in_state )

    Ax = sparse( v_Ax )                

    Bx = sparse( v_Bx )

    Cx = sparse( v_Cx )

    return (; Ax, Bx, Cx )                    
end


function get_only_gen_avr_plant_system_matrices!(
    ( Ax_view, Bx_view, Cx_view ),
    plant;
    with_control_alg_vars = false,
    ref_in_state = false  )

    v_Ax, v_Bx, v_Cx   =
        get_only_gen_avr_pure_or_im_plant_system_matrices(
            plant;
            with_control_alg_vars = with_control_alg_vars,
            ref_in_state = ref_in_state )
    
    Ax_view[:,:] .= v_Ax 

    Bx_view[:,:] .= v_Bx

    Cx_view[:,:] .= v_Cx         

    return nothing
    
    
end


function get_only_gen_avr_plant_system_matrices(
    Gen, Exc, dummy_vr1, dummy_vf_tilade;
    with_control_alg_vars = false,
    ref_in_state = false )

    v_Ax, v_Bx, v_Cx =
        get_only_gen_avr_pure_or_im_plant_system_matrices(
            plant;
            with_control_alg_vars = with_control_alg_vars,
            ref_in_state = ref_in_state )
    
    Ax = sparse( v_Ax  )

    Bx = sparse( v_Bx )

    Cx = sparse( v_Cx )

    return (; Ax, Bx, Cx)
    
end


# ------------------------------------------------------
# gen, gov and avr
# ------------------------------------------------------


function get_gen_gov_avr_pure_plant_system_matrices(
    Gen, Gov, Exc, dummy_vr1, dummy_vf_tilade;
    ω_ref0 = ωs, ref_in_state = false )

    #--------------------------------
        
    Ax_gen = Gen.Ax(
        Gen,
        Gen.param_matrix_values...)
    
    Bx_gen = Gen.Bx(
        Gen,
        Gen.param_matrix_values...)
    
    Cx_gen = Gen.Cx(
        Gen,
        Gen.param_matrix_values...)
        
    Ax_gen_S_gov_S = Gen.Ax_gen_S_gov_S(
        Gen, Gov)

    Ax_gen_S_avr_S = Gen.Ax_gen_S_avr_S(
        Gen, Exc)
    
    #--------------------------------

    gov_Ax_S_S = Gov.Ax_S_S(Gov )

    gov_Ax_S_Sgen = Gov.Ax_S_Sgen(Gov )

    gov_Bx_S = Gov.Bx_S(Gov )

    gov_Cx_S = Gov.Cx_S(Gov )

    #--------------------------------

    avr_Ax_S_S = Exc.Ax_S_S(
        Exc,
        dummy_vf_tilade,
        dummy_vr1)

    avr_Ax_S_Sgen = Exc.Ax_S_Sgen( Exc )

    avr_Bx_S = Exc.Bx_S( Exc )

    avr_Cx_S = Exc.Cx_S( Exc )

    #--------------------------------
    #--------------------------------

    z_gov_S_avr_S =
        zero_coupling_state_vars_to_state_vars(
            Gov, Exc )
    
    z_avr_S_gov_S =
        zero_coupling_state_vars_to_state_vars(
            Exc, Gov )
    
    #--------------------------------

    gen_m = hcat(Ax_gen, Ax_gen_S_gov_S, Ax_gen_S_avr_S)

    gov_m = hcat(gov_Ax_S_Sgen, gov_Ax_S_S, z_gov_S_avr_S )

    avr_m = hcat(avr_Ax_S_Sgen, z_avr_S_gov_S , avr_Ax_S_S)

    #--------------------------------

    Ax = vcat( gen_m,
               gov_m,
               avr_m )

    Bx = vcat(Bx_gen,
              gov_Bx_S,
              avr_Bx_S)

    Cx = vcat(Cx_gen,
              gov_Cx_S,
              avr_Cx_S )

    return (; Ax, Bx, Cx)        
end



function get_gen_gov_avr_im_plant_system_matrices(
    Gen, Gov, Exc, dummy_vr1, dummy_vf_tilade;
    ω_ref0 = ωs, ref_in_state = false )

    #--------------------------------
    
    Ax_gen = Gen.Ax(
        Gen,
        Gen.param_matrix_values...)
    
    Bx_gen = Gen.Bx(
        Gen,
        Gen.param_matrix_values...)
    
    Cx_gen = Gen.Cx(
        Gen,
        Gen.param_matrix_values...)

    #---------------------------------------

    Ax_gen_S_gov_S = Gen.Ax_gen_S_gov_S(
        Gen, Gov )

    Ax_gen_S_gov_A = Gen.Ax_gen_S_gov_A(
        Gen, Gov )

    #---------------------------------------
    
    Ax_gen_S_avr_S = Gen.Ax_gen_S_avr_S(
        Gen, Exc )

    Ax_gen_S_avr_A = Gen.Ax_gen_S_avr_A(
        Gen, Exc )

    #---------------------------------------

    gov_Ax_S_S    = Gov.Ax_S_S( Gov )

    gov_Ax_S_Sgen = Gov.Ax_S_Sgen( Gov )

    gov_Bx_S      = Gov.Bx_S( Gov )

    gov_Cx_S      = Gov.Cx_S( Gov )

    #--------------------------------

    avr_Ax_S_S = Exc.Ax_S_S(
        Exc,
        dummy_vf_tilade,
        dummy_vr1)

    avr_Ax_S_Sgen = Exc.Ax_S_Sgen( Exc )

    avr_Bx_S = Exc.Bx_S( Exc )

    avr_Cx_S = Exc.Cx_S( Exc )

    #--------------------------------

    gov_Ax_S_A = Gov.Ax_S_A(
        Gov; ω_ref0 = ωs,
        ref_in_state = ref_in_state  )

    gov_Ax_A_S = Gov.Ax_A_S(
        Gov; ω_ref0 = ωs,
        ref_in_state = ref_in_state)

    gov_Ax_A_A  = Gov.Ax_A_A(
        Gov; ω_ref0 = ωs,
        ref_in_state = ref_in_state)

    gov_Ax_A_Sgen = Gov.Ax_A_Sgen(
        Gov; ω_ref0 = ωs,
        ref_in_state = ref_in_state)

    #--------------------------------

    avr_Ax_S_A = Exc.Ax_S_A(
        Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_S = Exc.Ax_A_S(
        Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_A = Exc.Ax_A_A(
        Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_Sgen = Exc.Ax_A_Sgen(
        Exc;
        ref_in_state = ref_in_state )
    
    #--------------------------------

    z_gov_S_avr_A =
        zero_coupling_state_vars_to_im_alg_vars(
        Gov, Exc)

    z_gov_S_avr_S =
        zero_coupling_state_vars_to_state_vars(
        Gov, Exc )

    z_gov_A_avr_S =
        zero_coupling_im_alg_vars_to_state_vars(
        Gov, Exc )

    z_gov_A_avr_A =
        zero_coupling_im_alg_vars_to_im_alg_vars(
        Gov, Exc )

    #--------------------------------
    
    z_avr_S_gov_A =
        zero_coupling_state_vars_to_im_alg_vars(
         Exc, Gov)

    z_avr_S_gov_S =
        zero_coupling_state_vars_to_state_vars(
         Exc, Gov )


    z_avr_A_gov_S =
        zero_coupling_im_alg_vars_to_state_vars(
         Exc, Gov, )

    z_avr_A_gov_A =
        zero_coupling_im_alg_vars_to_im_alg_vars(
         Exc, Gov )

    #--------------------------------                

    gov_Bx_A = Gov.Bx_A(
        Gov; ω_ref0 = ωs,
        ref_in_state = ref_in_state )

    gov_Cx_A = Gov.Cx_A(
        Gov; ω_ref0 = ωs,
        ref_in_state = ref_in_state )

    #--------------------------------

    avr_Bx_A = Exc.Bx_A(
        Exc;
        ref_in_state = ref_in_state )

    avr_Cx_A = Exc.Cx_A(
        Exc;
        ref_in_state = ref_in_state )

    #--------------------------------

    gen_m = hcat(Ax_gen, Ax_gen_S_gov_S,
                 Ax_gen_S_avr_S, Ax_gen_S_gov_A,
                 Ax_gen_S_avr_A )

    gov_S_m = hcat( gov_Ax_S_Sgen, gov_Ax_S_S,
                    z_gov_S_avr_S, gov_Ax_S_A,
                    z_gov_S_avr_A )

    gov_A_m = hcat( gov_Ax_A_Sgen, gov_Ax_A_S,
                    z_gov_A_avr_S , gov_Ax_A_A,
                    z_gov_A_avr_A )

    avr_S_m = hcat( avr_Ax_S_Sgen, z_avr_S_gov_S,
                    avr_Ax_S_S, z_avr_S_gov_A,
                    avr_Ax_S_A )

    avr_A_m = hcat( avr_Ax_A_Sgen, z_avr_A_gov_S,
                    avr_Ax_A_S, z_avr_A_gov_A,
                    avr_Ax_A_A )

    #--------------------------------

    Ax = vcat(
        gen_m,
        gov_S_m,
        avr_S_m,
        gov_A_m,        
        avr_A_m)

   Bx = vcat(
       Bx_gen,
       gov_Bx_S,
       avr_Bx_S,
       gov_Bx_A,
       avr_Bx_A)

   Cx = vcat(
       Cx_gen,
       gov_Cx_S,
       avr_Cx_S,
       gov_Cx_A,
       avr_Cx_A)
    
    #--------------------------------

    return (; Ax, Bx, Cx)
    
end


function get_gen_gov_avr_pure_plant_system_matrices(
    plant; ω_ref0 = ωs, ref_in_state = false )

    #--------------------------------

    return get_gen_gov_avr_im_plant_system_matrices(
        plant.Gen,
        plant.Gov,
        plant.Exc,
        plant.dummy_vr1,
        plant.dummy_vf_tilade;
        ω_ref0 = ω_ref0 ,
        ref_in_state = ref_in_state )
    
end


function get_gen_gov_avr_im_plant_system_matrices(
    plant; ω_ref0 = ωs, ref_in_state = false )

        # plant.Exc,
        # plant.dummy_vf_tilade,
        # plant.dummy_vr1

    return  get_gen_gov_avr_im_plant_system_matrices(
        plant.Gen,
        plant.Gov,
        plant.Exc,
        plant.dummy_vr1,
        plant.dummy_vf_tilade;
        ω_ref0 = ω_ref0,
        ref_in_state = ref_in_state )
    
end



function get_only_gen_gov_avr_pure_or_im_plant_system_matrices(
    plant; with_control_alg_vars = false,
    ref_in_state = false  )

    if with_control_alg_vars == false
        
        return get_gen_gov_avr_pure_plant_system_matrices(
            plant )
        
    else
        
        return get_gen_gov_avr_im_plant_system_matrices(
            plant; ref_in_state = ref_in_state )
    end
    

end

# ref_in_state = ref_in_state
function get_only_gen_gov_avr_pure_or_im_plant_system_matrices(
    Gen, Gov, Exc, dummy_vr1, dummy_vf_tilade;
    with_control_alg_vars = false, ref_in_state = false )

    
    if with_control_alg_vars == false
        
        return get_gen_gov_avr_pure_plant_system_matrices(
            Gen, Gov, Exc, dummy_vr1, dummy_vf_tilade )
        
    else
        
        return get_gen_gov_avr_im_plant_system_matrices(
            Gen, Gov, Exc, dummy_vr1, dummy_vf_tilade;
            ref_in_state = ref_in_state )
        
    end

end



function get_gen_gov_avr_plant_system_matrices(
    plant; with_control_alg_vars = false,
    ref_in_state = false  )

    v_Ax, v_Bx, v_Cx =
        get_only_gen_gov_avr_pure_or_im_plant_system_matrices(
            plant;
            with_control_alg_vars =
                with_control_alg_vars,
            ref_in_state = ref_in_state  )
    
    #--------------------------------

    Ax = sparse( v_Ax )

    Bx = sparse( v_Bx )

    Cx = sparse( v_Cx )
    
    #--------------------------------
                
    return (; Ax, Bx, Cx )
        
end


function get_gen_gov_avr_plant_system_matrices!(
    ( Ax_view, Bx_view, Cx_view ),
    plant
    ; with_control_alg_vars =
        false,
    ref_in_state =
        false )

    #--------------------------------

    (v_Ax,
     v_Bx,
     v_Cx) =
        get_only_gen_gov_avr_pure_or_im_plant_system_matrices(
            plant;
            with_control_alg_vars =
                with_control_alg_vars,
            ref_in_state =
                ref_in_state  )
    
    #--------------------------------

    Ax_view[:,:] .= v_Ax

    Bx_view[:,:] .= v_Bx

    Cx_view[:,:] .= v_Cx

    return nothing
    
end


function get_gen_gov_avr_plant_system_matrices(
    Gen,
    Gov,
    Exc,
    dummy_vr1,
    dummy_vf_tilade
    ; with_control_alg_vars = false,
    ref_in_state = false )

    #--------------------------------

    (v_Ax,
     v_Bx,
     v_Cx) =
        get_only_gen_gov_avr_pure_or_im_plant_system_matrices(
            Gen,
            Gov,
            Exc,
            dummy_vr1,
            dummy_vf_tilade
            ;with_control_alg_vars = with_control_alg_vars,
            ref_in_state = ref_in_state )
        
    #--------------------------------

    Ax = sparse( v_Ax )

    Bx = sparse( v_Bx )

    Cx = sparse( v_Cx )

    #--------------------------------

     return (; Ax, Bx, Cx)

    
end


# ------------------------------------------------------
# ------------------------------------------------------

function get_a_plant_system_matrices(
    plant;
    only_gen = false,
    with_control_alg_vars = false  )


    if with_control_alg_vars == false
        
        if only_gen_in_plant(plant) || only_gen == true

            return  get_only_gen_plant_system_matrices(
                plant)

        elseif only_gen_avr_in_plant(plant)

            return get_only_gen_avr_plant_system_matrices(
                plant)


        elseif only_gen_gov_in_plant(plant)

            return get_only_gen_gov_plant_system_matrices(
                plant)

        elseif gen_gov_avr_in_plant(plant)

            return get_gen_gov_avr_plant_system_matrices(
                plant)

        else
            return nothing
        end
    else
        
        if only_gen_in_plant(plant) || only_gen == true

            return  get_only_gen_plant_system_matrices(
                plant; with_control_alg_vars = true )

        elseif only_gen_avr_in_plant(plant)

            return get_only_gen_avr_plant_system_matrices(
                plant ; with_control_alg_vars = true)


        elseif only_gen_gov_in_plant(plant)

            return get_only_gen_gov_plant_system_matrices(
                plant; with_control_alg_vars = true)

        elseif gen_gov_avr_in_plant(plant)

            return get_gen_gov_avr_plant_system_matrices(
                plant; with_control_alg_vars = true )

        else
            return nothing
        end        
    end
        
end



function get_a_plant_system_matrices!(
    ( Ax_view, Bx_view, Cx_view ), plant;
    only_gen = false,
    τm_vf_view = nothing,
    with_control_alg_vars = false )


    if with_control_alg_vars == false
        
        if only_gen_in_plant(plant) || only_gen == true

            get_only_gen_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ),
                plant; τm_vf_view = τm_vf_view  )

        elseif plant.isa_condenser || only_gen_avr_in_plant(
            plant)

            get_only_gen_avr_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ), plant )

        elseif only_gen_gov_in_plant(plant)

            get_only_gen_gov_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ), plant )

        elseif gen_gov_avr_in_plant(plant)

            get_gen_gov_avr_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ), plant )

        else
            nothing
        end
    else
        
        if only_gen_in_plant(plant) || only_gen == true

            get_only_gen_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ), plant;
                τm_vf_view = τm_vf_view,
                with_control_alg_vars = true  )

        elseif plant.isa_condenser || only_gen_avr_in_plant(
            plant)

            get_only_gen_avr_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ), plant;
                with_control_alg_vars = true )

        elseif only_gen_gov_in_plant(plant)

            get_only_gen_gov_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ), plant;
                with_control_alg_vars = true )

        elseif gen_gov_avr_in_plant(plant)

            get_gen_gov_avr_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ), plant;
                with_control_alg_vars = true )

        else
            nothing
        end        
    end
            
end


#--------------------------------------------------------


function instantiate_a_plant_system_matrices(
    plant; only_gen = true, with_control_alg_vars = false )

    if with_control_alg_vars == false
        
        if only_gen_in_plant( plant ) || only_gen == true

            return  get_only_gen_plant_system_matrices(
                plant.Gen )

        elseif only_gen_avr_in_plant( plant )

            return get_only_gen_avr_plant_system_matrices(
                plant.Gen,
                plant.Exc,
                plant.dummy_vr1,
                plant.dummy_vf_tilade  )

        elseif only_gen_gov_in_plant( plant )

            return get_only_gen_gov_plant_system_matrices(
                plant.Gen,
                plant.Gov )

        elseif gen_gov_avr_in_plant( plant )

            return get_gen_gov_avr_plant_system_matrices(
                plant.Gen, plant.Gov,
                plant.Exc,
                plant.dummy_vr1,
                plant.dummy_vf_tilade )

        else
            return nothing
        end
    else
        if only_gen_in_plant( plant ) || only_gen == true

            return  get_only_gen_plant_system_matrices(
                plant.Gen;
                with_control_alg_vars = true )

        elseif only_gen_avr_in_plant( plant )

            return get_only_gen_avr_plant_system_matrices(
                plant.Gen,
                plant.Exc,
                plant.dummy_vr1,
                plant.dummy_vf_tilade;
                with_control_alg_vars = true  )

        elseif only_gen_gov_in_plant( plant )

            return get_only_gen_gov_plant_system_matrices(
                plant.Gen,
                plant.Gov;
                with_control_alg_vars = true )

        elseif gen_gov_avr_in_plant( plant )

            return get_gen_gov_avr_plant_system_matrices(
                plant.Gen,
                plant.Gov,
                plant.Exc,
                plant.dummy_vr1,
                plant.dummy_vf_tilade;
                with_control_alg_vars = true )

        else
            return nothing
        end
        
    end
    
end



function instantiate_a_plant_system_matrices!(
    ( Ax_view, Bx_view, Cx_view ),
    plant;
    only_gen = false,
    with_control_alg_vars = false )


    if with_control_alg_vars == false
        
        if only_gen_in_plant( plant ) || only_gen == true

            get_only_gen_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ),
                plant.Gen )

        elseif only_gen_avr_in_plant( plant )

            get_only_gen_avr_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ),
                plant  )

        elseif only_gen_gov_in_plant( plant )

            get_only_gen_gov_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ),
                plant )

        elseif gen_gov_avr_in_plant( plant )

            return get_gen_gov_avr_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ),
                plant )

        else
           return  nothing
        end
    else
        if only_gen_in_plant( plant ) || only_gen == true

            get_only_gen_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ),
                plant.Gen;
                with_control_alg_vars = true )

        elseif only_gen_avr_in_plant( plant )

            get_only_gen_avr_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ),
                plant;
                with_control_alg_vars = true  )

        elseif only_gen_gov_in_plant( plant )

            get_only_gen_gov_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ),
                plant;
                with_control_alg_vars = true )

        elseif gen_gov_avr_in_plant( plant )

            return get_gen_gov_avr_plant_system_matrices!(
                ( Ax_view, Bx_view, Cx_view ),
                plant;
                with_control_alg_vars = true )

        else
           return  nothing
        end
        
    end
    
        
end



#--------------------------------------------------
###################################################
#--------------------------------------------------

#--------------------------------------------------
# im functions
#--------------------------------------------------



function create_flattened_type_im_vars_dims_offset_Idx_for_gens_nodes( comp_collection::OrderedDict; only_gen = false )

    gens_nodes_collection = get_gens_nodes(
        comp_collection )
    
    if only_gen == false
                
        return create_flattened_type_dims_offset_Idx(
            get_components_im_vars_syms,
            gens_nodes_collection)
    else
        sync_gens_collection = [
            a_plant.Gen
            for a_plant in gens_nodes_collection ]

        return create_flattened_type_dims_offset_Idx(
            get_components_state_vars,
            sync_gens_collection)
    end
    
        
end


         
function create_flattened_type_im_vars_dims_offset_Idx_for_gens_nodes(
    gens_nodes_collection; only_gen = false  )

    if only_gen == false

        return create_flattened_type_dims_offset_Idx(
            get_components_im_vars_syms,
            gens_nodes_collection )
        
    else
        
        sync_gens_collection = [
            a_plant.Gen
            for a_plant in gens_nodes_collection ]

        return create_flattened_type_dims_offset_Idx(
            get_components_state_vars,
            sync_gens_collection)

    end
            
end


function create_gens_nodes_im_aggregate_system_matrices_Idx_and_views(
    gens_nodes_collection )
    
    #--------------------------------------------------

    (nodes_dim_size,
     nodes_offset,
     nodes_state_Idx) = create_flattened_type_im_vars_dims_offset_Idx_for_gens_nodes( gens_nodes_collection )

    Ax_size =  sum(nodes_dim_size)
    
    #--------------------------------------------------

    no_nodes = length( gens_nodes_collection )
    
    id_iq_ph_vh_vector  = [:id, :iq , :ph,  :vh]

    id_iq_ph_vh_dims    = [
        length( id_iq_ph_vh_vector )
        for k in 1:no_nodes ]

    id_iq_ph_vh_offsets = create_offsets( id_iq_ph_vh_dims )
    
    id_iq_ph_vh_idxs    = create_idxs( id_iq_ph_vh_offsets,
                                    id_iq_ph_vh_dims )

    #--------------------------------------------------


    # ω_ref_τm_v_ref_vector
    ωs_ω_ref_v_ref_vector = [:ωs, :ω_ref0,
                             :v_ref0, :p_order0 ]

    ωs_ω_ref_v_ref_dims = [
        length( ωs_ω_ref_v_ref_vector )
        for k in 1:no_nodes ]

    ωs_ω_ref_v_ref_offsets = create_offsets(
        ωs_ω_ref_v_ref_dims)

    # ω_ref_τm_v_ref_idxs
    ωs_ω_ref_v_ref_p_order_idxs = create_idxs(
        ωs_ω_ref_v_ref_offsets,
        ωs_ω_ref_v_ref_dims )
        
    #--------------------------------------------------
    
    Bx_dims = [
        ( a_row, a_col )
        for ( a_row, a_col ) in zip(
            nodes_dim_size, id_iq_ph_vh_dims ) ]

    Cx_dims = [
        ( a_row, a_col )
        for ( a_row, a_col ) in zip(
            nodes_dim_size, ωs_ω_ref_v_ref_dims ) ]
    
    #--------------------------------------------------
    
    Bx_row_dims    = [ first( Bxi )
                       for Bxi in Bx_dims ]
    
    Bx_row_offsets = create_offsets(Bx_row_dims)
    
    Bx_row_idxs    = create_idxs(
        Bx_row_offsets, Bx_row_dims)

    Bx_col_dims    = [ last( Bxi )
                       for Bxi in Bx_dims ]
    
    Bx_col_offsets = create_offsets(Bx_col_dims)
    
    Bx_col_idxs    = create_idxs(
        Bx_col_offsets, Bx_col_dims)
    
    Bx_row_size    = sum( Bx_row_dims )
    
    Bx_col_size    = sum( Bx_col_dims )
    
    
    #--------------------------------------------------
    
    Cx_row_dims    = [ first( Cxi )
                       for Cxi in Cx_dims ]
    
    Cx_row_offsets = create_offsets(Cx_row_dims)
    
    Cx_row_idxs    = create_idxs(
        Cx_row_offsets, Cx_row_dims )

    Cx_col_dims    = [ last( Cxi )
                       for Cxi in Cx_dims ]
    
    Cx_col_offsets = create_offsets(Cx_col_dims)
    
    Cx_col_idxs    = create_idxs(
        Cx_col_offsets, Cx_col_dims )
    
    Cx_row_size    = sum( Cx_row_dims )
    
    Cx_col_size    = sum( Cx_col_dims )

    #--------------------------------------------------
      
    Ax_matrix = spzeros( Ax_size, Ax_size )
        
    Bx_matrix = spzeros( Bx_row_size, Bx_col_size )
    
    Cx_matrix = spzeros( Cx_row_size, Cx_col_size )

    #--------------------------------------------------

 
    vec_Ax_views = [
        view( Ax_matrix, node_state_Idx, node_state_Idx )
        for node_state_Idx in nodes_state_Idx ]

    vec_Bx_views = [
        view( Bx_matrix, Bx_row_idx, Bx_col_idx )
        for ( Bx_row_idx, Bx_col_idx ) in zip(
            Bx_row_idxs, Bx_col_idxs ) ]

    vec_Cx_views = [
        view( Cx_matrix, Cx_row_idx, Cx_col_idx )
        for ( Cx_row_idx, Cx_col_idx ) in zip(
            Cx_row_idxs, Cx_col_idxs ) ]
    

    #--------------------------------------------------
    
    #--------------------------------------------------

    Ax_Bx_Cx_views = (; vec_Ax_views,
                      vec_Bx_views,
                      vec_Cx_views )

    Ax_Bx_Cx_matrix = (; Ax_matrix,
                       Bx_matrix,
                       Cx_matrix )

    Bx_idxs = (; Bx_row_idxs,
               Bx_col_idxs )

    Cx_idxs = (; Cx_row_idxs,
               Cx_col_idxs )

    idxs = (; nodes_state_Idx, Bx_idxs, Cx_idxs,
            id_iq_ph_vh_idxs,
            ωs_ω_ref_v_ref_p_order_idxs  )


    return (; Ax_Bx_Cx_views,
            Ax_Bx_Cx_matrix,
            idxs )

        
end



function update_a_im_plant_Ax_system_matrices!(
    Ax, state, plant  )

    if only_gen_avr_in_plant(plant) || gen_gov_avr_in_plant(plant)
        dict_im_vars_syms = plant.dict_im_vars_syms

        vr1_idx  = dict_im_vars_syms[:vr1]
        
        vf_tilade_idx = dict_im_vars_syms[:vf_tilade]

        # (Ta, Te, Tf, Tr,
        #  Ka, Ke, Kf,
        #  V_R_max, V_R_min,
        #  Ae, Be = V_R_max) =
        #      plant.Exc.param_values

        (Ta, Te, Tf, Tr,
         Ka, Ke, Kf,
         V_R_max, V_R_min,
         Ae, Be) =
             plant.Exc.param_values

        # Be = V_R_max


        # It seems to me later that is not necessary to
        # update α, i.e the column of vr1 in Ax

        # vr1 = x[dict_state_vars_syms[:vr1]]
        # α = threshold_limits(vr1,V_R_max,V_R_min)/(vr1*Te)


        vf_tilade = state[ dict_im_vars_syms[:vf_tilade] ]


        α =  1/ Te

        γ = -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te

        # Ax[vr1_idx, vr1_idx]         = -α
        # Ax[vf_tilade_idx, vr1_idx]   = α
        
        Ax[vf_tilade_idx, vf_tilade_idx] = γ

    end

    return nothing        
end



function update_a_im_plant_system_matrices!(
    Ax, x, plant )

    dict_im_vars_syms = plant.dict_im_vars_syms

    vr1_idx       = dict_im_vars_syms[:vr1]
    vf_tilade_idx = dict_im_vars_syms[:vf_tilade]

    (Ta, Te, Tf, Tr, Ka, Ke, Kf,
     V_R_max, V_R_min,
     Ae, Be) =
        V_R_max = plant.Exc.param_values

    (Ta, Te,
     Tf, Tr,
     Ka, Ke, Kf,
     V_R_max, V_R_min,
     Ae, Be)  =
        plant.Exc.param_values

    # Be = V_R_max

    vf_tilade = x[ dict_im_vars_syms[:vf_tilade] ]

    # It seems to me later that is not necessary to
    # update α, i.e the column of vr1 in Ax

    # vr1 = x[dict_state_vars_syms[:vr1]]
    # α = threshold_limits(vr1, V_R_max, V_R_min)/(vr1 * Te)

    α =  1/ Te

    γ = -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te

    # Ax[vr1_idx, vr1_idx]             = -α
    # Ax[vf_tilade_idx, vr1_idx]       = α
    
    Ax[vf_tilade_idx, vf_tilade_idx] = γ

    return nothing        
end

  
function sd_update_a_im_plant_system_matrices!(
    Ax_sparse_nzvalues_DiffCache,
    sparse_row_idxs,
    sparse_col_idxs,
    a_stateDiffCache,
    plant )
    
    dict_im_vars_syms = plant.dict_im_vars_syms

    vr1_idx = dict_im_vars_syms[:vr1]
    
    vf_tilade_idx = dict_im_vars_syms[:vf_tilade]


    (Ta, Te, Tf, Tr,
     Ka, Ke, Kf,
     V_R_max, V_R_min,
     Ae, Be) =
         plant.Exc.param_values

    # Be = V_R_max

    vf_tilade = a_stateDiffCache[
        dict_im_vars_syms[:vf_tilade] ]

    # It seems to me later that is not necessary to
    # update α, i.e the column of vr1 in Ax

    # vr1  = a_gen_stateDiffCache[
    #     dict_state_vars_syms[:vr1]]
    
    # α    = threshold_limits(
    #         vr1, V_R_max, V_R_min)/(vr1 * Te)
    # α =  1/ Te

    # Ax[vr1_idx, vr1_idx]             = -α
    # Ax[vf_tilade_idx, vr1_idx]       = α
    
    γ = -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te

    vf_tilade_idx_in_sp_nzv =
        find_V_idx_in_sparse_matrix_IJV(
            vf_tilade_idx, vf_tilade_idx,
            sparse_row_idxs, sparse_col_idxs )
    
    Ax_sparse_nzvalues_DiffCache[vf_tilade_idx_in_sp_nzv] =
        γ

    return nothing
    
    # u_Ax = Matrix(Ax)
    # u_Ax[vf_tilade_idx, vf_tilade_idx]  = γ
    # return Matrix(Ax)
    
    
end



function update_a_im_plant_system_matrices_views!(
    (Ax_view, Bx_view, Cx_view), plant )

    plant.func_system_matrices[1](
        (Ax_view, Bx_view, Cx_view),
        plant;
        with_control_alg_vars = true )

    return nothing
    
end

"""
This function initialises system matrices views
` Ax_view, Bx_view, Cx_view` which are in vectors
`vec_Ax_views, vec_Bx_views, vec_Cx_views`.

These views are created by
`create_gens_nodes_im_aggregate_system_matrices_Idx_and_views`.

The intialisation is expected to be reflected in aggregate
system matrices `Ax_matrix, Bx_matrix, Cx_matrix `

"""
function init_im_plants_system_matrices_views!(
    (vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views),
    gens_nodes_collection )

    # update_a_plant_system_matrices_views!
    
    for (Ax_view, Bx_view, Cx_view, plant) in
        zip(vec_Ax_views, vec_Bx_views,
            vec_Cx_views, gens_nodes_collection)

        update_a_im_plant_system_matrices_views!(
            (Ax_view, Bx_view, Cx_view),
            plant )
    end
  
    return nothing        
end




function update_gens_nodes_im_Ax_system_matrices!(
    
    vec_Ax,
    state,
    nodes_state_Idx,
    gens_nodes_collection)

    # no_nodes = length( vec_Ax  )

    for (Ax_node_matrix, node_state_Idx, gen_node ) in zip(
        vec_Ax,
        nodes_state_Idx,
        gens_nodes_collection )

        update_a_im_plant_system_matrices!(
            Ax_node_matrix,
            state[node_state_Idx],
            gen_node )
    end    
    #--------------------------------------------------

    return nothing
end


function update_gens_nodes_im_system_matrices!(    
    vec_Ax,
    x,
    nodes_state_Idx,
    gens_nodes_collection )

    # no_nodes = length( vec_Ax  )

    for (Ax_node_matrix, node_state_Idx, gen_node ) in zip(
        vec_Ax,
        nodes_state_Idx,
        gens_nodes_collection )

        update_a_im_plant_system_matrices!(
            Ax_node_matrix,
            x[node_state_Idx],
            gen_node )
    end    
    #--------------------------------------------------

    return nothing
end

#--------------------------------------------------
###################################################
#--------------------------------------------------


#--------------------------------------------------------
# Generation plants system stability matrix
#--------------------------------------------------------


function get_only_gen_plant_system_stability_matrices(
    plant)

    return (
        stab_Ax = plant.Gen.stab_Ax(
            plant.Gen,
            plant.Gen.param_matrix_values...),
        stab_Bx = plant.Gen.stab_Bx(
            plant.Gen,
            plant.Gen.param_matrix_values...),
        stab_Cx = plant.Gen.stab_Cx(
            plant.Gen,
            plant.Gen.param_matrix_values...),
        stab_Ax_τm_vf = plant.Gen.stab_Ax_τm_vf(
            plant.Gen,
            plant.Gen.param_matrix_values...)) 

end


function get_only_gen_plant_system_stability_matrices(
    Gen; only_gen=true  )

    return (
        stab_Ax = Gen.stab_Ax(
            Gen,
            Gen.param_matrix_values...),
        stab_Bx = Gen.stab_Bx(
            Gen,
            Gen.param_matrix_values...),
        stab_Cx = Gen.stab_Cx(
            Gen,
            Gen.param_matrix_values...),
        stab_Ax_τm_vf = Gen.stab_Ax_τm_vf(
            Gen,
            Gen.param_matrix_values...)) 

end


function get_only_gen_plant_system_stability_matrices!(
    (stab_Ax_view, stab_Bx_view, stab_Cx_view),
    plant;
    only_gen=true,
    stab_τm_vf_view = nothing  )

    stab_Ax_view[:,:] .= plant.Gen.stab_Ax(
        plant.Gen,
        plant.Gen.param_matrix_values... )
    
    stab_Bx_view[:,:] .= plant.Gen.stab_Bx(
        plant.Gen,
        plant.Gen.param_matrix_values... )
    
    stab_Cx_view[:,:] .= plant.Gen.stab_Cx(
        plant.Gen,
        plant.Gen.param_matrix_values... )
    
    stab_τm_vf_view[:,:] .= plant.Gen.stab_Ax_τm_vf(
        plant.Gen,
        plant.Gen.param_matrix_values... )
    
    return nothing

end


function get_only_gen_gov_plant_system_stability_matrices(
    plant)

    stab_Ax_gen = plant.Gen.stab_Ax(
        plant.Gen,
        plant.Gen.param_matrix_values... )
    
    stab_Bx_gen = plant.Gen.stab_Bx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Cx_gen = plant.Gen.stab_Cx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Ax_gen_gov = plant.Gen.stab_Ax_gen_gov(
        plant.Gen,
        plant.Gov )

    Ax_gov = plant.Gov.Ax( plant.Gov )

    stab_Ac_gov = plant.Gov.stab_Ac( plant.Gov )

    Bx_gov = plant.Gov.Bx( plant.Gov )

    Cx_gov = plant.Gov.Cx( plant.Gov )

    stab_Ax = vcat(
        hcat(stab_Ax_gen, stab_Ax_gen_gov),
        hcat(stab_Ac_gov, Ax_gov))
    
    stab_Bx = vcat(stab_Bx_gen, Bx_gov)
    stab_Cx = vcat(stab_Cx_gen, Cx_gov)

    return (; stab_Ax, stab_Bx, stab_Cx)
    
end



function get_only_gen_gov_plant_system_stability_matrices(
    Gen, Gov )

    stab_Ax_gen = Gen.stab_Ax(
        Gen, Gen.param_matrix_values... )
    
    stab_Bx_gen = Gen.stab_Bx(
        Gen, Gen.param_matrix_values...)
    
    stab_Cx_gen = Gen.stab_Cx(
        Gen, Gen.param_matrix_values...)
    
    stab_Ax_gen_gov = Gen.stab_Ax_gen_gov(
        Gen, Gov )

    Ax_gov = Gov.Ax( Gov )

    stab_Ac_gov = Gov.stab_Ac( Gov )

    Bx_gov = Gov.Bx( Gov )

    Cx_gov = Gov.Cx( Gov )

    stab_Ax = vcat(
        hcat(stab_Ax_gen, stab_Ax_gen_gov),
        hcat(stab_Ac_gov, Ax_gov))
    
    stab_Bx = vcat(stab_Bx_gen, Bx_gov)
    
    stab_Cx = vcat(stab_Cx_gen, Cx_gov)

    return (; stab_Ax, stab_Bx, stab_Cx)
    
end



function get_only_gen_gov_plant_system_stability_matrices!(
    ( stab_Ax_view,
      stab_Bx_view,
      stab_Cx_view ),
    plant)

    stab_Ax_gen = plant.Gen.stab_Ax(
        plant.Gen,
        plant.Gen.param_matrix_values... )
    
    stab_Bx_gen = plant.Gen.stab_Bx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Cx_gen = plant.Gen.stab_Cx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Ax_gen_gov = plant.Gen.stab_Ax_gen_gov(
        plant.Gen, plant.Gov )

    Ax_gov = plant.Gov.Ax( plant.Gov )

    stab_Ac_gov = plant.Gov.stab_Ac( plant.Gov )

    Bx_gov = plant.Gov.Bx( plant.Gov )

    Cx_gov = plant.Gov.Cx( plant.Gov )

    stab_Ax_view[:,:] .= vcat(
        hcat(stab_Ax_gen, stab_Ax_gen_gov),
        hcat(stab_Ac_gov, Ax_gov))
    
    stab_Bx_view[:,:] .= vcat(stab_Bx_gen, Bx_gov)
    
    stab_Cx_view[:,:] .= vcat(stab_Cx_gen, Cx_gov)

    return nothing
    
end



function get_only_gen_avr_plant_system_stability_matrices(
    plant)

    stab_Ax_gen = plant.Gen.stab_Ax(
        plant.Gen,
        plant.Gen.param_matrix_values... )
    
    stab_Bx_gen = plant.Gen.stab_Bx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Cx_gen = plant.Gen.stab_Cx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Ax_gen_avr = plant.Gen.stab_Ax_gen_avr(
        plant.Gen, plant.Exc )

    Ax_avr = plant.Exc.Ax(
        plant.Exc,
        plant.dummy_vf_tilade,
        plant.dummy_vr1 )
    
    stab_Ac_avr = plant.Exc.stab_Ac( plant.Exc )
    
    Bx_avr = plant.Exc.Bx( plant.Exc )
    
    Cx_avr = plant.Exc.Cx( plant.Exc )

    stab_Ax = vcat(
        hcat(stab_Ax_gen, stab_Ax_gen_avr),
        hcat(stab_Ac_avr,Ax_avr))
    
    stab_Bx = vcat(stab_Bx_gen, Bx_avr)
    stab_Cx = vcat(stab_Cx_gen, Cx_avr)

    return (; stab_Ax, stab_Bx, stab_Cx)
    
    
end


function get_only_gen_avr_plant_system_stability_matrices(
    Gen, Exc, dummy_vr1, dummy_vf_tilade )

    stab_Ax_gen = Gen.stab_Ax(
        Gen, Gen.param_matrix_values... )
    
    stab_Bx_gen = Gen.stab_Bx(
        Gen, Gen.param_matrix_values...)
    
    stab_Cx_gen = Gen.stab_Cx(
        Gen, Gen.param_matrix_values...)
    
    stab_Ax_gen_avr = Gen.stab_Ax_gen_avr( Gen, Exc )

    Ax_avr = Exc.Ax( Exc, dummy_vf_tilade, dummy_vr1 )
    
    stab_Ac_avr = Exc.stab_Ac( Exc )
    
    Bx_avr      = Exc.Bx( Exc )
    
    Cx_avr      = Exc.Cx( Exc )

    stab_Ax = vcat(
        hcat(stab_Ax_gen, stab_Ax_gen_avr),
        hcat(stab_Ac_avr, Ax_avr))
    
    stab_Bx = vcat(stab_Bx_gen, Bx_avr)
    
    stab_Cx = vcat(stab_Cx_gen, Cx_avr)

    return (; stab_Ax, stab_Bx, stab_Cx)
    
    
end

# 


function get_only_gen_avr_plant_system_stability_matrices!(
    ( stab_Ax_view, stab_Bx_view, stab_Cx_view ),
    plant)

    stab_Ax_gen = plant.Gen.stab_Ax(
        plant.Gen,
        plant.Gen.param_matrix_values... )
    
    stab_Bx_gen = plant.Gen.stab_Bx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Cx_gen = plant.Gen.stab_Cx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Ax_gen_avr = plant.Gen.stab_Ax_gen_avr(
        plant.Gen, plant.Exc )

    Ax_avr = plant.Exc.Ax(
        plant.Exc,
        plant.dummy_vf_tilade,
        plant.dummy_vr1 )
    
    stab_Ac_avr = plant.Exc.stab_Ac( plant.Exc )
    Bx_avr = plant.Exc.Bx( plant.Exc )
    Cx_avr = plant.Exc.Cx( plant.Exc )

    stab_Ax_view[:,:] .= vcat(
        hcat(stab_Ax_gen, stab_Ax_gen_avr),
        hcat(stab_Ac_avr, Ax_avr))
    
    stab_Bx_view[:,:] .= vcat(stab_Bx_gen, Bx_avr)
    
    stab_Cx_view[:,:] .= vcat(stab_Cx_gen, Cx_avr)

    return nothing
    
    
end


function get_gen_gov_avr_plant_system_stability_matrices(
    plant)

    stab_Ax_gen = plant.Gen.stab_Ax(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Bx_gen = plant.Gen.stab_Bx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Cx_gen = plant.Gen.stab_Cx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Ax_gen_gov = plant.Gen.stab_Ax_gen_gov(
        plant.Gen, plant.Gov)
    
    stab_Ax_gen_avr = plant.Gen.stab_Ax_gen_avr(
        plant.Gen, plant.Exc)

    Ax_gov = plant.Gov.Ax(plant.Gov )
    
    stab_Ac_gov = plant.Gov.stab_Ac(plant.Gov )
    
    Bx_gov = plant.Gov.Bx(plant.Gov )
    
    Cx_gov = plant.Gov.Cx(plant.Gov )

    Ax_avr = plant.Exc.Ax(
        plant.Exc, plant.dummy_vf_tilade,
        plant.dummy_vr1)
    
    stab_Ac_avr = plant.Exc.stab_Ac(
        plant.Exc )
    
    Bx_avr = plant.Exc.Bx( plant.Exc )
    
    Cx_avr = plant.Exc.Cx( plant.Exc )

    gen_size = length(
        plant.Gen.stab_state_vars_syms )
    
    gov_size = length(
        plant.Gov.state_vars_syms )
    
    avr_size = length(
        plant.Exc.state_vars_syms )
    
    stab_Ax = vcat(
        hcat(stab_Ax_gen, stab_Ax_gen_gov, stab_Ax_gen_avr),
        hcat(stab_Ac_gov,  Ax_gov, zeros(gov_size, avr_size)),
        hcat(stab_Ac_avr, zeros(avr_size, gov_size), Ax_avr) )
    
    stab_Bx = vcat(stab_Bx_gen, Bx_gov, Bx_avr)
    
    stab_Cx = vcat(stab_Cx_gen, Cx_gov, Cx_avr)

    return (; stab_Ax, stab_Bx, stab_Cx)
    
end


function get_gen_gov_avr_plant_system_stability_matrices(
    Gen, Gov, Exc, dummy_vr1, dummy_vf_tilade )

    stab_Ax_gen = Gen.stab_Ax(Gen, Gen.param_matrix_values...)
    stab_Bx_gen = Gen.stab_Bx(Gen, Gen.param_matrix_values...)
    stab_Cx_gen = Gen.stab_Cx(Gen, Gen.param_matrix_values...)
    
    stab_Ax_gen_gov = Gen.stab_Ax_gen_gov(Gen, Gov)
    stab_Ax_gen_avr = Gen.stab_Ax_gen_avr(Gen, Exc )

    Ax_gov = Gov.Ax(Gov )
    stab_Ac_gov = Gov.stab_Ac(Gov )
    Bx_gov = Gov.Bx(Gov )
    Cx_gov = Gov.Cx(Gov )

    Ax_avr = Exc.Ax(Exc, dummy_vf_tilade, dummy_vr1)
    
    stab_Ac_avr = Exc.stab_Ac( Exc )
    
    Bx_avr = Exc.Bx( Exc )
    
    Cx_avr = Exc.Cx( Exc )

    gen_size = length(Gen.stab_state_vars_syms)
    
    gov_size = length(Gov.state_vars_syms)
    
    avr_size = length(Exc.state_vars_syms)
    
    stab_Ax = vcat(
        hcat(stab_Ax_gen, stab_Ax_gen_gov, stab_Ax_gen_avr),
        hcat(stab_Ac_gov, Ax_gov, zeros(gov_size, avr_size)),
        hcat(stab_Ac_avr, zeros(avr_size, gov_size), Ax_avr) )
    
    stab_Bx = vcat(stab_Bx_gen, Bx_gov, Bx_avr)
    
    stab_Cx = vcat(stab_Cx_gen, Cx_gov, Cx_avr)

    return (; stab_Ax, stab_Bx, stab_Cx)
    
end

# ( Ax_view, Bx_view, Cx_view ),


function get_gen_gov_avr_plant_system_stability_matrices!(
    ( stab_Ax_view, stab_Bx_view, stab_Cx_view ),
    plant )

    stab_Ax_gen = plant.Gen.stab_Ax(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Bx_gen = plant.Gen.stab_Bx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Cx_gen = plant.Gen.stab_Cx(
        plant.Gen,
        plant.Gen.param_matrix_values...)
    
    stab_Ax_gen_gov = plant.Gen.stab_Ax_gen_gov(
        plant.Gen, plant.Gov)
    
    stab_Ax_gen_avr = plant.Gen.stab_Ax_gen_avr(
        plant.Gen, plant.Exc)

    Ax_gov = plant.Gov.Ax(plant.Gov )
    
    stab_Ac_gov = plant.Gov.stab_Ac(plant.Gov )
    
    Bx_gov = plant.Gov.Bx(plant.Gov )
    
    Cx_gov = plant.Gov.Cx(plant.Gov )

    Ax_avr = plant.Exc.Ax(plant.Exc,
                          plant.dummy_vf_tilade,
                          plant.dummy_vr1)
    
    stab_Ac_avr = plant.Exc.stab_Ac(
        plant.Exc )
    
    Bx_avr = plant.Exc.Bx(
        plant.Exc )
    
    Cx_avr = plant.Exc.Cx(
        plant.Exc )

    gen_size = length(
        plant.Gen.stab_state_vars_syms)
    
    gov_size = length(
        plant.Gov.state_vars_syms)
    
    avr_size = length(
        plant.Exc.state_vars_syms)
    
    stab_Ax_view[:,:] .= vcat(
        hcat(stab_Ax_gen, stab_Ax_gen_gov, stab_Ax_gen_avr),
        hcat(stab_Ac_gov, Ax_gov, zeros(gov_size, avr_size)),
        hcat(stab_Ac_avr, zeros(avr_size, gov_size), Ax_avr) )
    
    stab_Bx_view[:,:] .= vcat(stab_Bx_gen, Bx_gov, Bx_avr)
    
    stab_Cx_view[:,:] .= vcat(stab_Cx_gen, Cx_gov, Cx_avr)

    return nothing
    
end



function get_a_plant_system_stability_matrices(
    plant; only_gen = false)


    if only_gen_in_plant(plant) || only_gen == true

        return  get_only_gen_plant_system_stability_matrices(
            plant)

    elseif only_gen_avr_in_plant(plant)

        return get_only_gen_avr_plant_system_stability_matrices(plant)
                
    elseif only_gen_gov_in_plant(plant)

        return get_only_gen_gov_plant_system_stability_matrices(plant)

    elseif gen_gov_avr_in_plant(plant)

        return get_gen_gov_avr_plant_system_stability_matrices(plant)
        
    else
        return nothing
    end
        
end



function get_a_plant_system_stability_matrices!(
    ( stab_Ax_view, stab_Bx_view, stab_Cx_view ),
    plant;
    only_gen = false,
    stab_τm_vf_view = nothing )


    if only_gen_in_plant(plant) || only_gen == true

        get_only_gen_plant_system_stability_matrices!(
            ( stab_Ax_view, stab_Bx_view, stab_Cx_view ),
            plant;
            stab_τm_vf_view = stab_τm_vf_view  )

    elseif only_gen_avr_in_plant(plant)

        get_only_gen_avr_plant_system_stability_matrices!(
            ( stab_Ax_view, stab_Bx_view, stab_Cx_view ),
            plant )
                
        
    elseif only_gen_gov_in_plant(plant)

        get_only_gen_gov_plant_system_stability_matrices!(
            ( stab_Ax_view, stab_Bx_view, stab_Cx_view ),
            plant )

    elseif gen_gov_avr_in_plant(plant)

        get_gen_gov_avr_plant_system_stability_matrices!(
            ( stab_Ax_view, stab_Bx_view, stab_Cx_view ),
            plant )
        
    else
        nothing
    end

    return nothing
        
end

#--------------------------------------------------------

function instantiate_a_plant_system_stability_matrices(
    plant;
    only_gen = true )


    if only_gen_in_plant( plant ) || only_gen == true

        return  get_only_gen_plant_system_stability_matrices( plant.Gen )

    elseif only_gen_avr_in_plant( plant )

        return get_only_gen_gov_plant_system_stability_matrices(
            plant.Gen,
            plant.Exc,
            plant.dummy_vr1,
            plant.dummy_vf_tilade  )
        
    elseif only_gen_gov_in_plant( plant )

        return get_only_gen_avr_plant_system_stability_matrices(
            plant.Gen,
            plant.Gov )

    elseif gen_gov_avr_in_plant( plant )

        return get_gen_gov_avr_plant_system_stability_matrices(
            plant.Gen,
            plant.Gov,
            plant.Exc,
            plant.dummy_vr1,
            plant.dummy_vf_tilade )
        
    else
        return nothing
    end
        
end



function instantiate_a_plant_system_stability_matrices!(
    ( stab_Ax_view, stab_Bx_view, stab_Cx_view ),
    plant;
    only_gen = false )


    if only_gen_in_plant( plant ) || only_gen == true

        get_only_gen_plant_system_stability_matrices!(
            ( stab_Ax_view, stab_Bx_view, stab_Cx_view ),
            plant.Gen )

    elseif only_gen_avr_in_plant( plant )

        get_only_gen_gov_plant_system_stability_matrices!(
            ( stab_Ax_view, stab_Bx_view, stab_Cx_view ),
            plant  )
        
    elseif only_gen_gov_in_plant( plant )

        get_only_gen_avr_plant_system_stability_matrices!(
            ( stab_Ax_view, stab_Bx_view, stab_Cx_view ),
            plant )

    elseif gen_gov_avr_in_plant( plant )

        return get_gen_gov_avr_plant_system_stability_matrices!( ( stab_Ax_view, stab_Bx_view, stab_Cx_view ), plant )
        
    else
        nothing
    end

    return nothing
        
end



#--------------------------------------------------------
#--------------------------------------------------------
#  create aggregate system matrices Idx and views
#--------------------------------------------------------
#--------------------------------------------------------

#--------------------------------------------------------
#--------------------------------------------------------
# post network-idx.jl, relocated from network-idx.jl
#--------------------------------------------------------
#--------------------------------------------------------


function create_flattened_type_states_dims_offset_Idx_for_gens_nodes( comp_collection::OrderedDict; only_gen = false )

    gens_nodes_collection = get_gens_nodes( comp_collection )
    
    if only_gen == false
                
        return create_flattened_type_dims_offset_Idx(
            get_components_state_vars,
            gens_nodes_collection)
    else
        sync_gens_collection = [
            a_plant.Gen
            for a_plant in gens_nodes_collection ]

        return create_flattened_type_dims_offset_Idx(
            get_components_state_vars, sync_gens_collection)
    end
    
        
end



function create_flattened_type_stab_states_dims_offset_Idx_for_gens_nodes(
    comp_collection::OrderedDict )

    gens_nodes_collection = get_gens_nodes(
        comp_collection )

    return create_flattened_type_dims_offset_Idx(
        get_components_stab_state_vars,
        gens_nodes_collection)    
        
end

         
function create_flattened_type_states_dims_offset_Idx_for_gens_nodes(
    gens_nodes_collection; only_gen = false  )

    if only_gen == false

        return create_flattened_type_dims_offset_Idx(
            get_components_state_vars,
            gens_nodes_collection )
        
    else
        
        sync_gens_collection = [
            a_plant.Gen
            for a_plant in gens_nodes_collection ]

        return create_flattened_type_dims_offset_Idx(
            get_components_state_vars,
            sync_gens_collection)

    end
            
end

function create_flattened_type_stab_states_dims_offset_Idx_for_gens_nodes(
    gens_nodes_collection )

    return create_flattened_type_dims_offset_Idx(
        get_components_stab_state_vars,
        gens_nodes_collection)    
        
end




function create_gens_nodes_aggregate_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = false )
    
    #--------------------------------------------------

    if only_gen == false 

        (nodes_dim_size,
         nodes_offset,
         nodes_state_Idx) = create_flattened_type_states_dims_offset_Idx_for_gens_nodes( gens_nodes_collection )
    
        Ax_size =  sum(nodes_dim_size)
    else

        (nodes_dim_size,
         nodes_offset,
         nodes_state_Idx) = create_flattened_type_states_dims_offset_Idx_for_gens_nodes(
            gens_nodes_collection; only_gen = true  )
    
        Ax_size =  sum(nodes_dim_size)
        

    end
    
    
    #--------------------------------------------------

    no_nodes = length( gens_nodes_collection )
    
    id_iq_ph_vh_vector  = [:id, :iq , :ph,  :vh]

    id_iq_ph_vh_dims    = [
        length( id_iq_ph_vh_vector )
        for k in 1:no_nodes ]

    id_iq_ph_vh_offsets = create_offsets(
        id_iq_ph_vh_dims )
    
    id_iq_ph_vh_idxs = create_idxs(
        id_iq_ph_vh_offsets,
        id_iq_ph_vh_dims )

    #--------------------------------------------------


    ω_ref_τm_v_ref_vector = [
        :ω_ref, :τm,  :v_ref, :p_order ]

    ω_ref_τm_v_ref_dims = [
        length( ω_ref_τm_v_ref_vector )
        for k in 1:no_nodes ]

    ω_ref_τm_v_ref_offsets = create_offsets(
        ω_ref_τm_v_ref_dims)
    
    ω_ref_τm_v_ref_idxs = create_idxs(
        ω_ref_τm_v_ref_offsets,
        ω_ref_τm_v_ref_dims )

    #--------------------------------------------------
    
    Bx_dims = [
        ( a_row, a_col )
        for ( a_row, a_col ) in zip(
            nodes_dim_size, id_iq_ph_vh_dims ) ]

    Cx_dims = [
        ( a_row, a_col )
        for ( a_row, a_col ) in zip(
            nodes_dim_size, ω_ref_τm_v_ref_dims ) ]
    
    #--------------------------------------------------
    
    Bx_row_dims    = [ first( Bxi ) for Bxi in Bx_dims ]
    
    Bx_row_offsets = create_offsets(Bx_row_dims)
    
    Bx_row_idxs    = create_idxs(Bx_row_offsets, Bx_row_dims)

    Bx_col_dims    = [ last( Bxi )
                       for Bxi in Bx_dims ]
    
    Bx_col_offsets = create_offsets(Bx_col_dims)
    
    Bx_col_idxs    = create_idxs(
        Bx_col_offsets, Bx_col_dims)
    
    Bx_row_size    = sum( Bx_row_dims )
    
    Bx_col_size    = sum( Bx_col_dims )
    
    
    #--------------------------------------------------
    
    Cx_row_dims    = [ first( Cxi )
                       for Cxi in Cx_dims ]
    
    Cx_row_offsets = create_offsets(Cx_row_dims)
    
    Cx_row_idxs    = create_idxs(
        Cx_row_offsets, Cx_row_dims )

    Cx_col_dims    = [ last( Cxi )
                       for Cxi in Cx_dims ]
    
    Cx_col_offsets = create_offsets(Cx_col_dims)
    
    Cx_col_idxs    = create_idxs(
        Cx_col_offsets, Cx_col_dims )
    
    Cx_row_size    = sum( Cx_row_dims )
    
    Cx_col_size    = sum( Cx_col_dims )

    #--------------------------------------------------
      
    Ax_matrix = spzeros( Ax_size, Ax_size )
    
    # Bx_matrix = sparse( 1.0I, Bx_row_size, Bx_row_size )
    
    Bx_matrix = spzeros( Bx_row_size, Bx_col_size )
    
    Cx_matrix = spzeros( Cx_row_size, Cx_col_size )

    #--------------------------------------------------

 
    vec_Ax_views = [
        view( Ax_matrix, node_state_Idx, node_state_Idx )
        for node_state_Idx in nodes_state_Idx ]

    vec_Bx_views = [
        view( Bx_matrix, Bx_row_idx, Bx_col_idx )
        for ( Bx_row_idx, Bx_col_idx ) in zip(
            Bx_row_idxs, Bx_col_idxs ) ]

    vec_Cx_views = [
        view( Cx_matrix, Cx_row_idx, Cx_col_idx )
        for ( Cx_row_idx, Cx_col_idx ) in zip(
            Cx_row_idxs, Cx_col_idxs ) ]
    

    #--------------------------------------------------

    if only_gen == true
        
        τm_vf_vector = [:τm, :vf ]

        τm_vf_dims = [ length( τm_vf_vector )
                       for k in 1:no_nodes ]

        Ax_τm_vf_dims = [
            (a_row, a_col) for (a_row, a_col) in zip(
                nodes_dim_size, τm_vf_dims ) ]

        Ax_τm_vf_row_dims = [
            first( Ax_τm_vfi )
            for Ax_τm_vfi in Ax_τm_vf_dims ]

        Ax_τm_vf_row_offsets = create_offsets(
            Ax_τm_vf_row_dims)

        Ax_τm_vf_row_idxs = create_idxs(
            Ax_τm_vf_row_offsets,
            Ax_τm_vf_row_dims)

        Ax_τm_vf_col_dims = [
            last( Ax_τm_vfi )
            for Ax_τm_vfi in Ax_τm_vf_dims ]

        Ax_τm_vf_col_offsets = create_offsets(
            Ax_τm_vf_col_dims)

        Ax_τm_vf_col_idxs = create_idxs(
            Ax_τm_vf_col_offsets,
            Ax_τm_vf_col_dims)

        Ax_τm_vf_row_size = sum( Ax_τm_vf_row_dims )

        Ax_τm_vf_col_size = sum( Ax_τm_vf_col_dims )
        

        Ax_τm_vf_matrix =  spzeros(
            Ax_τm_vf_row_size, Ax_τm_vf_col_size )
        
        # vec_τm_vf_views
        vec_Ax_τm_vf_views = [
            view(Ax_τm_vf_matrix,τm_vf_row_idx,τm_vf_col_idx)
            for ( τm_vf_row_idx, τm_vf_col_idx ) in zip(
                Ax_τm_vf_row_idxs, Ax_τm_vf_col_idxs ) ]

    end
    
    #--------------------------------------------------

    if only_gen == false

        Ax_Bx_Cx_views  = (; vec_Ax_views, vec_Bx_views,
                           vec_Cx_views )
        
        Ax_Bx_Cx_matrix = (; Ax_matrix, Bx_matrix, Cx_matrix )

        Bx_idxs         = (; Bx_row_idxs, Bx_col_idxs )

        Cx_idxs         = (; Cx_row_idxs, Cx_col_idxs )

        idxs            = (;nodes_state_Idx, Bx_idxs, Cx_idxs)
        
        
        return (; Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs )
        
    else

        Ax_Bx_Cx_views = (; vec_Ax_views, vec_Bx_views,
                          vec_Cx_views, vec_Ax_τm_vf_views )

        Ax_Bx_Cx_matrix = (; Ax_matrix, Bx_matrix,
                           Cx_matrix, Ax_τm_vf_matrix )

        Bx_idxs         = (; Bx_row_idxs, Bx_col_idxs )

        Cx_idxs         = (; Cx_row_idxs, Cx_col_idxs )

        τm_vf_idxs      = (; Ax_τm_vf_row_idxs,
                           Ax_τm_vf_col_idxs )

        idxs            = (; nodes_state_Idx,
                           Bx_idxs, Cx_idxs, τm_vf_idxs  )        
        
        return (; Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs )
    end
    
        
end

function create_gens_nodes_aggregate_system_stability_matrices_Idx_and_views(
    gens_nodes_collection; only_gen = false )
    
    #--------------------------------------------------

    (nodes_dim_size,
     nodes_offset,
     nodes_state_Idx) =
         create_flattened_type_stab_states_dims_offset_Idx_for_gens_nodes( gens_nodes_collection )
    
    Ax_size =  sum( nodes_dim_size )
    
    #--------------------------------------------------

    no_nodes            = length( gens_nodes_collection )
    
    id_iq_ph_vh_vector  = [ :id, :iq , :ph,  :vh ]

    id_iq_ph_vh_dims    = [
        length( id_iq_ph_vh_vector )
        for k in 1:no_nodes ]

    id_iq_ph_vh_offsets = create_offsets( id_iq_ph_vh_dims )
    
    id_iq_ph_vh_idxs    = create_idxs( id_iq_ph_vh_offsets,
                                    id_iq_ph_vh_dims )

    #--------------------------------------------------


    ω_ref_τm_v_ref_vector = [ :ω_ref, :τm,  :v_ref, :p_order ]

    ω_ref_τm_v_ref_dims = [ length( ω_ref_τm_v_ref_vector )
                            for k in 1:no_nodes ]

    ω_ref_τm_v_ref_offsets = create_offsets(
        ω_ref_τm_v_ref_dims )
    
    ω_ref_τm_v_ref_idxs = create_idxs(
        ω_ref_τm_v_ref_offsets,
        ω_ref_τm_v_ref_dims )

    #--------------------------------------------------
    
    Bx_dims = [ ( a_row, a_col )
                for ( a_row, a_col ) in
                    zip( nodes_dim_size,
                         id_iq_ph_vh_dims)]

    Cx_dims = [ ( a_row, a_col )
                for ( a_row, a_col ) in
                    zip(nodes_dim_size,
                        ω_ref_τm_v_ref_dims)]
    
    #--------------------------------------------------
    
    Bx_row_dims    = [ first( Bxi )
                       for Bxi in Bx_dims ]
    
    Bx_row_offsets = create_offsets( Bx_row_dims )
    
    Bx_row_idxs    = create_idxs(Bx_row_offsets, Bx_row_dims)

    Bx_col_dims    = [ last( Bxi ) for Bxi in Bx_dims ]
    
    Bx_col_offsets = create_offsets( Bx_col_dims )
    
    Bx_col_idxs    = create_idxs(Bx_col_offsets, Bx_col_dims )
    
    Bx_row_size    = sum( Bx_row_dims )
    
    Bx_col_size    = sum( Bx_col_dims )
    
    
    #--------------------------------------------------
    
    Cx_row_dims    = [ first( Cxi ) for Cxi in Cx_dims ]
    
    Cx_row_offsets = create_offsets( Cx_row_dims )
    
    Cx_row_idxs    = create_idxs(Cx_row_offsets, Cx_row_dims)

    Cx_col_dims    = [ last( Cxi ) for Cxi in Cx_dims ]
    
    Cx_col_offsets = create_offsets( Cx_col_dims )
    
    Cx_col_idxs    = create_idxs(Cx_col_offsets, Cx_col_dims)
    
    Cx_row_size    = sum( Cx_row_dims )
    
    Cx_col_size    = sum( Cx_col_dims )
    
    #--------------------------------------------------
      
    Ax_matrix = zeros( Ax_size, Ax_size )
    
    # Bx_matrix = sparse( 1.0I, Bx_row_size, Bx_row_size )

    Bx_matrix = zeros( Bx_row_size, Bx_col_size )
    
    Cx_matrix = zeros( Cx_row_size, Cx_col_size )

    #--#
    
    # Ax_matrix = sparse( 1.0I, Ax_size, Ax_size )
    
    # Bx_matrix = sparse( 1.0I, Bx_row_size, Bx_col_size )
    
    # Cx_matrix = sparse( 1.0I, Cx_row_size, Cx_col_size )    
    
    #--------------------------------------------------

 
    vec_Ax_views = [
        view( Ax_matrix, node_state_Idx, node_state_Idx )
        for node_state_Idx in nodes_state_Idx ]

    vec_Bx_views = [
        view( Bx_matrix, Bx_row_idx, Bx_col_idx )
        for ( Bx_row_idx, Bx_col_idx ) in
            zip( Bx_row_idxs, Bx_col_idxs ) ]

    vec_Cx_views = [
        view( Cx_matrix, Cx_row_idx, Cx_col_idx )
        for ( Cx_row_idx, Cx_col_idx ) in
            zip( Cx_row_idxs, Cx_col_idxs ) ]
    
     
    #--------------------------------------------------

    if only_gen == true
        
        τm_vf_vector  = [:τm, :vf ]

        τm_vf_dims    = [
            length( τm_vf_vector )
            for k in 1:no_nodes ]

        Ax_τm_vf_dims = [
            ( a_row, a_col )
            for ( a_row, a_col ) in
                zip(nodes_dim_size, τm_vf_dims ) ]


        Ax_τm_vf_row_dims = [
            first( Ax_τm_vfi )
            for Ax_τm_vfi in Ax_τm_vf_dims ]

        Ax_τm_vf_row_offsets = create_offsets(
            Ax_τm_vf_row_dims )

        Ax_τm_vf_row_idxs = create_idxs(
            Ax_τm_vf_row_offsets,
            Ax_τm_vf_row_dims )

        Ax_τm_vf_col_dims = [
            last( Ax_τm_vfi )
            for Ax_τm_vfi in Ax_τm_vf_dims ]

        Ax_τm_vf_col_offsets = create_offsets(
            Ax_τm_vf_col_dims )

        Ax_τm_vf_col_idxs = create_idxs(
            Ax_τm_vf_col_offsets, Ax_τm_vf_col_dims )

        Ax_τm_vf_row_size = sum( Ax_τm_vf_row_dims )

        Ax_τm_vf_col_size = sum( Ax_τm_vf_col_dims )        

        Ax_τm_vf_matrix =  zeros(
            Ax_τm_vf_row_size, Ax_τm_vf_col_size )
        
        vec_Ax_τm_vf_views = [
            view( Ax_τm_vf_matrix, τm_vf_row_idx,
                  τm_vf_col_idx )
            for ( τm_vf_row_idx, τm_vf_col_idx ) in
                zip( Ax_τm_vf_row_idxs, Ax_τm_vf_col_idxs ) ]

    end
    
    #--------------------------------------------------

    if only_gen == false

        Ax_Bx_Cx_views  = (; vec_Ax_views, vec_Bx_views,
                           vec_Cx_views )
        
        Ax_Bx_Cx_matrix = (; Ax_matrix, Bx_matrix, Cx_matrix )

        Bx_idxs         = (; Bx_row_idxs, Bx_col_idxs )

        Cx_idxs         = (; Cx_row_idxs, Cx_col_idxs )

        idxs            = (; nodes_state_Idx, Bx_idxs, Cx_idxs  )

        #--------------------------------------------------
        
        
        return (; Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs )
        
    else

        Ax_Bx_Cx_views  = (; vec_Ax_views, vec_Bx_views,
                           vec_Cx_views, vec_Ax_τm_vf_views )

        Ax_Bx_Cx_matrix = (; Ax_matrix, Bx_matrix,
                           Cx_matrix, Ax_τm_vf_matrix )

        Bx_idxs         = (; Bx_row_idxs, Bx_col_idxs )

        Cx_idxs         = (; Cx_row_idxs, Cx_col_idxs )

        τm_vf_idxs      = (; Ax_τm_vf_row_idxs,
                           Ax_τm_vf_col_idxs )

        idxs            = (; nodes_state_Idx, Bx_idxs,
                           Cx_idxs, τm_vf_idxs  )

        #--------------------------------------------------
        
        return (; Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs )
        
    end
    
        
end



function create_gens_nodes_system_matrices_Idx_and_views(
    gens_nodes_collection; only_gen = false )

    #--------------------------------------------------

    # nodes_dim_size, nodes_offset, nodes_state_Idx = create_flattened_type_stab_states_dims_offset_Idx_for_gens_nodes( gens_nodes_collection )

    nodes_dim_size, nodes_offset, nodes_state_Idx = create_flattened_type_states_dims_offset_Idx_for_gens_nodes(
        gens_nodes_collection )
    
    Ax_size =  sum(nodes_dim_size)
    
    #--------------------------------------------------

    no_nodes = length( gens_nodes_collection )
    
    id_iq_ph_vh_vector    = [:id, :iq , :ph,  :vh]

    id_iq_ph_vh_dims = [
        length( id_iq_ph_vh_vector )
        for k in 1:no_nodes ]

    id_iq_ph_vh_offsets = create_offsets(
        id_iq_ph_vh_dims )
    
    id_iq_ph_vh_idxs = create_idxs(
        id_iq_ph_vh_offsets,
        id_iq_ph_vh_dims )

    #--------------------------------------------------


    ω_ref_τm_v_ref_vector = [
        :ω_ref, :τm,  :v_ref, :p_order ]

    ω_ref_τm_v_ref_dims = [
        length( ω_ref_τm_v_ref_vector )
        for k in 1:no_nodes ]

    ω_ref_τm_v_ref_offsets = create_offsets(
        ω_ref_τm_v_ref_dims)
    
    ω_ref_τm_v_ref_idxs = create_idxs(
        ω_ref_τm_v_ref_offsets,
        ω_ref_τm_v_ref_dims )

    #--------------------------------------------------
    
    Bx_dims = [
        ( a_row, a_col )
        for ( a_row, a_col ) in zip(
            nodes_dim_size, id_iq_ph_vh_dims ) ]

    Cx_dims = [
        ( a_row, a_col )
        for ( a_row, a_col ) in zip(
            nodes_dim_size, ω_ref_τm_v_ref_dims ) ]
    
    #--------------------------------------------------
    
    Bx_row_dims    = [
        first( Bxi ) for Bxi in Bx_dims ]
    
    Bx_row_offsets = create_offsets(Bx_row_dims)
    
    Bx_row_idxs    = create_idxs(
        Bx_row_offsets, Bx_row_dims)

    Bx_col_dims    = [
        last( Bxi ) for Bxi in Bx_dims ]
    
    Bx_col_offsets = create_offsets(Bx_col_dims)
    
    Bx_col_idxs    = create_idxs(
        Bx_col_offsets, Bx_col_dims)
    
    Bx_row_size = sum( Bx_row_dims )
    
    Bx_col_size = sum( Bx_col_dims )
    
    
    #--------------------------------------------------
    
    Cx_row_dims = [
        first( Cxi )
        for Cxi in Cx_dims ]
    
    Cx_row_offsets = create_offsets(
        Cx_row_dims)
    
    Cx_row_idxs = create_idxs(
        Cx_row_offsets, Cx_row_dims )

    Cx_col_dims = [
        last( Cxi ) for Cxi in Cx_dims ]
    
    Cx_col_offsets = create_offsets(
        Cx_col_dims)
    
    Cx_col_idxs = create_idxs(
        Cx_col_offsets, Cx_col_dims )
    
    Cx_row_size = sum( Cx_row_dims )
    
    Cx_col_size = sum( Cx_col_dims )
    
    #--------------------------------------------------
      
    Ax_matrix = zeros( Ax_size, Ax_size )
    
    # Bx_matrix = zeros( Bx_row_size, Bx_row_size )
    
    Bx_matrix = zeros( Bx_row_size, Bx_col_size )
    
    Cx_matrix = zeros( Cx_row_size, Cx_col_size )

    #--------------------------------------------------
 
    vec_Ax_views = [
        view( Ax_matrix, node_state_Idx, node_state_Idx )
        for node_state_Idx in nodes_state_Idx ]

    vec_Bx_views = [
        view( Bx_matrix, Bx_row_idx, Bx_col_idx )
        for ( Bx_row_idx, Bx_col_idx ) in zip(
            Bx_row_idxs, Bx_col_idxs ) ]

    vec_Cx_views = [
        view( Cx_matrix, Cx_row_idx, Cx_col_idx )
        for ( Cx_row_idx, Cx_col_idx ) in zip(
            Cx_row_idxs, Cx_col_idxs ) ]
    
    #--------------------------------------------------
    # Stability indices and matrices
    # The element (δ) of the first index is filtered out
    #--------------------------------------------------

    filter_stab_nodes_state_Idx = [
        first(Idx)+1:last(Idx)
        for Idx in nodes_state_Idx ]

    filter_stab_Bx_row_idxs     = [
        first(Idx)+1:last(Idx)
        for Idx in Bx_row_idxs ]
    
    filter_stab_Cx_row_idxs     = [
        first(Idx)+1:last(Idx)
        for Idx in Cx_row_idxs ]

    
    stab_Ax_dims     = length.(
        filter_stab_nodes_state_Idx )
    
    stab_Bx_row_dims = length.(
        filter_stab_Bx_row_idxs )
    
    stab_Cx_row_dims = length.(
        filter_stab_Cx_row_idxs )

    stab_Ax_offsets = create_offsets(
        stab_Ax_dims)
    
    stab_Ax_idxs    = create_idxs(
        stab_Ax_offsets, stab_Ax_dims )      

    stab_Bx_row_offsets = create_offsets(
        stab_Bx_row_dims)
    
    stab_Bx_row_idxs    = create_idxs(
        stab_Bx_row_offsets, stab_Bx_row_dims )     

    stab_Cx_row_offsets = create_offsets(
        stab_Cx_row_dims )
    
    stab_Cx_row_idxs    = create_idxs(
        stab_Cx_row_offsets, stab_Cx_row_dims )    

   #--------------------------------------------------

    stab_Ax_size     = sum( stab_Ax_dims )
    stab_Bx_row_size = sum( stab_Bx_row_dims )
    stab_Cx_row_size = sum( stab_Cx_row_dims )
    
    stab_Ax_matrix = zeros(
        stab_Ax_size, stab_Ax_size )
    
    stab_Bx_matrix = zeros(
        stab_Bx_row_size, Bx_col_size )
    
    stab_Cx_matrix = zeros(
        stab_Cx_row_size, Cx_col_size )

    #--------------------------------------------------
    # create views that can be used to link `stab_Ax_matrix`
    # to `Ax_matrix`

    vec_stab_link_Ax_views = [
        view( stab_Ax_matrix, Ax_idx, Ax_idx )
        for Ax_idx in stab_Ax_idxs ]

    vec_stab_link_Bx_views = [
        view( stab_Bx_matrix, stab_Bx_row_idx, Bx_col_idx )
        for ( stab_Bx_row_idx, Bx_col_idx ) in zip(
            stab_Bx_row_idxs, Bx_col_idxs ) ]

    vec_stab_link_Cx_views = [
        view( stab_Cx_matrix, stab_Cx_row_idx, Cx_col_idx )
        for ( stab_Cx_row_idx, Cx_col_idx ) in zip(
            stab_Cx_row_idxs, Cx_col_idxs ) ]
     
    #--------------------------------------------------

    vec_stab_Ax_views = [
        view( Ax_matrix, node_stab_state_Idx,
              node_stab_state_Idx )
        for node_stab_state_Idx in
            filter_stab_nodes_state_Idx ]

    vec_stab_Bx_views = [
        view( Bx_matrix, Bx_stab_row_idx, Bx_col_idx )
        for ( Bx_stab_row_idx, Bx_col_idx ) in
            zip( filter_stab_Bx_row_idxs, Bx_col_idxs ) ]

    vec_stab_Cx_views = [
        view( Cx_matrix, Cx_stab_row_idx, Cx_col_idx )
        for ( Cx_stab_row_idx, Cx_col_idx ) in zip(
                filter_stab_Cx_row_idxs, Cx_col_idxs ) ]
     
    #--------------------------------------------------

    if only_gen == true
        
        τm_vf_vector = [:τm, :vf ]

        τm_vf_dims = [
            length( τm_vf_vector )
            for k in 1:no_nodes ]

        Ax_τm_vf_dims        = [
            ( a_row, a_col )
            for ( a_row, a_col ) in zip(
                nodes_dim_size, τm_vf_dims ) ]

        Ax_τm_vf_row_dims    = [
            first( Ax_τm_vfi )
            for Ax_τm_vfi in Ax_τm_vf_dims ]

        Ax_τm_vf_row_offsets = create_offsets(
            Ax_τm_vf_row_dims )

        Ax_τm_vf_row_idxs    = create_idxs(
            Ax_τm_vf_row_offsets, Ax_τm_vf_row_dims )

        Ax_τm_vf_col_dims    = [
            last( Ax_τm_vfi )
            for Ax_τm_vfi in Ax_τm_vf_dims ]

        Ax_τm_vf_col_offsets = create_offsets(
            Ax_τm_vf_col_dims )

        Ax_τm_vf_col_idxs    = create_idxs(
            Ax_τm_vf_col_offsets, Ax_τm_vf_col_dims )

        Ax_τm_vf_row_size    = sum( Ax_τm_vf_row_dims )

        Ax_τm_vf_col_size    = sum( Ax_τm_vf_col_dims )
        

        Ax_τm_vf_matrix      = zeros(
            Ax_τm_vf_row_size, Ax_τm_vf_col_size )
        
        vec_Ax_τm_vf_views      = [
            view( Ax_τm_vf_matrix, τm_vf_row_idx,
                  τm_vf_col_idx )
            for ( τm_vf_row_idx, τm_vf_col_idx ) in zip(
                Ax_τm_vf_row_idxs, Ax_τm_vf_col_idxs ) ]

        #--------------------------------------------------
        # Stability indices and matrices
        #--------------------------------------------------

        filter_stab_Ax_τm_vf_row_idxs = [
            first(Idx)+1:last(Idx)
            for Idx in Ax_τm_vf_row_idxs ]
        
        vec_stab_τm_vf_views   = [
            view( Ax_τm_vf_matrix, stab_τm_vf_row_idx,
                  τm_vf_col_idx )
            for ( stab_τm_vf_row_idx, τm_vf_col_idx ) in zip(
                filter_stab_Ax_τm_vf_row_idxs,
                Ax_τm_vf_col_idxs ) ]

        stab_Ax_τm_vf_row_dims  = length.(
            filter_stab_Ax_τm_vf_row_idxs )
        
        stab_Ax_τm_vf_row_offsets = create_offsets(
            stab_Ax_τm_vf_row_dims)
        
        stab_Ax_τm_vf_row_idxs  = create_idxs(
            stab_Ax_τm_vf_row_offsets,
            stab_Ax_τm_vf_row_dims )

        stab_Ax_τm_vf_row_size  = sum(
            stab_Ax_τm_vf_row_dims )

        stab_Ax_τm_vf_matrix  = zeros(
            stab_Ax_τm_vf_row_size,
            Ax_τm_vf_col_size ) 

        vec_stab_link_τm_vf_views   = [
            view( stab_Ax_τm_vf_matrix,
                  stab_τm_vf_row_idx,
                  τm_vf_col_idx )
            for ( stab_τm_vf_row_idx, τm_vf_col_idx ) in zip(
                stab_Ax_τm_vf_row_idxs, Ax_τm_vf_col_idxs ) ]
        
    end
    
    #--------------------------------------------------

    if only_gen == false

        Ax_Bx_Cx_views = (
            ; vec_Ax_views, vec_Bx_views,
            vec_Cx_views )
        
        Ax_Bx_Cx_matrix = (
            ; Ax_matrix, Bx_matrix, Cx_matrix )

        Bx_idxs = (
            ; Bx_row_idxs, Bx_col_idxs )

        Cx_idxs = (
            ; Cx_row_idxs, Cx_col_idxs )

        idxs = (
            ; nodes_state_Idx, Bx_idxs, Cx_idxs  )

        #--------------------------------------------------
        
        # Stability
        
        stab_Ax_Bx_Cx_views      = (
            ; vec_stab_Ax_views,
            vec_stab_Bx_views,
            vec_stab_Cx_views )
        
        stab_link_Ax_Bx_Cx_views = (
            ; vec_stab_link_Ax_views,
            vec_stab_link_Bx_views,
            vec_stab_link_Cx_views )
        
        
        stab_Ax_Bx_Cx_matrix = (
            ; stab_Ax_matrix,
            stab_Bx_matrix,
            stab_Cx_matrix )


        stab_Bx_idxs     = (
            ; stab_Bx_row_idxs,
            Bx_col_idxs )

        stab_Cx_idxs     = (
            ; stab_Cx_row_idxs,
            Cx_col_idxs )

        stab_idxs        = (
            ; stab_Ax_idxs,
            stab_Bx_idxs,
            stab_Cx_idxs  )

        filter_stab_idxs = (
            ; filter_stab_nodes_state_Idx,
            filter_stab_Bx_row_idxs,
            filter_stab_Cx_row_idxs  )
        
        #--------------------------------------------------
        
        return (; Ax_Bx_Cx_views, stab_Ax_Bx_Cx_views, stab_link_Ax_Bx_Cx_views ), (; Ax_Bx_Cx_matrix,  stab_Ax_Bx_Cx_matrix ), (; idxs, stab_idxs,  filter_stab_idxs )
        
    else

        Ax_Bx_Cx_views = (
            ; vec_Ax_views, vec_Bx_views, vec_Cx_views,
            vec_Ax_τm_vf_views )

        Ax_Bx_Cx_matrix = (
            ; Ax_matrix, Bx_matrix, Cx_matrix,
            Ax_τm_vf_matrix )

        Bx_idxs    = (; Bx_row_idxs, Bx_col_idxs )

        Cx_idxs    = (; Cx_row_idxs, Cx_col_idxs )

        τm_vf_idxs = (; Ax_τm_vf_row_idxs, Ax_τm_vf_col_idxs )

        idxs       = (; nodes_state_Idx, Bx_idxs,
                      Cx_idxs, τm_vf_idxs  )

        #--------------------------------------------------

        # Stability
        
        stab_Ax_Bx_Cx_views = (; vec_stab_Ax_views,
                               vec_stab_Bx_views,
                               vec_stab_Cx_views,
                               vec_stab_τm_vf_views )

        stab_link_Ax_Bx_Cx_views = (
            ; vec_stab_link_Ax_views,
            vec_stab_link_Bx_views,
            vec_stab_link_Cx_views,
            vec_stab_link_τm_vf_views )
        
        stab_Ax_Bx_Cx_matrix = (
            ; stab_Ax_matrix, stab_Bx_matrix,
            stab_Cx_matrix, stab_Ax_τm_vf_matrix )

        stab_Bx_idxs    = (
            ; stab_Bx_row_idxs, Bx_col_idxs )

        stab_Cx_idxs    = (
            ; stab_Cx_row_idxs, Cx_col_idxs )

        stab_τm_vf_idxs = (
            ; stab_Ax_τm_vf_row_idxs,
            Ax_τm_vf_col_idxs )

        stab_idxs        = (
            ; stab_Ax_idxs,
            stab_Bx_idxs, stab_Cx_idxs,
            stab_τm_vf_idxs  )

        filter_stab_idxs = (
            ; filter_stab_nodes_state_Idx,
            filter_stab_Bx_row_idxs,
            filter_stab_Cx_row_idxs,
            filter_stab_Ax_τm_vf_row_idxs  )
                
        #--------------------------------------------------

        return (; Ax_Bx_Cx_views, stab_Ax_Bx_Cx_views,
                stab_link_Ax_Bx_Cx_views ),
        (; Ax_Bx_Cx_matrix,  stab_Ax_Bx_Cx_matrix ),
        (; idxs, stab_idxs,  filter_stab_idxs )
        
        # return ( Ax_Bx_Cx_views, Ax_Bx_Cx_matrix, idxs ), ( stab_Ax_Bx_Cx_views, stab_Ax_Bx_Cx_matrix, stab_idxs )
        
    end
    
    
        
end



function create_gens_nodes_system_blockdiag_matrices_Idx_and_views(
    gens_nodes_collection; only_gen = false )

    #--------------------------------------------------

    # nodes_dim_size, nodes_offset, nodes_state_Idx = create_flattened_type_stab_states_dims_offset_Idx_for_gens_nodes( gens_nodes_collection )

    (nodes_dim_size,
     nodes_offset,
     nodes_state_Idx) = create_flattened_type_states_dims_offset_Idx_for_gens_nodes( gens_nodes_collection )
    
    Ax_size =  sum(nodes_dim_size)
    
    #--------------------------------------------------

    no_nodes = length( gens_nodes_collection )
    
    id_iq_ph_vh_vector    = [:id, :iq , :ph,  :vh]

    id_iq_ph_vh_dims = [ length( id_iq_ph_vh_vector )
                         for k in 1:no_nodes ]

    id_iq_ph_vh_offsets = create_offsets( id_iq_ph_vh_dims )
    
    id_iq_ph_vh_idxs = create_idxs( id_iq_ph_vh_offsets,
                                    id_iq_ph_vh_dims )

    #--------------------------------------------------


    ω_ref_τm_v_ref_vector = [:ω_ref, :τm,  :v_ref, :p_order ]

    ω_ref_τm_v_ref_dims = [ length( ω_ref_τm_v_ref_vector )
                            for k in 1:no_nodes ]

    ω_ref_τm_v_ref_offsets = create_offsets(
        ω_ref_τm_v_ref_dims)
    
    ω_ref_τm_v_ref_idxs = create_idxs(
        ω_ref_τm_v_ref_offsets,
        ω_ref_τm_v_ref_dims )

    #--------------------------------------------------
    
    Bx_dims = [ ( a_row, a_col )
                for ( a_row, a_col ) in
                    zip(nodes_dim_size,
                        id_iq_ph_vh_dims ) ]

    Cx_dims = [ ( a_row, a_col )
                for ( a_row, a_col ) in
                    zip(nodes_dim_size,
                        ω_ref_τm_v_ref_dims ) ]
    
    #--------------------------------------------------
    
    Bx_row_dims    = [ first( Bxi ) for Bxi in Bx_dims ]
    
    Bx_row_offsets = create_offsets(Bx_row_dims)
    
    Bx_row_idxs    = create_idxs(Bx_row_offsets, Bx_row_dims)

    Bx_col_dims    = [ last( Bxi ) for Bxi in Bx_dims ]
    
    Bx_col_offsets = create_offsets(Bx_col_dims)
    
    Bx_col_idxs    = create_idxs(Bx_col_offsets, Bx_col_dims)
    
    Bx_row_size = sum( Bx_row_dims )
    
    Bx_col_size = sum( Bx_col_dims )
    
    
    #--------------------------------------------------
    
    Cx_row_dims = [ first( Cxi ) for Cxi in Cx_dims ]
    
    Cx_row_offsets = create_offsets(Cx_row_dims)
    
    Cx_row_idxs = create_idxs( Cx_row_offsets, Cx_row_dims )

    Cx_col_dims = [ last( Cxi ) for Cxi in Cx_dims ]
    
    Cx_col_offsets = create_offsets(Cx_col_dims)
    
    Cx_col_idxs = create_idxs( Cx_col_offsets, Cx_col_dims )
    
    Cx_row_size = sum( Cx_row_dims )
    
    Cx_col_size = sum( Cx_col_dims )
    
    #--------------------------------------------------
      
    Ax_matrix = zeros( Ax_size, Ax_size )
    
    # Bx_matrix = zeros( Bx_row_size, Bx_row_size )
    
    Bx_matrix = zeros( Bx_row_size, Bx_col_size )
    
    Cx_matrix = zeros( Cx_row_size, Cx_col_size )

    #--------------------------------------------------
 
    vec_Ax_views = [
        view( Ax_matrix,
             node_state_Idx,
             node_state_Idx )
        for node_state_Idx in
                        nodes_state_Idx ]

    vec_Bx_views = [
        view(Bx_matrix,
             Bx_row_idx,
             Bx_col_idx )
        for ( Bx_row_idx, Bx_col_idx ) in
                        zip( Bx_row_idxs,
                             Bx_col_idxs ) ]

    vec_Cx_views = [
        view( Cx_matrix,
             Cx_row_idx,
             Cx_col_idx )
        for ( Cx_row_idx, Cx_col_idx ) in
                        zip( Cx_row_idxs,
                             Cx_col_idxs ) ]
    
    #--------------------------------------------------
    # Stability indices and matrices
    # The element (δ) of the first index is filtered out
    #--------------------------------------------------

    filter_stab_nodes_state_Idx = [
        first(Idx)+1:last(Idx)
        for Idx in nodes_state_Idx ]

    filter_stab_Bx_row_idxs = [
        first(Idx)+1:last(Idx)
        for Idx in Bx_row_idxs ]
    
    filter_stab_Cx_row_idxs = [
        first(Idx)+1:last(Idx)
        for Idx in Cx_row_idxs ]

    
    stab_Ax_dims     = length.( filter_stab_nodes_state_Idx ) 
    stab_Bx_row_dims = length.( filter_stab_Bx_row_idxs ) 
    stab_Cx_row_dims = length.( filter_stab_Cx_row_idxs )

    stab_Ax_offsets  = create_offsets(stab_Ax_dims)    
    stab_Ax_idxs = create_idxs(
        stab_Ax_offsets, stab_Ax_dims )      

    stab_Bx_row_offsets = create_offsets(stab_Bx_row_dims)    
    stab_Bx_row_idxs = create_idxs(
        stab_Bx_row_offsets, stab_Bx_row_dims )     

    stab_Cx_row_offsets = create_offsets(
        stab_Cx_row_dims )
    
    stab_Cx_row_idxs = create_idxs(
        stab_Cx_row_offsets, stab_Cx_row_dims )    

   #--------------------------------------------------

    vec_stab_Ax_views = [
        view(Ax_matrix,
            node_stab_state_Idx,
            node_stab_state_Idx )
        for node_stab_state_Idx in
            filter_stab_nodes_state_Idx ]

    vec_stab_Bx_views = [
        view( Bx_matrix, Bx_stab_row_idx, Bx_col_idx )
        for ( Bx_stab_row_idx, Bx_col_idx ) in
            zip( filter_stab_Bx_row_idxs, Bx_col_idxs ) ]

    vec_stab_Cx_views = [
        view( Cx_matrix, Cx_stab_row_idx, Cx_col_idx )
        for ( Cx_stab_row_idx, Cx_col_idx ) in
            zip( filter_stab_Cx_row_idxs, Cx_col_idxs ) ]
    
    #--------------------------------------------------

    stab_Ax_size     = sum( stab_Ax_dims )
    stab_Bx_row_size = sum( stab_Bx_row_dims )
    stab_Cx_row_size = sum( stab_Cx_row_dims )
    
    # stab_Ax_matrix = zeros(stab_Ax_size, stab_Ax_size )    
    # stab_Bx_matrix = zeros(stab_Bx_row_size,Bx_col_size)    
    # stab_Cx_matrix = zeros( stab_Cx_row_size, Cx_col_size )
    
    #--------------------------------------------------
    
    stab_Ax_matrix = blockdiag( vec_stab_Ax_views... )
    
    stab_Bx_matrix = blockdiag( vec_stab_Bx_views... )
    
    stab_Cx_matrix = blockdiag( vec_stab_Cx_views )
     
    #--------------------------------------------------

    if only_gen == true
        
        τm_vf_vector = [:τm, :vf ]

        τm_vf_dims = [ length( τm_vf_vector ) for k in 1:no_nodes ]

        Ax_τm_vf_dims        = [ ( a_row, a_col ) for ( a_row, a_col ) in zip(nodes_dim_size, τm_vf_dims ) ]

        Ax_τm_vf_row_dims    = [ first( Ax_τm_vfi ) for Ax_τm_vfi in Ax_τm_vf_dims ]

        Ax_τm_vf_row_offsets = create_offsets( Ax_τm_vf_row_dims )

        Ax_τm_vf_row_idxs    = create_idxs( Ax_τm_vf_row_offsets, Ax_τm_vf_row_dims )

        Ax_τm_vf_col_dims    = [ last( Ax_τm_vfi ) for Ax_τm_vfi in Ax_τm_vf_dims ]

        Ax_τm_vf_col_offsets = create_offsets( Ax_τm_vf_col_dims )

        Ax_τm_vf_col_idxs    = create_idxs( Ax_τm_vf_col_offsets, Ax_τm_vf_col_dims )

        Ax_τm_vf_row_size    = sum( Ax_τm_vf_row_dims )

        Ax_τm_vf_col_size    = sum( Ax_τm_vf_col_dims )
        

        Ax_τm_vf_matrix      = zeros( Ax_τm_vf_row_size, Ax_τm_vf_col_size )
        
        vec_Ax_τm_vf_views      = [ view( Ax_τm_vf_matrix, τm_vf_row_idx, τm_vf_col_idx ) for ( τm_vf_row_idx, τm_vf_col_idx ) in zip( Ax_τm_vf_row_idxs, Ax_τm_vf_col_idxs ) ]

        #--------------------------------------------------
        # Stability indices and matrices
        #--------------------------------------------------

        filter_stab_Ax_τm_vf_row_idxs = [ first(Idx)+1:last(Idx) for Idx in Ax_τm_vf_row_idxs ]
        vec_stab_τm_vf_views   = [ view( Ax_τm_vf_matrix, stab_τm_vf_row_idx, τm_vf_col_idx ) for ( stab_τm_vf_row_idx, τm_vf_col_idx ) in zip( filter_stab_Ax_τm_vf_row_idxs, Ax_τm_vf_col_idxs ) ]

        stab_Ax_τm_vf_row_dims  = length.( filter_stab_Ax_τm_vf_row_idxs )
        
        stab_Ax_τm_vf_row_offsets = create_offsets(stab_Ax_τm_vf_row_dims)
        
        stab_Ax_τm_vf_row_idxs  = create_idxs( stab_Ax_τm_vf_row_offsets, stab_Ax_τm_vf_row_dims )

        stab_Ax_τm_vf_row_size  = sum( stab_Ax_τm_vf_row_dims )

        stab_Ax_τm_vf_matrix  = blockdiag( vec_stab_τm_vf_views  ) 
        
    end
    
    #--------------------------------------------------

    if only_gen == false

        Ax_Bx_Cx_views = (; vec_Ax_views,
                          vec_Bx_views, vec_Cx_views )
        
        Ax_Bx_Cx_matrix = (; Ax_matrix,
                           Bx_matrix, Cx_matrix )

        Bx_idxs = (; Bx_row_idxs, Bx_col_idxs )

        Cx_idxs = (; Cx_row_idxs, Cx_col_idxs )

        idxs = (; nodes_state_Idx,
                Bx_idxs, Cx_idxs  )

        #--------------------------------------------------
        
        # Stability
        
        stab_Ax_Bx_Cx_views  =
            (; vec_stab_Ax_views,
             vec_stab_Bx_views, vec_stab_Cx_views )
        
        stab_Ax_Bx_Cx_matrix =
            (; stab_Ax_matrix,
             stab_Bx_matrix, stab_Cx_matrix )


        stab_Bx_idxs =
            (; stab_Bx_row_idxs, Bx_col_idxs )

        stab_Cx_idxs  =
            (; stab_Cx_row_idxs, Cx_col_idxs )

        stab_idxs  =
            (; stab_Ax_idxs, stab_Bx_idxs,
             stab_Cx_idxs  )

        filter_stab_idxs =
            (; filter_stab_nodes_state_Idx,
             filter_stab_Bx_row_idxs,
             filter_stab_Cx_row_idxs  )

        
        #--------------------------------------------------
        
        return (; Ax_Bx_Cx_views, stab_Ax_Bx_Cx_views  ),
        (; Ax_Bx_Cx_matrix,  stab_Ax_Bx_Cx_matrix ),
        (; idxs, stab_idxs,  filter_stab_idxs )
        
    else

        Ax_Bx_Cx_views =
            (; vec_Ax_views, vec_Bx_views,
             vec_Cx_views, vec_Ax_τm_vf_views )

        Ax_Bx_Cx_matrix =
            (; Ax_matrix, Bx_matrix,
             Cx_matrix, Ax_τm_vf_matrix )

        Bx_idxs = (; Bx_row_idxs, Bx_col_idxs )

        Cx_idxs = (; Cx_row_idxs, Cx_col_idxs )

        τm_vf_idxs = (; Ax_τm_vf_row_idxs,
                      Ax_τm_vf_col_idxs )

        idxs = (; nodes_state_Idx, Bx_idxs,
                Cx_idxs, τm_vf_idxs  )

        #--------------------------------------------------

        # Stability
        
        stab_Ax_Bx_Cx_views =
            (; vec_stab_Ax_views,
             vec_stab_Bx_views, vec_stab_Cx_views,
             vec_stab_τm_vf_views )
        
        stab_Ax_Bx_Cx_matrix =
            (; stab_Ax_matrix, stab_Bx_matrix,
             stab_Cx_matrix, stab_Ax_τm_vf_matrix )

        stab_Bx_idxs =
            (; stab_Bx_row_idxs, Bx_col_idxs )

        stab_Cx_idxs =
            (; stab_Cx_row_idxs, Cx_col_idxs )

        stab_τm_vf_idxs =
            (; stab_Ax_τm_vf_row_idxs,
             Ax_τm_vf_col_idxs )

        stab_idxs  =
            (; stab_Ax_idxs, stab_Bx_idxs,
             stab_Cx_idxs, stab_τm_vf_idxs  )

        filter_stab_idxs =
            (; filter_stab_nodes_state_Idx,
             filter_stab_Bx_row_idxs,
             filter_stab_Cx_row_idxs,
             filter_stab_Ax_τm_vf_row_idxs  )
                
        #--------------------------------------------------

         return (; Ax_Bx_Cx_views, stab_Ax_Bx_Cx_views ), (; Ax_Bx_Cx_matrix,  stab_Ax_Bx_Cx_matrix ), (; idxs, stab_idxs,  filter_stab_idxs )
        
        
    end
    
        
end

#--------------------------------------------------------
#--------------------------------------------------------



function make_system_matrices_with_link_matrices( gens_nodes_collection ;only_gen = false )
    
    """
    update is done by `update_gens_nodes_linked_system_matrices!( system_views, system_matrices, system_idxs, state, gens_nodes_collection ; only_gen = only_gen )`
    """

    system_views, system_matrices, system_idxs =  create_gens_nodes_system_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )

        
    Ax_Bx_Cx_views, stab_Ax_Bx_Cx_views, stab_link_Ax_Bx_Cx_views = system_views
    
    Ax_Bx_Cx_matrix,  stab_Ax_Bx_Cx_matrix = system_matrices
    
    idxs, stab_idxs,  filter_stab_idxs = system_idxs

    if only_gen == false
    
        vec_Ax_views, vec_Bx_views, vec_Cx_views = Ax_Bx_Cx_views
        
        vec_stab_Ax_views, vec_stab_Bx_views, vec_stab_Cx_views = stab_Ax_Bx_Cx_views
        
        vec_stab_link_Ax_views, vec_stab_link_Bx_views, vec_stab_link_Cx_views = stab_link_Ax_Bx_Cx_views

        Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix
        
        stab_Ax_matrix, stab_Bx_matrix, stab_Cx_matrix = stab_Ax_Bx_Cx_matrix

        nodes_state_Idx, Bx_idxs, Cx_idxs = idxs
        
        stab_Ax_idxs, stab_Bx_idxs, stab_Cx_idxs = stab_idxs
        
        filter_stab_nodes_state_Idx, filter_stab_Bx_row_idxs, filter_stab_Cx_row_idxs = filter_stab_idxs

        init_plants_system_matrices_views!(( vec_Ax_views, vec_Bx_views, vec_Cx_views ), gens_nodes_collection; only_gen = only_gen, vec_Ax_τm_vf_views = nothing )
        
    else
        
        vec_Ax_views, vec_Bx_views, vec_Cx_views, vec_Ax_τm_vf_views = Ax_Bx_Cx_views
        
        vec_stab_Ax_views, vec_stab_Bx_views, vec_stab_Cx_views, vec_stab_τm_vf_views  = stab_Ax_Bx_Cx_views
        
        vec_stab_link_Ax_views, vec_stab_link_Bx_views, vec_stab_link_Cx_views, vec_stab_link_τm_vf_views = stab_link_Ax_Bx_Cx_views

        Ax_matrix, Bx_matrix, Cx_matrix, Ax_τm_vf_matrix = Ax_Bx_Cx_matrix
        
        stab_Ax_matrix, stab_Bx_matrix, stab_Cx_matrix, stab_Ax_τm_vf_matrix = stab_Ax_Bx_Cx_matrix

        nodes_state_Idx, Bx_idxs, Cx_idxs, τm_vf_idxs  = idxs
        
        stab_Bx_row_idxs, Bx_col_idxs = stab_idxs
        
        filter_stab_nodes_state_Idx, filter_stab_Bx_row_idxs, filter_stab_Cx_row_idxs, filter_stab_Ax_τm_vf_row_idxs = filter_stab_idxs

        init_plants_system_matrices_views!(( vec_Ax_views, vec_Bx_views, vec_Cx_views ), gens_nodes_collection; only_gen = only_gen, vec_Ax_τm_vf_views = vec_Ax_τm_vf_views )
        
    end

    return (; system_views, system_matrices, system_idxs )
    
end


function make_industrial_model_system_matrices(
    system_views,
    system_matrices,
    system_idxs,
    state,
    gens_nodes_collection;
    only_gen = false )

    update_gens_nodes_linked_system_matrices!(
        system_views,
        system_matrices,
        system_idxs,
        state,
        gens_nodes_collection;
        only_gen = only_gen )
        
    (Ax_Bx_Cx_views,
     stab_Ax_Bx_Cx_views,
     stab_link_Ax_Bx_Cx_views) =
         system_views
    
    (Ax_Bx_Cx_matrix,
     stab_Ax_Bx_Cx_matrix) =
         system_matrices
    
    (idxs,
     stab_idxs,
     filter_stab_idxs) =
         system_idxs

    if only_gen == false
    
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views) =
             Ax_Bx_Cx_views
        
        (vec_stab_Ax_views,
         vec_stab_Bx_views,
         vec_stab_Cx_views) =
             stab_Ax_Bx_Cx_views
        
        (vec_stab_link_Ax_views,
         vec_stab_link_Bx_views,
         vec_stab_link_Cx_views) =
             stab_link_Ax_Bx_Cx_views

        (Ax_matrix,
         Bx_matrix,
         Cx_matrix) =
             Ax_Bx_Cx_matrix
        
        (stab_Ax_matrix,
         stab_Bx_matrix,
         stab_Cx_matrix) =
             stab_Ax_Bx_Cx_matrix

        (nodes_state_Idx,
         Bx_idxs,
         Cx_idxs) = idxs
        
        (stab_Ax_idxs,
         stab_Bx_idxs,
         stab_Cx_idxs) =
             stab_idxs
        
        (filter_stab_nodes_state_Idx,
         filter_stab_Bx_row_idxs,
         filter_stab_Cx_row_idxs) =
             filter_stab_idxs
        
    else
        
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         vec_Ax_τm_vf_views) =
             Ax_Bx_Cx_views
        
        (vec_stab_Ax_views,
         vec_stab_Bx_views,
         vec_stab_Cx_views,
         vec_stab_τm_vf_views)  =
             stab_Ax_Bx_Cx_views
        
        (vec_stab_link_Ax_views,
         vec_stab_link_Bx_views,
         vec_stab_link_Cx_views,
         vec_stab_link_τm_vf_views) =
            stab_link_Ax_Bx_Cx_views

        (Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         Ax_τm_vf_matrix) =
            Ax_Bx_Cx_matrix
        
        (stab_Ax_matrix,
         stab_Bx_matrix,
         stab_Cx_matrix,
         stab_Ax_τm_vf_matrix) =
             stab_Ax_Bx_Cx_matrix

        (nodes_state_Idx,
         Bx_idxs,
         Cx_idxs,
         τm_vf_idxs)  =
             idxs
        
        (stab_Bx_row_idxs,
         Bx_col_idxs) =
             stab_idxs
        
        (filter_stab_nodes_state_Idx,
         filter_stab_Bx_row_idxs,
         filter_stab_Cx_row_idxs,
         filter_stab_Ax_τm_vf_row_idxs) =
             filter_stab_idxs
        
    end

    if only_gen == false

        return (; vec_Ax_views,
                Ax_matrix,
                Bx_matrix,
                Cx_matrix )
    else
        return (; vec_Ax_views,
                Ax_matrix,
                Bx_matrix,
                Cx_matrix,
                Ax_τm_vf_matrix,
                vec_Ax_τm_vf_views )
    end
    
    
end



function make_system_matrices_inputs(
    netd, bus_dict_init )


    nodes_init_data_from_pf =
        get_nodes_init_data_from_pf(
            netd, bus_dict_init )

   #---------------------------------------------------

    gens_δ_id_iq_vd_vq_ed_dash_eq_dash_etc_from_pf = get_gens_δ_id_iq_vd_vq_ed_dash_eq_dash_etc_from_pf( netd, bus_dict_init )

    #---------------------------------------------------

    gens_δ_ω_ed_dash_eq_dash_and_τm_etc_from_pf = get_gens_δ_ω_ed_dash_eq_dash_and_τm_etc_from_pf( netd, bus_dict_init )

    gens_δ_ω_ed_dash_eq_dash = first.( gens_δ_ω_ed_dash_eq_dash_and_τm_etc_from_pf )

    gens_τe_Vf = second.( gens_δ_ω_ed_dash_eq_dash_and_τm_etc_from_pf )
    
    gens_id_iq_vd_vq = third.( gens_δ_ω_ed_dash_eq_dash_and_τm_etc_from_pf )

    """ convert tuples to arrays,  tuples to arrarys
    """
    
    vec_gens_δ_ω_ed_dash_eq_dash = collect( Iterators.flatten( gens_δ_ω_ed_dash_eq_dash ))
    
    vec_gens_τe_Vf = collect( Iterators.flatten( gens_τe_Vf )) 

    vec_gens_id_iq_vd_vq = collect( Iterators.flatten( gens_id_iq_vd_vq ))

    return (; vec_gens_δ_ω_ed_dash_eq_dash,
            vec_gens_τe_Vf,
            vec_gens_id_iq_vd_vq )
    
end


function make_system_matrices_with_blockdiag_matrices(
    gens_nodes_collection ;only_gen = false)
    
    """
    update is done by `update_gens_nodes_Ax_system_matrices!( vec_Ax_views, state, nodes_state_Idx, gens_nodes_collection )`
    """

    
    system_views, system_matrices, system_idxs =  create_gens_nodes_system_blockdiag_matrices_Idx_and_views( gens_nodes_collection; only_gen = only_gen )

        
    Ax_Bx_Cx_views, stab_Ax_Bx_Cx_views  = system_views
    
    Ax_Bx_Cx_matrix,  stab_Ax_Bx_Cx_matrix = system_matrices
    
    idxs, stab_idxs,  filter_stab_idxs = system_idxs

    if only_gen == false
    
        vec_Ax_views, vec_Bx_views, vec_Cx_views = Ax_Bx_Cx_views
        
        vec_stab_Ax_views, vec_stab_Bx_views, vec_stab_Cx_views = stab_Ax_Bx_Cx_views
        
        Ax_matrix, Bx_matrix, Cx_matrix = Ax_Bx_Cx_matrix
        
        stab_Ax_matrix, stab_Bx_matrix, stab_Cx_matrix = stab_Ax_Bx_Cx_matrix

        nodes_state_Idx, Bx_idxs, Cx_idxs = idxs
        
        stab_Ax_idxs, stab_Bx_idxs, stab_Cx_idxs = stab_idxs
        
        filter_stab_nodes_state_Idx, filter_stab_Bx_row_idxs, filter_stab_Cx_row_idxs = filter_stab_idxs

        init_plants_system_matrices_views!(( vec_Ax_views, vec_Bx_views, vec_Cx_views ), gens_nodes_collection; only_gen = only_gen, vec_Ax_τm_vf_views = nothing )
        
    else
        
        vec_Ax_views, vec_Bx_views, vec_Cx_views, vec_Ax_τm_vf_views = Ax_Bx_Cx_views
        
        vec_stab_Ax_views, vec_stab_Bx_views, vec_stab_Cx_views, vec_stab_τm_vf_views  = stab_Ax_Bx_Cx_views

        Ax_matrix, Bx_matrix, Cx_matrix, Ax_τm_vf_matrix = Ax_Bx_Cx_matrix
        
        stab_Ax_matrix, stab_Bx_matrix, stab_Cx_matrix, stab_Ax_τm_vf_matrix = stab_Ax_Bx_Cx_matrix

        nodes_state_Idx, Bx_idxs, Cx_idxs, τm_vf_idxs  = idxs
        
        stab_Bx_row_idxs, Bx_col_idxs = stab_idxs
        
        filter_stab_nodes_state_Idx, filter_stab_Bx_row_idxs, filter_stab_Cx_row_idxs, filter_stab_Ax_τm_vf_row_idxs = filter_stab_idxs

        init_plants_system_matrices_views!(( vec_Ax_views, vec_Bx_views, vec_Cx_views ), gens_nodes_collection; only_gen = only_gen, vec_Ax_τm_vf_views = vec_Ax_τm_vf_views )
        
    end
    
    return (; system_views, system_matrices, system_idxs)
    
end



########################################################
# ------------------------------------------------------
#  Relocated from  core.jl
# ------------------------------------------------------
########################################################


function get_gens_nodes_aggregate_system_matrices(
    vec_Ax, vec_Bx, vec_Cx,
    gens_nodes_collection;
    only_gen = false,
    vec_Ax_τm_vf = nothing  )

    no_nodes = length( gens_nodes_collection )

    nodes_dim_size, nodes_offset, nodes_state_Idx = create_flattened_type_states_dims_offset_Idx_for_gens_nodes( gens_nodes_collection )
    

    Ax_size =  sum(nodes_dim_size)

    #--------------------------------------------------
    
    # id_iq_ph_vh_vector    = [:id, :iq , :ph,  :vh]

    # id_iq_ph_vh_dims = [ length( id_iq_ph_vh_vector ) for k in 1:no_nodes ]

    # id_iq_ph_vh_offsets = create_offsets( id_iq_ph_vh_dims )
    
    # id_iq_ph_vh_idxs = create_idxs( id_iq_ph_vh_offsets,
    #                                 id_iq_ph_vh_dims )


    # ω_ref_τm_v_ref_vector = [:ω_ref, :τm,  :v_ref ]

    # ω_ref_τm_v_ref_dims = [ length( ω_ref_τm_v_ref_vector ) for k in 1:no_nodes ]

    # ω_ref_τm_v_ref_offsets = create_offsets(
    #     ω_ref_τm_v_ref_dims)
    
    # ω_ref_τm_v_ref_idxs = create_idxs(
    #     ω_ref_τm_v_ref_offsets,
    #     ω_ref_τm_v_ref_dims )
    
    #--------------------------------------------------
    
    Bx_row_dims    = [first(size(Bxi))
                      for Bxi in vec_Bx ]
    
    Bx_row_offsets = create_offsets(
        Bx_row_dims)
    
    Bx_row_idxs    = create_idxs(
        Bx_row_offsets,Bx_row_dims)

    Bx_col_dims    = [ last(size(Bxi))
                       for Bxi in vec_Bx ]
    
    Bx_col_offsets = create_offsets(Bx_col_dims)
    
    Bx_col_idxs    = create_idxs(
        Bx_col_offsets,Bx_col_dims)
    
    Bx_row_size = sum( Bx_row_dims )
    
    Bx_col_size = sum( Bx_col_dims )
    
    
    #--------------------------------------------------
    
    Cx_row_dims = [first(size(Cxi))
                   for Cxi in vec_Cx ]
    
    Cx_row_offsets = create_offsets(
        Cx_row_dims)
    
    Cx_row_idxs = create_idxs(
        Cx_row_offsets, Cx_row_dims )

    Cx_col_dims = [ last(size(Cxi))
                    for Cxi in vec_Cx ]
    
    Cx_col_offsets = create_offsets(
        Cx_col_dims)
    
    Cx_col_idxs = create_idxs(
        Cx_col_offsets, Cx_col_dims )
    
    Cx_row_size = sum(
        Cx_row_dims )
    
    Cx_col_size = sum(
        Cx_col_dims )

    #--------------------------------------------------
      
    Ax_matrix = sparse(
        1.0I, Ax_size, Ax_size )
    
    Bx_matrix = sparse(
        1.0I, Bx_row_size, Bx_row_size )
    
    Cx_matrix = sparse(
        1.0I, Cx_row_size, Cx_col_size )

    #--------------------------------------------------

    for (i, Ax_node_matrix ) in enumerate( vec_Ax )
        
        node_state_Idx = nodes_state_Idx[i]
        
        copyto!( @view( Ax_matrix[
            node_state_Idx, node_state_Idx ]),
                 Ax_node_matrix )
        
    end


    for (Bx_row_idx, Bx_col_idx, Bx_node_matrix ) in
        zip( Bx_row_idxs, Bx_col_idxs, vec_Bx )
        
        copyto!( @view( Bx_matrix[
            Bx_row_idx, Bx_col_idx ]),
                 Bx_node_matrix )
        
    end

    for (Cx_row_idx, Cx_col_idx, Cx_node_matrix ) in
        zip( Cx_row_idxs, Cx_col_idxs, vec_Cx )
        
        copyto!( @view( Cx_matrix[
            Cx_row_idx, Cx_col_idx ]),
                 Cx_node_matrix )
        
    end

    #--------------------------------------------------

    if only_gen == true

        τm_vf_vector = [:τm, :vf ]

        τm_vf_dims = [
            length( τm_vf_vector )
            for k in 1:no_nodes ]

        τm_vf_offsets = create_offsets(
            τm_vf_dims )

        τm_vf_idxs = create_idxs(
            τm_vf_offsets, τm_vf_dims )

        τm_vf_size =  sum( τm_vf_dims )

        Ax_τm_vf_row_dims = [
            first(size(Ax_τm_vfi))
            for Ax_τm_vfi in Ax_τm_vf ]
        
        Ax_τm_vf_row_offsets = create_offsets(
            Ax_τm_vf_row_dims)
        
        Ax_τm_vf_row_idxs = create_idxs(
            Ax_τm_vf_row_offsets,
            Ax_τm_vf_row_dims )

        Ax_τm_vf_col_dims = [
            last(size(Ax_τm_vfi)) for
                Ax_τm_vfi in vec_Ax_τm_vf ]
        
        Ax_τm_vf_col_offsets = create_offsets(
            Ax_τm_vf_col_dims)
        
        Ax_τm_vf_col_idxs = create_idxs(
            Ax_τm_vf_col_offsets,
            Ax_τm_vf_col_dims )

        Ax_τm_vf_row_size = sum(
            Ax_τm_vf_row_dims )
        
        Ax_τm_vf_col_size = sum(
            Ax_τm_vf_col_dims )

        Ax_τm_vf_matrix = sparse(
            1.0I,
            Ax_τm_vf_row_size,
            Ax_τm_vf_col_size )

        for (Ax_τm_vf_row_idx, Ax_τm_vf_col_idx, Ax_τm_vf_node_matrix ) in zip( Ax_τm_vf_row_idxs, Ax_τm_vf_col_idxs, vec_Ax_τm_vf )
            copyto!( @view( Ax_τm_vf_matrix[
                Ax_τm_vf_row_idx, Ax_τm_vf_col_idx ]),
                     Ax_τm_vf_node_matrix )

        end
        

    end
    
    #--------------------------------------------------

    if only_gen == false

        return Ax_matrix, Bx_matrix, Cx_matrix
        
    else

        return (Ax_matrix,
                Bx_matrix,
                Cx_matrix,
                τm_vf_matrix)
    end
    
        
end


function get_gens_nodes_aggregate_system_stability_matrices(
    vec_Ax,
    vec_Bx,
    vec_Cx,
    gens_nodes_collection;
    only_gen = false,
    vec_stab_Ax_τm_vf = nothing  )

    no_nodes = length( gens_nodes_collection )

    (nodes_dim_size,
     nodes_offset,
     nodes_state_Idx) = create_flattened_type_states_dims_offset_Idx_for_gens_nodes( gens_nodes_collection )
    

    Ax_size =  sum(nodes_dim_size)

    #--------------------------------------------------
    
    # id_iq_ph_vh_vector    = [:id, :iq , :ph,  :vh]

    # id_iq_ph_vh_dims = [ length( id_iq_ph_vh_vector ) for k in 1:no_nodes ]

    # id_iq_ph_vh_offsets = create_offsets( id_iq_ph_vh_dims )
    
    # id_iq_ph_vh_idxs = create_idxs( id_iq_ph_vh_offsets,
    #                                 id_iq_ph_vh_dims )


    # ω_ref_τm_v_ref_vector = [:ω_ref, :τm,  :v_ref ]

    # ω_ref_τm_v_ref_dims = [ length( ω_ref_τm_v_ref_vector ) for k in 1:no_nodes ]

    # ω_ref_τm_v_ref_offsets = create_offsets(
    #     ω_ref_τm_v_ref_dims)
    
    # ω_ref_τm_v_ref_idxs = create_idxs(
    #     ω_ref_τm_v_ref_offsets,
    #     ω_ref_τm_v_ref_dims )
    
    #--------------------------------------------------
    
    Bx_row_dims    = [first(size(Bxi))
                      for Bxi in vec_Bx ]
    
    Bx_row_offsets = create_offsets(
        Bx_row_dims)
    
    Bx_row_idxs    = create_idxs(
        Bx_row_offsets,Bx_row_dims)

    Bx_col_dims    = [ last(size(Bxi))
                       for Bxi in vec_Bx ]
    
    Bx_col_offsets = create_offsets(
        Bx_col_dims)
    
    Bx_col_idxs    = create_idxs(
        Bx_col_offsets,Bx_col_dims)
    
    Bx_row_size = sum( Bx_row_dims )
    
    Bx_col_size = sum( Bx_col_dims )
    
    
    #--------------------------------------------------
    
    Cx_row_dims = [first(size(Cxi))
                   for Cxi in vec_Cx ]
    
    Cx_row_offsets = create_offsets(
        Cx_row_dims)
    
    Cx_row_idxs = create_idxs(
        Cx_row_offsets, Cx_row_dims )

    Cx_col_dims = [ last(size(Cxi))
                    for Cxi in vec_Cx ]
    
    Cx_col_offsets = create_offsets(
        Cx_col_dims)
    
    Cx_col_idxs = create_idxs(
        Cx_col_offsets, Cx_col_dims )
    
    Cx_row_size = sum( Cx_row_dims )
    
    Cx_col_size = sum( Cx_col_dims )

    #--------------------------------------------------
      
    Ax_matrix = sparse( 1.0I, Ax_size, Ax_size )
    
    Bx_matrix = sparse( 1.0I, Bx_row_size, Bx_row_size )
    
    Cx_matrix = sparse( 1.0I, Cx_row_size, Cx_col_size )

    #--------------------------------------------------

    for (i, Ax_node_matrix ) in enumerate( vec_Ax )
        
        node_state_Idx = nodes_state_Idx[i]
        
        copyto!( @view( Ax_matrix[
            node_state_Idx, node_state_Idx ]),
                 Ax_node_matrix )
        
    end


    for (Bx_row_idx, Bx_col_idx, Bx_node_matrix ) in
        zip( Bx_row_idxs, Bx_col_idxs, vec_Bx )
        
        copyto!( @view( Bx_matrix[
            Bx_row_idx, Bx_col_idx ]),
                 Bx_node_matrix )
        
    end

    for (Cx_row_idx, Cx_col_idx, Cx_node_matrix ) in
        zip( Cx_row_idxs, Cx_col_idxs, vec_Cx )
        
        copyto!( @view( Cx_matrix[
            Cx_row_idx, Cx_col_idx ]),
                 Cx_node_matrix )
        
    end

    #--------------------------------------------------

    if only_gen == true

        τm_vf_vector = [:τm, :vf ]

        τm_vf_dims = [
            length( τm_vf_vector )
            for k in 1:no_nodes ]

        τm_vf_offsets = create_offsets(
            τm_vf_dims )

        τm_vf_idxs = create_idxs(
            τm_vf_offsets, τm_vf_dims )

        τm_vf_size =  sum( τm_vf_dims )

        Ax_τm_vf_row_dims = [
            first(size(Ax_τm_vfi))
            for Ax_τm_vfi in Ax_τm_vf ]
        
        Ax_τm_vf_row_offsets = create_offsets(
            Ax_τm_vf_row_dims)
        
        Ax_τm_vf_row_idxs = create_idxs(
            Ax_τm_vf_row_offsets, Ax_τm_vf_row_dims )

        Ax_τm_vf_col_dims = [
            last(size(Ax_τm_vfi))
            for Ax_τm_vfi in
                vec_stab_Ax_τm_vf ]
        
        Ax_τm_vf_col_offsets = create_offsets(
            Ax_τm_vf_col_dims)
        
        Ax_τm_vf_col_idxs = create_idxs(
            Ax_τm_vf_col_offsets,
            Ax_τm_vf_col_dims )

        Ax_τm_vf_row_size = sum( Ax_τm_vf_row_dims )
        
        Ax_τm_vf_col_size = sum( Ax_τm_vf_col_dims )

        Ax_τm_vf_matrix = sparse(
            1.0I, Ax_τm_vf_row_size,
            Ax_τm_vf_col_size )

        for (Ax_τm_vf_row_idx, Ax_τm_vf_col_idx, Ax_τm_vf_node_matrix ) in zip( Ax_τm_vf_row_idxs, Ax_τm_vf_col_idxs, vec_stab_Ax_τm_vf )
            
            copyto!( @view( Ax_τm_vf_matrix[
                Ax_τm_vf_row_idx, Ax_τm_vf_col_idx ]),
                     Ax_τm_vf_node_matrix )

        end
        

    end
    
    #--------------------------------------------------

    if only_gen == false

        return (Ax_matrix,
                Bx_matrix,
                Cx_matrix)
        
    else

        return (Ax_matrix,
                Bx_matrix,
                Cx_matrix,
                τm_vf_matrix)
    end
    
        
end


#--------------------------------------------------
#--------------------------------------------------


function update_a_plant_Ax_system_matrices!(
    Ax, state, plant;
    only_gen =false  )

    if no_control_device == false
        
        if only_gen_avr_in_plant(plant) || gen_gov_avr_in_plant(plant)
            dict_state_vars_syms = plant.dict_state_vars_syms

            vr1_idx = dict_state_vars_syms[:vr1]
            
            vf_tilade_idx = dict_state_vars_syms[:vf_tilade]

            # Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be = V_R_max = plant.Exc.param_values

            (Ta, Te, Tf, Tr, Ka, Ke, Kf,
             V_R_max, V_R_min, Ae, Be) =
                plant.Exc.param_values

            #  Be = V_R_max


            # It seems to me later that is not necessary to
            # update α, i.e the column of vr1 in Ax

            # vr1  = x[dict_state_vars_syms[:vr1]]
            # α = threshold_limits(vr1,V_R_max,V_R_min)/(vr1 * Te)


            vf_tilade = state[ dict_state_vars_syms[:vf_tilade] ]


            α =  1/ Te

            γ = -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te

            # Ax[vr1_idx, vr1_idx]             = -α
            # Ax[vf_tilade_idx, vr1_idx]       = α
            
            Ax[vf_tilade_idx, vf_tilade_idx] = γ

        end
    else
        nothing
    end
    

    return nothing        
end



function update_a_plant_system_matrices!(
    Ax, x, plant;
    only_gen = false )

    if only_gen == false
        
        if only_gen_avr_in_plant(plant) || gen_gov_avr_in_plant(plant)
            dict_state_vars_syms = plant.dict_state_vars_syms

            vr1_idx       = dict_state_vars_syms[:vr1]
            vf_tilade_idx = dict_state_vars_syms[:vf_tilade]

            # Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be = V_R_max = plant.Exc.param_values

            (Ta, Te, Tf, Tr, Ka, Ke, Kf,
             V_R_max, V_R_min,
             Ae, Be) =
                plant.Exc.param_values
            
            # Be = V_R_max

            vf_tilade = x[ dict_state_vars_syms[:vf_tilade] ]

            # It seems to me later that is not necessary to
            # update α, i.e the column of vr1 in Ax

            # vr1  = x[dict_state_vars_syms[:vr1]]
            # α    =  threshold_limits(vr1, V_R_max, V_R_min)/(vr1 * Te)

            α =  1/ Te

            γ = -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te

            # Ax[vr1_idx, vr1_idx]             = -α
            # Ax[vf_tilade_idx, vr1_idx]       = α
            
            Ax[vf_tilade_idx, vf_tilade_idx] = γ

        end
    else
        nothing
    end
    

    return nothing        
end


function update_a_plant_system_matrices_views!(
    (Ax_view, Bx_view, Cx_view), plant;
    τm_vf_view = nothing,
    only_gen = false )

    if only_gen == false

        plant.func_system_matrices[1](
            (Ax_view, Bx_view, Cx_view), plant)
    else
        plant.func_system_matrices[1](
            (Ax_view, Bx_view, Cx_view), plant;
            only_gen=only_gen, τm_vf_view = τm_vf_view )
    end

    return nothing
    
end

"""
This function initialises system matrices views
` Ax_view, Bx_view, Cx_view` which are in vectors
`vec_Ax_views, vec_Bx_views, vec_Cx_views`.

These views are created by
`create_gens_nodes_aggregate_system_matrices_Idx_and_views`.

The intialisation is expected to be reflected in aggregate
system matrices `Ax_matrix, Bx_matrix, Cx_matrix `

"""
function init_plants_system_matrices_views!(
    (vec_Ax_views,
     vec_Bx_views,
     vec_Cx_views),
    gens_nodes_collection;
    vec_Ax_τm_vf_views = nothing,
    only_gen = false )

    # update_a_plant_system_matrices_views!

    if only_gen == false
        for (Ax_view, Bx_view, Cx_view, plant) in
            zip(vec_Ax_views, vec_Bx_views,
                vec_Cx_views, gens_nodes_collection)
            
            update_a_plant_system_matrices_views!(
                (Ax_view, Bx_view, Cx_view), plant)
        end
    else
        for (Ax_view, Bx_view,Cx_view,τm_vf_view, plant) in
            zip(vec_Ax_views, vec_Bx_views, vec_Cx_views,
                vec_Ax_τm_vf_views, gens_nodes_collection)
            
            update_a_plant_system_matrices_views!(
                (Ax_view, Bx_view, Cx_view), plant;
                only_gen = only_gen,
                τm_vf_view = τm_vf_view )
        end
    end
  
    return nothing        
end


#--------------------------------------------------


function update_a_plant_system_stability_matrices!(
    stab_Ax, x, plant ;
    only_gen = false )

    if only_gen == false
        
        if only_gen_avr_in_plant(plant) || gen_gov_avr_in_plant(plant)
            dict_stab_state_syms = plant.dict_stab_state_syms

            vr1_idx       = dict_stab_state_syms[:vr1]
            vf_tilade_idx = dict_stab_state_syms[:vf_tilade]

            # Ta, Te, Tf, Tr, Ka, Ke, Kf, V_R_max, V_R_min, Ae, Be = V_R_max = plant.Exc.param_values

            (Ta, Te, Tf, Tr, Ka, Ke, Kf,
             V_R_max, V_R_min,
             Ae, Be) =
                 plant.Exc.param_values
            
            # Be = V_R_max

            vf_tilade = x[dict_stab_state_syms[:vf_tilade]]

            # It seems to me later that is not necessary to
            # update α, i.e the column of vr1 in Ax

            # vr1  = x[dict_state_vars_syms[:vr1]]
            # α    =  threshold_limits(vr1, V_R_max, V_R_min)/(vr1 * Te)

            α =  1/ Te

            γ = -1*(Ke + Sevf(Ae, Be, vf_tilade))/Te

            # stab_Ax[vr1_idx, vr1_idx]             = -α
            # stab_Ax[vf_tilade_idx, vr1_idx]       = α
            
            stab_Ax[vf_tilade_idx, vf_tilade_idx] = γ

        end
    else
        nothing
    end
    

    return nothing        
end


function update_a_plant_system_stability_matrices_views!(
    (stab_Ax_view,  stab_Bx_view,  stab_Cx_view),
    plant;
    only_gen = false,
    stab_τm_vf_view = nothing )

    if only_gen == false

        plant.func_system_stab_matrices[1](
            ( stab_Ax_view,
              stab_Bx_view,
              stab_Cx_view), plant)
    else
        plant.func_system_stab_matrices[1](
            ( stab_Ax_view,
              stab_Bx_view,
              stab_Cx_view), plant;
            only_gen=true,
            stab_τm_vf_view =
                stab_τm_vf_view )
    end

    return nothing
    
end

"""
This function initialises system stability matrices views
` stab_Ax_view, stab_Bx_view, stab_Cx_view` which are in vectors
`vec_stab_Ax_views, vec_stab_Bx_views, vec_stab_Cx_views`.

These views are created by
`create_gens_nodes_aggregate_system_stability_matrices_Idx_and_views`.

The intialisation is expected to be reflected in aggregate
system matrices `stab_Ax_matrix, stab_Bx_matrix, stab_Cx_matrix `
"""
function init_plants_system_stability_matrices_views!(
    ( vec_stab_Ax_views,
      vec_stab_Bx_views,
      vec_stab_Cx_views ),
    gens_nodes_collection;
    only_gen = false,
    vec_stab_τm_vf_views = nothing )

    # update_a_plant_system_matrices_views!

    if only_gen == false
        for (stab_Ax_view,stab_Bx_view,stab_Cx_view,plant) in
            zip( vec_stab_Ax_views, vec_stab_Bx_views,
                vec_stab_Cx_views, gens_nodes_collection )
            
            update_a_plant_system_stability_matrices_views!(
                ( stab_Ax_view,
                  stab_Bx_view,
                  stab_Cx_view ), plant)
        end
    else
        for ( stab_Ax_view, stab_Bx_view, stab_Cx_view, stab_τm_vf_view, plant ) in
            zip(
                vec_stab_Ax_views,
                vec_stab_Bx_views,
                vec_stab_Cx_views,
                vec_stab_τm_vf_views,
                gens_nodes_collection)
            
            update_a_plant_system_stability_matrices_views!(
                ( stab_Ax_view, stab_Bx_view, stab_Cx_view ),
                plant;
                only_gen = true,
                stab_τm_vf_view = stab_τm_vf_view )
        end
    end
  
    return nothing        
end



function update_gens_nodes_Ax_system_matrices!(    
    vec_Ax,
    state,
    nodes_state_Idx,
    gens_nodes_collection;
    only_gen = false )

    # nodes_x_view  = [view(x, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    if only_gen == false
        no_nodes = length( vec_Ax  )

        for (Ax_node_matrix, node_state_Idx, gen_node) in
            zip(
            vec_Ax,
            nodes_state_Idx,
            gens_nodes_collection )

            update_a_plant_system_matrices!(
                Ax_node_matrix,
                state[node_state_Idx],
                gen_node )
        end
    else

        no_nodes = length( vec_Ax  )

        for (Ax_node_matrix, node_state_Idx, gen_node) in
            zip(
            vec_Ax,
            nodes_state_Idx,
            gens_nodes_collection )

            update_a_plant_system_matrices!(
                Ax_node_matrix,
                state[node_state_Idx],
                gen_node,
                only_gen = only_gen )
        end                
        
    end
    

    #--------------------------------------------------

    return nothing
end


function update_gens_nodes_system_matrices!(    
    vec_Ax, x,
    nodes_state_Idx,
    gens_nodes_collection;
    only_gen = false )

    # nodes_x_view  = [view(x, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    if only_gen == false
        
        no_nodes = length( vec_Ax  )

        for (Ax_node_matrix, node_state_Idx, gen_node) in
            zip(
            vec_Ax,
            nodes_state_Idx,
            gens_nodes_collection )

            update_a_plant_system_matrices!(
                Ax_node_matrix,
                x[node_state_Idx],
                gen_node )
        end
    else
        no_nodes = length( vec_Ax  )

        for (Ax_node_matrix, node_state_Idx, gen_node) in
            zip(
            vec_Ax,
            nodes_state_Idx,
            gens_nodes_collection )

            update_a_plant_system_matrices!(
                Ax_node_matrix,
                x[node_state_Idx],
                gen_node;
                only_gen = only_gen )            
        end
        
    end
        
    #--------------------------------------------------

    return nothing
end


function update_gens_nodes_system_stability_matrices!(    
    vec_Ax,
    x,
    nodes_state_Idx,
    gens_nodes_collection;
    only_gen = false )

    # nodes_x_view  = [view(x, nodes_state_Idx[Ind]) for Ind in collect(1:length(nodes_state_Idx))]

    if only_gen == false
        
        no_nodes = length( vec_Ax  )

        for (Ax_node_matrix, node_state_Idx, gen_node) in
            zip(
            vec_Ax, nodes_state_Idx, gens_nodes_collection )

            update_a_plant_system_stability_matrices!(
                Ax_node_matrix,
                x[node_state_Idx],
                gen_node )
        end
    else

        
        no_nodes = length( vec_Ax  )

        for (Ax_node_matrix, node_state_Idx, gen_node) in
            zip(
            vec_Ax,
            nodes_state_Idx,
            gens_nodes_collection )

            update_a_plant_system_stability_matrices!(
                Ax_node_matrix,
                x[node_state_Idx],
                gen_node;
                only_gen =only_gen )
        end
    end
    
    #--------------------------------------------------

    return nothing
end


function update_gens_nodes_linked_system_matrices!(
    system_views,
    system_matrices,
    system_idxs,
    state,
    gens_nodes_collection;
    only_gen = false )


    (Ax_Bx_Cx_views,
     stab_Ax_Bx_Cx_views,
     stab_link_Ax_Bx_Cx_views) =
         system_views
    
    (Ax_Bx_Cx_matrix,
     stab_Ax_Bx_Cx_matrix) =
         system_matrices
    
    (idxs,
     stab_idxs,
     filter_stab_idxs) =
         system_idxs


    if only_gen == false
    
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views) =
             Ax_Bx_Cx_views
        
        (vec_stab_Ax_views,
         vec_stab_Bx_views,
         vec_stab_Cx_views) =
             stab_Ax_Bx_Cx_views
        
        (vec_stab_link_Ax_views,
         vec_stab_link_Bx_views,
         vec_stab_link_Cx_views) =
             stab_link_Ax_Bx_Cx_views

        (Ax_matrix,
         Bx_matrix,
         Cx_matrix) =
             Ax_Bx_Cx_matrix
        
        (stab_Ax_matrix,
         stab_Bx_matrix,
         stab_Cx_matrix) =
             stab_Ax_Bx_Cx_matrix

        (nodes_state_Idx,
         Bx_idxs,
         Cx_idxs) =
             idxs
        
        (stab_Ax_idxs,
         stab_Bx_idxs,
         stab_Cx_idxs) =
             stab_idxs
        
        (filter_stab_nodes_state_Idx,
         filter_stab_Bx_row_idxs,
         filter_stab_Cx_row_idxs) =
             filter_stab_idxs

        (Bx_row_idxs,
         Bx_col_idxs) =
             Bx_idxs
        
        (Cx_row_idxs,
         Cx_col_idxs) =
             Cx_idxs

        (stab_Bx_row_idxs,
         Bx_col_idxs) =
             stab_Bx_idxs
        
        (stab_Cx_row_idxs,
         Cx_col_idxs) =
             stab_Cx_idxs

        """
        update vec_Ax_views, this will automatically
        update Ax_matrix
        """
        
        update_gens_nodes_Ax_system_matrices!(
            vec_Ax_views,
            state,
            nodes_state_Idx,
            gens_nodes_collection )

        """
        copy the contents of vec_stab_Ax_views,
        vec_stab_Bx_views
        vec_stab_Cx_views to vec_stab_link_Ax_views,
        vec_stab_link_Bx_views, vec_stab_link_Cx_views
        respectively. This will automatically update
        stab_Ax_matrix, stab_Bx_matrix and stab_Cx_matrix
        """
        
        for (stab_link_Ax_view, stab_Ax_view ) in
            zip(vec_stab_link_Ax_views, vec_stab_Ax_views )

            stab_link_Ax_view[:, :] .=
                stab_Ax_view[:,:]
        end

        for ( stab_link_Bx_view, stab_Bx_view ) in
            zip(vec_stab_link_Bx_views, vec_stab_Bx_views )

            stab_link_Bx_view[:, :] .=
                stab_Bx_view[:,:]
        end
        
        for ( stab_link_Cx_view, stab_Cx_view ) in
            zip(vec_stab_link_Cx_views, vec_stab_Cx_views )

            stab_link_Cx_view[:, :] .=
                stab_Cx_view[:,:]
        end
        
    else
        
        (vec_Ax_views,
         vec_Bx_views,
         vec_Cx_views,
         vec_Ax_τm_vf_views) =
             Ax_Bx_Cx_views
        
        (vec_stab_Ax_views,
         vec_stab_Bx_views,
         vec_stab_Cx_views,
         vec_stab_τm_vf_views) =
             stab_Ax_Bx_Cx_views
        
        (vec_stab_link_Ax_views,
         vec_stab_link_Bx_views,
         vec_stab_link_Cx_views,
         vec_stab_link_τm_vf_views) =
             stab_link_Ax_Bx_Cx_views

        (Ax_matrix,
         Bx_matrix,
         Cx_matrix,
         Ax_τm_vf_matrix) =
             Ax_Bx_Cx_matrix
        
        (stab_Ax_matrix,
         stab_Bx_matrix,
         stab_Cx_matrix,
         stab_Ax_τm_vf_matrix) =
             stab_Ax_Bx_Cx_matrix

        (nodes_state_Idx,
         Bx_idxs,
         Cx_idxs,
         τm_vf_idxs) =
             idxs
        
        (stab_Ax_idxs,
         stab_Bx_row_idxs,
         Bx_col_idxs,
         stab_τm_vf_idxs) =
             stab_idxs
        
        (filter_stab_nodes_state_Idx,
         filter_stab_Bx_row_idxs,
         filter_stab_Cx_row_idxs,
         filter_stab_Ax_τm_vf_row_idxs) =
             filter_stab_idxs

        (Bx_row_idxs,
         Bx_col_idxs) =
             Bx_idxs
        
        (Cx_row_idxs,
         Cx_col_idxs) =
             Cx_idxs
        
        (Ax_τm_vf_row_idxs,
         Ax_τm_vf_col_idxs) =
             τm_vf_idxs

        (stab_Bx_row_idxs,
         Bx_col_idxs) =
             stab_Bx_idxs
        
        (stab_Cx_row_idxs,
         Cx_col_idxs) =
             stab_Cx_idxs

        (stab_Ax_τm_vf_row_idxs,
         Ax_τm_vf_col_idxs) =
             stab_τm_vf_idxs

        """
        update vec_Ax_views, this will automatically
        update Ax_matrix
        """
        
        update_gens_nodes_Ax_system_matrices!(
            vec_Ax_views,
            state,
            nodes_state_Idx,
            gens_nodes_collection )

        """
        copy the contents of vec_stab_Ax_views,
        vec_stab_Bx_views
        vec_stab_Cx_views, and vec_stab_link_τm_vf_views to
        vec_stab_link_Ax_views, vec_stab_link_Bx_views,
        vec_stab_link_Cx_views and vec_stab_τm_vf_views
        respectively.
        
        This will automatically update stab_Ax_matrix,
        stab_Bx_matrix, stab_Cx_matrix and Ax_τm_vf_matrix
        """        
        
        for ( stab_link_Ax_view, stab_Ax_view ) in
            zip( vec_stab_link_Ax_views,
                 vec_stab_Ax_views )

            stab_link_Ax_view[:, :] .=
                stab_Ax_view[:,:]
        end

        for ( stab_link_Bx_view, stab_Bx_view ) in
            zip( vec_stab_link_Bx_views,
                 vec_stab_Bx_views )

            stab_link_Bx_view[:, :] .=
                stab_Bx_view[:,:]
        end
        
        for ( stab_link_Cx_view, stab_Cx_view ) in
            zip( vec_stab_link_Cx_views,
                 vec_stab_Cx_views )

            stab_link_Cx_view[:, :] .=
                stab_Cx_view[:,:]
        end

        
        for ( stab_link_τm_vf_view, stab_τm_vf_view ) in
            zip( vec_stab_link_τm_vf_views,
                 vec_stab_τm_vf_views )

            stab_link_τm_vf_view[:, :] .=
                stab_τm_vf_view[:,:]
        end        
        
    end
    
    #--------------------------------------------------

    return nothing
end



#--------------------------------------------------
###################################################
#--------------------------------------------------

#--------------------------------------------------------
# Diagnosis
#--------------------------------------------------------

# M-x flush-lines "^$"

# To resent an emacs register

# (set-register ?a nil)

#--------------------------------------------------------
# Diagnosis
#--------------------------------------------------------

"""
https://www.masteringemacs.org/article/removing-blank-lines-buffer
https://whatacold.io/blog/2023-06-12-emacs-join-lines/

"""


"""

 SC = ( "bus3", "bus6", "bus8" )

 ω_ref0 = ωs

 plant1 = netd.nodes["bus1"]

 plant3 = netd.nodes["bus3"]


t_Ax_gen, t_Bx_gen, t_Cx_gen, t_Ax_gen_S_avr_S, t_Ax_gen_S_avr_A , t_avr_Ax_S_S, t_avr_Ax_S_Sgen, t_avr_Bx_S, t_avr_Cx_S, t_avr_Ax_S_A, t_avr_Ax_A_S, t_avr_Ax_A_Sgen, t_avr_Ax_A_A, t_avr_Bx_A, t_avr_Cx_A, t_gen_m, t_avr_S_m, t_avr_A_m, t_Ax, t_Bx, t_Cx  =  diagnosis_get_only_gen_avr_im_plant_system_matrices( plant3; ref_in_state = false  )


"""

function diagnosis_get_only_gen_avr_im_plant_system_matrices(
    Gen, Exc, dummy_vr1, dummy_vf_tilade;
    ref_in_state = false )

    Ax_gen = Gen.Ax(
        Gen,
        Gen.param_matrix_values... )
    
    Bx_gen = Gen.Bx(
        Gen,
        Gen.param_matrix_values...)
    
    Cx_gen = Gen.Cx(
        Gen,
        Gen.param_matrix_values...)

    Ax_gen_S_avr_S = Gen.Ax_gen_S_avr_S(
        Gen,
        Exc )

    Ax_gen_S_avr_A = Gen.Ax_gen_S_avr_A(
        Gen,
        Exc )

    #---------------------------------------

    avr_Ax_S_S = Exc.Ax_S_S(
        Exc,
        dummy_vf_tilade,
        dummy_vr1)

    avr_Ax_S_Sgen = Exc.Ax_S_Sgen( Exc )

    avr_Bx_S = Exc.Bx_S( Exc )

    avr_Cx_S = Exc.Cx_S( Exc )

    #---------------------------------------

    avr_Ax_S_A = Exc.Ax_S_A(
        Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_S = Exc.Ax_A_S(
        Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_Sgen = Exc.Ax_A_Sgen(
        Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_A  = Exc.Ax_A_A(
        Exc;
        ref_in_state = ref_in_state )

    avr_Bx_A = Exc.Bx_A(
        Exc;
        ref_in_state = ref_in_state )

    avr_Cx_A = Exc.Cx_A(
        Exc;
        ref_in_state = ref_in_state )

    gen_m = hcat(Ax_gen, Ax_gen_S_avr_S, Ax_gen_S_avr_A )

    avr_S_m = hcat(avr_Ax_S_Sgen, avr_Ax_S_S,  avr_Ax_S_A )

    avr_A_m =hcat(avr_Ax_A_Sgen, avr_Ax_A_S, avr_Ax_A_A )

    Ax = vcat(gen_m,
              avr_S_m,
              avr_A_m)              

    Bx = vcat(Bx_gen,
              avr_Bx_S,
              avr_Bx_A)

    Cx = vcat(Cx_gen,
              avr_Cx_S,
              avr_Cx_A )

    return (; Ax_gen,Bx_gen,Cx_gen,
            Ax_gen_S_avr_S,
            Ax_gen_S_avr_A, avr_Ax_S_S,
            avr_Ax_S_Sgen, avr_Bx_S, avr_Cx_S,
            avr_Ax_S_A, avr_Ax_A_S, avr_Ax_A_Sgen,
            avr_Ax_A_A, avr_Bx_A, avr_Cx_A,
            gen_m, avr_S_m, avr_A_m,
            Ax, Bx, Cx  )                
    
end

function diagnosis_get_only_gen_avr_im_plant_system_matrices(
    plant; ref_in_state = false )

    return diagnosis_get_only_gen_avr_im_plant_system_matrices(
    plant.Gen, plant.Exc, plant.dummy_vr1, plant.dummy_vf_tilade;
    ref_in_state = ref_in_state )                
            
end






"""

 SM = ( "bus1", "bus2" )

 ω_ref0 = ωs

 plant1 = netd.nodes["bus1"]


diag_gen_gov_avr_plant = diagnosis_get_gen_gov_avr_im_plant_system_matrices(  plant1.Gen, plant1.Gov, plant1.Exc, plant1.dummy_vr1, plant1.dummy_vf_tilade; ω_ref0 = ωs, ref_in_state = false )


tt_Ax_gen,tt_Bx_gen,tt_Cx_gen, tt_Ax_gen_S_gov_S, tt_Ax_gen_S_gov_A, tt_Ax_gen_S_avr_S, tt_Ax_gen_S_avr_A, tt_gov_Ax_S_S, tt_gov_Ax_S_Sgen, tt_gov_Bx_S, tt_gov_Cx_S, tt_avr_Ax_S_S, tt_avr_Ax_S_Sgen, tt_avr_Bx_S, tt_avr_Cx_S, tt_gov_Ax_S_A, tt_gov_Ax_A_S, tt_gov_Ax_A_A, tt_gov_Ax_A_Sgen, tt_avr_Ax_S_A, tt_avr_Ax_A_S, tt_avr_Ax_A_A, tt_avr_Ax_A_Sgen, tt_z_gov_S_avr_A, tt_z_gov_S_avr_S, tt_z_gov_A_avr_S, tt_z_gov_A_avr_A, tt_z_avr_S_gov_A, tt_z_avr_S_gov_S, tt_z_avr_A_gov_S, tt_z_avr_A_gov_A, tt_gov_Bx_A, tt_gov_Cx_A, tt_avr_Bx_A, tt_avr_Cx_A, tt_gen_m, tt_gov_S_m, tt_gov_A_m, tt_avr_S_m, tt_avr_A_m, tt_Ax, tt_Bx, tt_Cx = diag_gen_gov_avr_plant


tt_Ax_gen
tt_Bx_gen
tt_Cx_gen
tt_Ax_gen_S_gov_S
tt_Ax_gen_S_gov_A
tt_Ax_gen_S_avr_S
tt_Ax_gen_S_avr_A

tt_gov_Ax_S_S
tt_gov_Ax_S_Sgen
tt_gov_Bx_S
tt_gov_Cx_S

tt_avr_Ax_S_S
tt_avr_Ax_S_Sgen
tt_avr_Bx_S
tt_avr_Cx_S

tt_gov_Ax_S_A
tt_gov_Ax_A_S
tt_gov_Ax_A_A
tt_gov_Ax_A_Sgen

tt_avr_Ax_S_A
tt_avr_Ax_A_S
tt_avr_Ax_A_A
tt_avr_Ax_A_Sgen

tt_z_gov_S_avr_A
tt_z_gov_S_avr_S
tt_z_gov_A_avr_S
tt_z_gov_A_avr_A

tt_z_avr_S_gov_A
tt_z_avr_S_gov_S
tt_z_avr_A_gov_S
tt_z_avr_A_gov_A

tt_gov_Bx_A
tt_gov_Cx_A
tt_avr_Bx_A
tt_avr_Cx_A

tt_gen_m

tt_gov_S_m
tt_gov_A_m

tt_avr_S_m
tt_avr_A_m

tt_Ax
tt_Bx
tt_Cx 


    """


function diagnosis_get_gen_gov_avr_im_plant_system_matrices(
    Gen, Gov, Exc, dummy_vr1, dummy_vf_tilade;
    ω_ref0 = ωs, ref_in_state = false )

    #--------------------------------
    
    Ax_gen = Gen.Ax(
        Gen,
        Gen.param_matrix_values...)
    
    Bx_gen = Gen.Bx(
        Gen,
        Gen.param_matrix_values...)
    
    Cx_gen = Gen.Cx(
        Gen,
        Gen.param_matrix_values...)

    #---------------------------------------

    Ax_gen_S_gov_S = Gen.Ax_gen_S_gov_S(
        Gen, Gov )

    Ax_gen_S_gov_A = Gen.Ax_gen_S_gov_A(
        Gen, Gov )

    #---------------------------------------
    
    Ax_gen_S_avr_S = Gen.Ax_gen_S_avr_S(
        Gen, Exc )

    Ax_gen_S_avr_A = Gen.Ax_gen_S_avr_A(
        Gen, Exc )

    #---------------------------------------

    gov_Ax_S_S    = Gov.Ax_S_S( Gov )

    gov_Ax_S_Sgen = Gov.Ax_S_Sgen( Gov )

    gov_Bx_S      = Gov.Bx_S( Gov )

    gov_Cx_S      = Gov.Cx_S( Gov )

    #--------------------------------

    avr_Ax_S_S = Exc.Ax_S_S(
        Exc,
        dummy_vf_tilade,
        dummy_vr1)

    avr_Ax_S_Sgen = Exc.Ax_S_Sgen( Exc )

    avr_Bx_S = Exc.Bx_S( Exc )

    avr_Cx_S = Exc.Cx_S( Exc )

    #--------------------------------

    gov_Ax_S_A = Gov.Ax_S_A(
        Gov; ω_ref0 = ωs,
        ref_in_state = ref_in_state  )

    gov_Ax_A_S = Gov.Ax_A_S(
        Gov; ω_ref0 = ωs,
        ref_in_state = ref_in_state)

    gov_Ax_A_A  = Gov.Ax_A_A(
        Gov; ω_ref0 = ωs,
        ref_in_state = ref_in_state)

    gov_Ax_A_Sgen = Gov.Ax_A_Sgen(
        Gov; ω_ref0 = ωs,
        ref_in_state = ref_in_state)

    #--------------------------------

    avr_Ax_S_A = Exc.Ax_S_A(
        Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_S = Exc.Ax_A_S(
        Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_A = Exc.Ax_A_A(
        Exc;
        ref_in_state = ref_in_state )

    avr_Ax_A_Sgen = Exc.Ax_A_Sgen(
        Exc;
        ref_in_state = ref_in_state )
    
    #--------------------------------

    z_gov_S_avr_A =
        zero_coupling_state_vars_to_im_alg_vars(
        Gov, Exc)

    z_gov_S_avr_S =
        zero_coupling_state_vars_to_state_vars(
        Gov, Exc )

    z_gov_A_avr_S =
        zero_coupling_im_alg_vars_to_state_vars(
        Gov, Exc )

    z_gov_A_avr_A =
        zero_coupling_im_alg_vars_to_im_alg_vars(
        Gov, Exc )

    #--------------------------------
    
    z_avr_S_gov_A =
        zero_coupling_state_vars_to_im_alg_vars(
         Exc, Gov)

    z_avr_S_gov_S =
        zero_coupling_state_vars_to_state_vars(
         Exc, Gov )


    z_avr_A_gov_S =
        zero_coupling_im_alg_vars_to_state_vars(
         Exc, Gov, )

    z_avr_A_gov_A =
        zero_coupling_im_alg_vars_to_im_alg_vars(
         Exc, Gov )

    #--------------------------------                

    gov_Bx_A = Gov.Bx_A(
        Gov; ω_ref0 = ωs,
        ref_in_state = ref_in_state )

    gov_Cx_A = Gov.Cx_A(
        Gov; ω_ref0 = ωs,
        ref_in_state = ref_in_state )

    #--------------------------------

    avr_Bx_A = Exc.Bx_A(
        Exc;
        ref_in_state = ref_in_state )

    avr_Cx_A = Exc.Cx_A(
        Exc;
        ref_in_state = ref_in_state )

    #--------------------------------

    gen_m = hcat(Ax_gen, Ax_gen_S_gov_S,
                 Ax_gen_S_avr_S, Ax_gen_S_gov_A,
                 Ax_gen_S_avr_A )

    gov_S_m = hcat( gov_Ax_S_Sgen, gov_Ax_S_S,
                    z_gov_S_avr_S, gov_Ax_S_A,
                    z_gov_S_avr_A )

    gov_A_m = hcat( gov_Ax_A_Sgen, gov_Ax_A_S,
                    z_gov_A_avr_S , gov_Ax_A_A,
                    z_gov_A_avr_A )

    avr_S_m = hcat( avr_Ax_S_Sgen, z_avr_S_gov_S,
                    avr_Ax_S_S, z_avr_S_gov_A,
                    avr_Ax_S_A )

    avr_A_m = hcat( avr_Ax_A_Sgen, z_avr_A_gov_S,
                    avr_Ax_A_S, z_avr_A_gov_A,
                    avr_Ax_A_A )

    #--------------------------------

    Ax = vcat(
        gen_m,
        gov_S_m,
        avr_S_m,
        gov_A_m,        
        avr_A_m)

   Bx = vcat(
       Bx_gen,
       gov_Bx_S,
       avr_Bx_S,
       gov_Bx_A,
       avr_Bx_A)

   Cx = vcat(
       Cx_gen,
       gov_Cx_S,
       avr_Cx_S,
       gov_Cx_A,
       avr_Cx_A)
    
    #--------------------------------

    return (; Ax_gen, Bx_gen, Cx_gen,
            Ax_gen_S_gov_S, Ax_gen_S_gov_A,
            Ax_gen_S_avr_S, Ax_gen_S_avr_A,
            gov_Ax_S_S, gov_Ax_S_Sgen,
            gov_Bx_S,   gov_Cx_S,
            avr_Ax_S_S, avr_Ax_S_Sgen,
            avr_Bx_S,   avr_Cx_S,
            gov_Ax_S_A, gov_Ax_A_S,
            gov_Ax_A_A, gov_Ax_A_Sgen,
            avr_Ax_S_A, avr_Ax_A_S,
            avr_Ax_A_A, avr_Ax_A_Sgen,
            z_gov_S_avr_A, z_gov_S_avr_S,
            z_gov_A_avr_S, z_gov_A_avr_A,
            z_avr_S_gov_A, z_avr_S_gov_S,
            z_avr_A_gov_S, z_avr_A_gov_A,
            gov_Bx_A,  gov_Cx_A,
            avr_Bx_A,  avr_Cx_A,
            gen_m,
            gov_S_m, gov_A_m,
            avr_S_m, avr_A_m,
            Ax, Bx, Cx)
    
end


function diagnosis_get_gen_gov_avr_im_plant_system_matrices(
    plant;
    ω_ref0 = ωs,
    ref_in_state = false )

    return diagnosis_get_gen_gov_avr_im_plant_system_matrices(
        plant.Gen,
        plant.Gov,
        plant.Exc,
        plant.dummy_vr1,
        plant.dummy_vf_tilade;
        ω_ref0 =
            ω_ref0,
        ref_in_state =
            ref_in_state )
    
end

