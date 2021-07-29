module Singularities


export LinearStrengthLineVortex,LinearStrengthLineSource,induced_velocity

abstract type Singularity end
abstract type PointSingularity <: Singularity end 
abstract type LineSingularity <: Singularity end 


struct LinearStrengthLineVortex{A,B,C,D,E,F} <: LineSingularity
    """
    A linear strength line vortex starting at (x_start,y_start)
    with strength γ_start and ending at (x_end,y_end) with strength
    γ_end. 
    """
    x_start::A
    y_start::B
    x_end::C
    y_end::D
    γ_start::E 
    γ_end::F
end 

struct LinearStrengthLineSource{A,B,C,D,E,F} <: LineSingularity
    """
    A linear strength line vortex starting at (x_start,y_start)
    with strength γ_start and ending at (x_end,y_end) with strength
    γ_end. 
    """
    x_start::A
    y_start::B
    x_end::C
    y_end::D
    σ_start::E
    σ_end::F
end

Base.Broadcast.broadcastable(a::LinearStrengthLineVortex) = Ref(a)
Base.Broadcast.broadcastable(a::LinearStrengthLineSource) = Ref(a)

singularity_length(singularity::T where T<:LineSingularity) = hypot((singularity.x_end - singularity.x_start) ,  (singularity.y_end - singularity.y_start))

macro strengths(singularity,symbol)
    ex = "($(singularity).$(symbol)"*"_start"*",$(singularity).$(symbol)"*"_end)"
    return Meta.parse(ex)
end

function get_strengths(singularity::T) where T<:LineSingularity
    symbol = fieldnames(T)[5:6]
    symbol_start,symbol_end = @strengths(singularity,symbol)
    return symbol_start, symbol_end
end 

function induced_velocity_panel_coordinates(xₚ_target,
                                            yₚ_target,
                                            strength_start,
                                            strength_end,
                                            panel_length,
                                            ln_r₂_r₁,
                                            dθ,
                                            is_on_panel,
                                            is_on_endpoint,
                                            source::LinearStrengthLineVortex) 
    dγ = strength_end - strength_start
    yₚ_target_regularized = is_on_panel ? 1 : yₚ_target

    uₚ_term1 = yₚ_target / (2π) * dγ / panel_length * ln_r₂_r₁
    uₚ_term2 = (strength_start * panel_length + dγ * xₚ_target) / (2π * panel_length) * dθ
    uₚ = is_on_panel ? 0 : uₚ_term1 + uₚ_term2

    vₚ_term1 = (strength_start * panel_length + dγ * xₚ_target) / (2π * panel_length) * ln_r₂_r₁
    vₚ_term2 = is_on_panel ? dγ / (2π) : yₚ_target / (2π) * dγ / panel_length * (panel_length / yₚ_target_regularized - dθ)
    vₚ = vₚ_term1 + vₚ_term2

    uₚ = is_on_endpoint ? 0 : uₚ
    vₚ = is_on_endpoint ? 0 : vₚ
    
    return uₚ,vₚ
end

function induced_velocity_panel_coordinates(xₚ_target,
                                            yₚ_target,
                                            strength_start,
                                            strength_end,
                                            panel_length,
                                            ln_r₂_r₁,
                                            dθ,
                                            is_on_panel,
                                            is_on_endpoint,
                                            source::LinearStrengthLineSource) 
    dσ = strength_start - strength_end
    yₚ_target_regularized = is_on_panel ? 1 : yₚ_target

    vₚ_term1 = yₚ_target / (2π) * dσ / panel_length * ln_r₂_r₁
    vₚ_term2 = (strength_start * panel_length + dσ * xₚ_target) / (2π * panel_length) * dθ
    vₚ = is_on_panel ? 0 : vₚ_term1 + vₚ_term2

    uₚ_term1 = -(strength_start * panel_length + dσ * xₚ_target) / (2π * panel_length) * ln_r₂_r₁
    uₚ_term2 = is_on_panel ? -dσ / (2π) : - yₚ_target / (2π) * dσ / panel_length * (panel_length / yₚ_target_regularized - dθ)
    uₚ = uₚ_term1 + uₚ_term2

    uₚ = is_on_endpoint ? 0 : uₚ
    vₚ = is_on_endpoint ? 0 : vₚ

    return uₚ,vₚ
end  



function induced_velocity_panel_coordinates(xₚ_target,yₚ_target,source::T where T<:LineSingularity)
    panel_length = singularity_length(source)
    is_on_panel = abs(yₚ_target) <= 1e-8
    r₁ = hypot(xₚ_target, yₚ_target)
    r₂ = hypot((xₚ_target - panel_length), yₚ_target)
    is_on_endpoint = (r₂ == 0) || (r₁ == 0)
    r₁ = r₁ == 0 ? 1 : r₁
    r₂ = r₂ == 0 ? 1 : r₂
    ln_r₂_r₁ = log(r₂/r₁)
    θ₁ = atan(yₚ_target,xₚ_target)
    θ₂ = atan(yₚ_target,xₚ_target-panel_length)
    dθ = θ₂ - θ₁
    strength_start,strength_end = get_strengths(source) 
    
    return induced_velocity_panel_coordinates(xₚ_target,
                                            yₚ_target,
                                            strength_start,
                                            strength_end,
                                            panel_length,
                                            ln_r₂_r₁,
                                            dθ,
                                            is_on_panel,
                                            is_on_endpoint,
                                            source::LinearStrengthLineVortex) 
end 


function induced_velocity(x_target,y_target,source::T where T<:LineSingularity)
    dx_panel = source.x_end - source.x_start
    dy_panel = source.y_end - source.y_start
    panel_length = hypot(dx_panel , dy_panel)
    xₚ_hat_x = dx_panel/panel_length
    xₚ_hat_y = dy_panel/panel_length
    yₚ_hat_x = - xₚ_hat_y
    yₚ_hat_y = xₚ_hat_x
    x_target_relative = x_target - source.x_start
    y_target_relative = y_target - source.y_start
    xₚ_target = x_target_relative * xₚ_hat_x + y_target_relative * xₚ_hat_y
    yₚ_target = y_target_relative * yₚ_hat_x + y_target_relative * yₚ_hat_y
    uₚ,vₚ = induced_velocity_panel_coordinates(xₚ_target,yₚ_target,source)
    u = uₚ * xₚ_hat_x + vₚ * yₚ_hat_x
    v = uₚ * xₚ_hat_y + vₚ * yₚ_hat_y
    return u,v
end


function induced_velocity(x_target,y_target,system::Array{<:Singularity})
    u,v = 0,0
    for singularity in system
        u_tmp,v_tmp = induced_velocity(x_target,y_target,singularity)
        u += u_tmp
        v += v_tmp
    end 
    return u,v
end 



end 





